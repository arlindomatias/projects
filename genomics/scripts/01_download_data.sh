#!/usr/bin/env bash

# ============================================================
# 01_download_data.sh
# Create project directories and download FASTQ files from ENA
# Accession: ERR1341793
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

ACCESSION="ERR1341793"

# Project directory
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Raw data directory
# Uses DATA environment variable if available.
# Otherwise, saves data inside the project.

if [[ -n "${DATA:-}" ]]; then
    RAW_DATA_DIR="$DATA"
    echo "Using DATA environment variable: $RAW_DATA_DIR"
else
    RAW_DATA_DIR="${PROJECT_DIR}/data/raw"
    echo "DATA variable not defined. Using project directory."
fi

RAW_DIR="${RAW_DATA_DIR}/${ACCESSION}"

# Other project directories
TRIMMED_DIR="${PROJECT_DIR}/data/trimmed"
REFERENCE_DIR="${PROJECT_DIR}/reference"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs"

LOG_FILE="${LOG_DIR}/01_download_data.log"

# ENA API
ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"

# ------------------------------------------------------------
# Create directory structure
# ------------------------------------------------------------

echo "Creating project directories..."

mkdir -p \
    "$RAW_DIR" \
    "$TRIMMED_DIR" \
    "$REFERENCE_DIR" \
    "$RESULTS_DIR/qc" \
    "$RESULTS_DIR/alignment" \
    "$RESULTS_DIR/coverage" \
    "$RESULTS_DIR/variants" \
    "$LOG_DIR"

# ------------------------------------------------------------
# Logging
# ------------------------------------------------------------

exec > >(tee -a "$LOG_FILE") 2>&1

echo "============================================================"
echo "Starting data download"
echo "Accession: $ACCESSION"
echo "Date: $(date)"
echo "============================================================"

# ------------------------------------------------------------
# Retrieve metadata from ENA
# ------------------------------------------------------------

echo "Retrieving FASTQ metadata from ENA..."

METADATA=$(mktemp)
trap 'rm -f "$METADATA"' EXIT

curl -fsSL --retry 3 \
    "${ENA_API}?accession=${ACCESSION}&result=read_run&fields=run_accession,library_layout,fastq_ftp,fastq_md5&format=tsv" \
    -o "$METADATA"

if [[ ! -s "$METADATA" ]]; then
    echo "ERROR: Empty response from ENA."
    exit 1
fi


IFS=$'\t' read -r RUN LAYOUT FASTQ_FTP FASTQ_MD5 < <(tail -n 1 "$METADATA")

if [[ "$RUN" != "$ACCESSION" ]]; then
    echo "ERROR: Accession not found or invalid ENA response."
    exit 1
fi

if [[ "$LAYOUT" != "PAIRED" ]]; then
    echo "ERROR: Expected paired-end data. Found: $LAYOUT"
    exit 1
fi

if [[ -z "$FASTQ_FTP" || -z "$FASTQ_MD5" ]]; then
    echo "ERROR: Missing FASTQ paths or MD5 checksums."
    exit 1
fi

IFS=';' read -r -a FASTQ_FILES <<< "$FASTQ_FTP"
IFS=';' read -r -a MD5_FILES <<< "$FASTQ_MD5"

if [[ "${#FASTQ_FILES[@]}" -ne "${#MD5_FILES[@]}" ]]; then
    echo "ERROR: FASTQ and MD5 entries do not match."
    exit 1
fi

echo "Run: $RUN"
echo "Layout: $LAYOUT"
echo "FASTQ entries returned by ENA: ${#FASTQ_FILES[@]}"

# ------------------------------------------------------------
# Select paired-end FASTQ files
# ------------------------------------------------------------

SELECTED_FILES=()
SELECTED_MD5=()

R1_COUNT=0
R2_COUNT=0

for i in "${!FASTQ_FILES[@]}"; do

    FILE_NAME="${FASTQ_FILES[$i]##*/}"

    case "$FILE_NAME" in
        *_1.fastq.gz)
            SELECTED_FILES+=("${FASTQ_FILES[$i]}")
            SELECTED_MD5+=("${MD5_FILES[$i]}")
            ((R1_COUNT += 1))
            ;;
        *_2.fastq.gz)
            SELECTED_FILES+=("${FASTQ_FILES[$i]}")
            SELECTED_MD5+=("${MD5_FILES[$i]}")
            ((R2_COUNT += 1))
            ;;
        *)
            echo "Ignoring additional FASTQ representation: $FILE_NAME"
            ;;
    esac

done

if [[ "$R1_COUNT" -ne 1 || "$R2_COUNT" -ne 1 ]]; then
    echo "ERROR: Expected exactly one R1 and one R2 FASTQ file."
    exit 1
fi

echo "Selected R1 and R2 files successfully."


echo "============================================================"
echo "DOWNLOAD COMPLETED"
echo "============================================================"


for i in "${!SELECTED_FILES[@]}"; do

    FTP_PATH="${SELECTED_FILES[$i]}"
    MD5="${SELECTED_MD5[$i]}"

    FILE_NAME="${FTP_PATH##*/}"
    OUTPUT="${RAW_DIR}/${FILE_NAME}"

    echo ""
    echo "Processing: $FILE_NAME"

    # Check if file already exists
    if [[ -f "$OUTPUT" ]]; then

        echo "File already exists. Checking integrity..."

        if gzip -t "$OUTPUT" 2>/dev/null &&
           [[ -n "$MD5" ]] &&
           [[ "$(md5sum "$OUTPUT" | awk '{print $1}')" == "$MD5" ]]; then

            echo "File is valid. Skipping download."
            continue

        else

            echo "File is corrupted or MD5 mismatch."
            echo "Removing invalid file..."
            rm -f "$OUTPUT"

        fi

    fi

    # Download missing or invalid file
    echo "Downloading $FILE_NAME..."

    wget -c \
        --tries=5 \
        --timeout=60 \
        "https://${FTP_PATH}" \
        -P "$RAW_DIR"

    # Validate gzip
    echo "Checking gzip integrity..."

    gzip -t "$OUTPUT"

    # Validate MD5
    echo "Checking MD5 checksum..."

    ACTUAL_MD5=$(md5sum "$OUTPUT" | awk '{print $1}')

    if [[ "$ACTUAL_MD5" != "$MD5" ]]; then
        echo "ERROR: MD5 checksum mismatch!"
        exit 1
    fi

    echo "File successfully downloaded and validated."

done

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo "DOWNLOAD COMPLETED"
echo "============================================================"

echo "Raw data directory: $RAW_DIR"

ls -lh "$RAW_DIR"

echo ""
echo "End date: $(date)"
echo "Log: $LOG_FILE"