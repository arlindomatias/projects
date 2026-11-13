#!/usr/bin/env bash
# Remove primers from one paired sample whose read orientation is unknown.
# In E-MTAB-13426 some runs have R1 starting with the forward primer (standard) and others
# with the reverse primer (swapped). Pass A takes swapped pairs, pass B takes standard pairs
# from what pass A left untrimmed; both primers must be found (--pair-filter=any).
# Outputs are oriented so that *_fwd always starts after the forward primer.
# Usage: trim_primers.sh R1 R2 FWD_PRIMER REV_PRIMER OUT_FWD OUT_REV OUT_JSON THREADS
set -euo pipefail
r1=$1; r2=$2; pf=$3; pr=$4; out_fwd=$5; out_rev=$6; out_json=$7; threads=$8

tmp=$(mktemp -d "$(dirname "$out_fwd")/tmp.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
opts=(--minimum-length 100 --pair-filter=any -j "$threads")

# Pass A: R1 starts with the reverse primer -> R1 is the "rev" read.
cutadapt -g "$pr" -G "$pf" "${opts[@]}" --json "$tmp/a.json" \
  --untrimmed-output "$tmp/u1.fq.gz" --untrimmed-paired-output "$tmp/u2.fq.gz" \
  -o "$tmp/a_rev.fq.gz" -p "$tmp/a_fwd.fq.gz" "$r1" "$r2"

# Pass B: R1 starts with the forward primer -> R1 is the "fwd" read.
cutadapt -g "$pf" -G "$pr" "${opts[@]}" --json "$tmp/b.json" --discard-untrimmed \
  -o "$tmp/b_fwd.fq.gz" -p "$tmp/b_rev.fq.gz" "$tmp/u1.fq.gz" "$tmp/u2.fq.gz"

cat "$tmp/a_fwd.fq.gz" "$tmp/b_fwd.fq.gz" > "$out_fwd"
cat "$tmp/a_rev.fq.gz" "$tmp/b_rev.fq.gz" > "$out_rev"

python3 - "$tmp/a.json" "$tmp/b.json" "$out_json" <<'PY'
import json, sys
a, b = (json.load(open(f))["read_counts"] for f in sys.argv[1:3])
swapped, standard = a["output"], b["output"]
summary = {
    "read_counts": {"input": a["input"], "output": swapped + standard},
    "swapped_pairs": swapped,
    "standard_pairs": standard,
    "orientation": "swapped" if swapped > standard else "standard",
}
json.dump(summary, open(sys.argv[3], "w"), indent=2)
print(summary)
PY
