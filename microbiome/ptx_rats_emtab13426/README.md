# Gut microbiota of paclitaxel-treated rats (E-MTAB-13426)

Snakemake + R Markdown re-analysis of 16S rRNA V3-V4 amplicon data from
Zhong et al. (2023), *Blockade of CCR5 suppresses paclitaxel-induced peripheral neuropathic
pain caused by increased deoxycholic acid*, Cell Reports 42:113386
(ArrayExpress [E-MTAB-13426](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13426),
ENA [PRJEB67467](https://www.ebi.ac.uk/ena/browser/view/PRJEB67467)).

The authors deposited only relative-abundance tables. This pipeline starts from the raw reads
and produces ASV **counts**, which allow count-based methods (Chao1, DESeq2, ANCOM-BC2,
rarefaction).

## Design

| Group | Rats | Samples | Sequencing run |
|---|---|---|---|
| Control (vehicle) | C1-C10 | 10 | AG (MiSeq 2×300) |
| PTX week 2 (day 14) | P1-P8 | 8 | AG (MiSeq 2×300) |
| PTX week 4 (day 28) | P1-P8 (same rats) | 8 | BD (MiSeq 2×250) |

**Caveat:** all week-4 samples come from a separate sequencing run, so run effects are
confounded with time. Only Control vs PTX week 2 is free of this confounding.

## Workflow

```
ENA FASTQ ──> cutadapt (both mate orientations) ──> DADA2 filterAndTrim
          ──> DADA2 learnErrors/dada per run x orientation batch ──> mergePairs ──> removeBimeraDenovo
          ──> SILVA 138.2 assignTaxonomy + addSpecies
          ──> MAFFT + FastTree (midpoint-rooted) ──> phyloseq
          ──> report.Rmd: composition, alpha/beta diversity, PERMANOVA, DESeq2 + ANCOM-BC2
```

Two quirks of the deposited data are handled explicitly:

- **Mixed mate orientation:** in 14 of the 26 runs R1 starts with the reverse primer (806R).
  `scripts/trim_primers.sh` tries both orientations and re-orients the reads, keeping 97.5-99.3%
  of pairs per sample.
- **Two read lengths:** the 2×250 run gives ~225/231 bp reads after primer removal. Truncation
  is set to 220/229 bp for all samples, which leaves ~20 bp of overlap on the ~429 bp amplicon.

## Layout

- `Snakefile`: all rules. Samples, groups and runs are parsed from the ArrayExpress SDRF.
- `config.yaml`: data paths (under `$DATA`), primers, DADA2 truncation, seed.
- `envs/`: conda environments (`cutadapt`, `dada2`, `tree`, `report`).
- `scripts/`: primer trimming and DADA2/phyloseq R scripts called by the rules.
- `report.Rmd`: statistics and figures, rendered by the `report` rule.
- `results/`, `logs/`: generated, git-ignored. Raw reads and SILVA are stored under `$DATA`.

## Running

```bash
conda activate genomic-variants       # Snakemake 9
export DATA=/mnt/sda2/HD/Datasets     # holds the ArrayExpress files; FASTQ and SILVA are downloaded here
snakemake -n --use-conda              # dry run
snakemake --use-conda --cores 8 --resources downloads=4
```

The report is written to `results/report/report.html`, and its tables to
`results/report/tables/`. The phyloseq object (`results/phyloseq/ps.rds`) and the ASV count
table (`results/phyloseq/asv_counts.tsv`) can be reused for other analyses.
