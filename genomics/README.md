# giab-variant-pipeline

Minimal Snakemake pipeline for germline variant calling, benchmarked against the GIAB HG002 truth set.

**Steps:** fastp (QC/trimming) -> bwa mem -> samtools (sort, fixmate, markdup) -> GATK HaplotypeCaller -> bcftools filter -> rtg vcfeval.

## Setup

1. Put the GRCh38 reference in `ref/` and the paired FASTQs in `data/`.
2. Download the GIAB HG002 GRCh38 benchmark VCF (+ `.tbi`) and BED into `truth/`.
3. Edit `config.yaml` (paths, region, filter).

## Run

```bash
snakemake -n                       # dry run
snakemake --use-conda --cores 8
```

Benchmark metrics (precision, recall, F-measure for SNPs and indels) end up in `results/eval/<sample>/summary.txt`.

## Notes

- Reference contig names must match the truth set (`chr20` vs `20`).
- bcftools is used for speed and simplicity; swap in GATK HaplotypeCaller or DeepVariant for production-grade calling.
