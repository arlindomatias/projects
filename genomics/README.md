# giab-variant-pipeline

Minimal Snakemake pipeline for germline variant calling, benchmarked against the GIAB HG002 truth set.

**Steps:** fastp (QC/trimming) -> bwa mem -> samtools (sort, fixmate, markdup) -> GATK HaplotypeCaller -> bcftools filter -> rtg vcfeval.

**Extras:** SnpEff annotation of the filtered VCF, and DeepVariant (Docker) on the same BAM, benchmarked the same way.

## Setup

1. Put the GRCh38 reference in `ref/` and the paired FASTQs under `$DATA` (see `config.yaml`).
2. Download the GIAB HG002 GRCh38 benchmark VCF (+ `.tbi`) and BED into `truth/`.
3. Edit `config.yaml` (paths, region, filter).

## Run

```bash
snakemake -n                       # dry run
snakemake --use-conda --cores 8
```

Requires the `DATA` environment variable and Docker (for DeepVariant).

Outputs:
- GATK benchmark: `results/eval/<sample>/summary.txt`
- DeepVariant benchmark: `results/eval_deepvariant/<sample>/summary.txt`
- SnpEff annotation: `results/annot/<sample>.snpeff.vcf.gz` (+ `.html` report)
- Logs: `logs/<rule>/<sample>.log`

## Notes

- Reference contig names must match the truth set (`chr20` vs `20`).
- The SnpEff database (`snpeff_db`) is downloaded once into `snpeff_data_dir`.
