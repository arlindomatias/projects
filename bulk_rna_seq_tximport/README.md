# Bulk RNA-seq from Salmon/tximport: differential expression and enrichment

A Snakemake generalization of the [bulk-rna-seq](https://github.com/arlindomatias/bulk-rna-seq) pipeline
(GSE178812). The original pipeline reads a gene-symbol count matrix from GEO; this one starts from
**Salmon quantifications summarized with tximport** and is configured per dataset in `config.yaml`.

Current dataset: **GSE295281**: hippocampus of 12-week-old female C57BL/6J mice, one month after AC-T
chemotherapy (doxorubicin + cyclophosphamide, then paclitaxel) vs saline (n = 5 vs 6). The input comes from
[`../salmon`](../salmon) (Salmon 2.11.0, GENCODE M35).

## What changed relative to bulk-rna-seq

| | bulk-rna-seq | this project |
|---|---|---|
| Input | GEO xlsx of counts by gene symbol | `txi_gene.rds` (tximport) + `samples.tsv` + `tx2gene.tsv` |
| DESeq2 object | `DESeqDataSetFromMatrix` | `DESeqDataSetFromTximport` (uses transcript-length offsets) |
| Gene IDs | symbols, duplicates dropped | Ensembl gene IDs (version stripped), symbol kept as annotation |
| GSEA keys | `SYMBOL` | `ENSEMBL` for GO, Entrez for KEGG (ties resolved by largest \|stat\|) |
| Groups, filter, cutoffs | hard-coded | `config.yaml` (filter default: smallest group size) |
| KEGG | `gseKEGG`, downloads KEGG at every run | KEGG REST downloaded once to `cache/kegg/` (rule `kegg_annotation`), then `GSEA(TERM2GENE)` |
| Pathview | fixed pathway | configured list, or top N significant KEGG pathways |
| Runner | Snakemake calling `Rscript` | Snakemake `script:` + per-rule conda env and logs |
| Extra QC | | sample distance heatmap, labelled PCA, marker-gene plot |

## Usage

```bash
conda activate genomic-variants        # Snakemake 9
snakemake -n --use-conda               # dry-run
snakemake --use-conda --cores 2
```

The R environment is built from `envs/deseq2.linux-64.pin.txt` (exact package versions, copied from the
bulk-rna-seq lock file); `envs/deseq2.yaml` is the unpinned fallback. Step 03 needs internet access (STRING, pathview).
The KEGG REST API can be very slow (about 10 min for the 1 MB mouse gene-pathway list on 2026-11-12), so the
annotation is cached in `cache/kegg/`; delete that folder to refresh it.

### Another dataset

Copy `config.yaml`, point `input:` to the tximport object and sample table, set the condition column and
group labels (syntactic R names, e.g. `ACT` instead of `AC-T`) and run with
`snakemake --use-conda --cores 2 --configfile other.yaml`. Outputs go to `results/` and `plots/`, so use
a separate checkout or working directory per dataset.

## Pipeline

| Rule | Script | Main outputs |
|---|---|---|
| `prepare` | `01_prepare.R` | `results/dds_raw.rds`, `results/colData.csv` |
| `differential` | `02_differential_expression_and_qc.R` | `results/top_table.csv`, `results/de_summary.tsv`, `results/normalized_counts.csv`, QC/DE plots |
| `kegg_annotation` | (curl) | `cache/kegg/{org}_link_pathway.tsv`, `cache/kegg/{org}_list_pathway.tsv` |
| `enrichment` | `03_enrichment.R` | `results/GSEA_GO_all_full.csv`, `results/GSEA_KEGG_full.csv`, `results/enrichment_STRING_full.csv`, `results/{edges,nodes}_PPI.csv`, enrichment plots |

`00_common.R` is sourced by every step (logging, config, colors). Logs are written to `logs/{rule}.log`.

## Methods

- **Filtering:** genes with ≥ 10 counts in at least as many samples as the smallest group (5).
- **Differential expression:** DESeq2, design `~ group` with the reference group from `config.yaml`,
  `fitType = "local"`. DEGs: `padj < 0.05` and `|log2FC| > 0.5`. `top_table.csv` also has the apeglm-shrunken
  log2FC (used in the MA plot).
- **GSEA:** genes ranked by the Wald statistic; GO (BP, CC, MF) and KEGG with `clusterProfiler`,
  BH adjustment, fixed seed. Positive NES = enriched in the treatment group.
- **STRING:** v12.0, mouse, combined score ≥ 700; one gene per symbol (lowest padj); background = all tested genes.

## Notes for GSE295281

- The first digit of the mouse ID (2/3 = saline, 5/6 = AC-T) matches the `ACTGolgi_N` suffix in the GEO count
  matrix: probably cage or processing batch. It is **nested within the group**, so it cannot be added as a
  covariate.
- SRR33247678 (saline, mouse 300) has almost no *Xist*, although all mice are female and no Y-linked gene is
  expressed. Check `plots/marker_genes.png`, `plots/pca_vst.png` and `plots/sample_distances.png`; it can be
  removed with `exclude_samples` in a sensitivity analysis.

## Data

GEO [GSE295281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295281) (PMID 41519984). Please cite the
original authors.
