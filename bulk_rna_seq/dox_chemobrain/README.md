---
title: "README"
output: html_document
---
# Chemobrain

## License
- Code (`scripts/`): MIT License, see `LICENSE`.
- Figures (`figures/`): All rights reserved.

## Structure
chemo-brain-rnaseq/
├── README.md
├── LICENSE
├── Makefile
├── .gitignore
├── renv.lock           # environment.yml / requirements.txt
│
├── data/
│   ├── README.md            # de onde vêm os dados (GSEXXXXX) e como baixar
│   └── metadata/
│       └── samples.csv      # amostra, grupo, tratamento, região do cérebro, etc.
│   
├── scripts/
│   ├── 01_download_and_prepare.R
│   ├── 02_qc_exploratory.R      # PCA, correlação entre amostras, outliers
│   ├── 03_differential_expression.R
│   └── 04_enrichment.R          # GO/KEGG/GSEA, se você fizer
│
├── results/
│   ├── tables/              # DEGs (CSV), resultados de enriquecimento
│   └── figures/             # PNG/PDF: PCA, volcano, heatmap, enriquecimento
│
└── report/                  # opcional
    │
    ├── analysis.qmd
    └── analysis.html        # versionado ou publicado no GitHub Pages