# Instalando pacotes

library(dada2)
library(phyloseq)
library(DESeq2)
library(vegan)
library(microbiome)
library(tidyverse)
library(readr)

# Defina o diretório
path <- "C:/Users/AsRock/Documents/feces"

otu_table <- read.table("otu_table_genus.txt",
                        header = TRUE, 
                        row.names = 1, 
                        sep = "\t")

metadata <- read_tsv("E-MTAB-13426.sdrf.txt")

otu_phylum <- read.delim("otu_table_phylum.txt", row.names = 1, check.names = FALSE)
otu_class <- read.delim("otu_table_class.txt", row.names = 1, check.names = FALSE)
otu_order <- read.delim("otu_table_order.txt", row.names = 1, check.names = FALSE)
otu_family <- read.delim("otu_table_family.txt", row.names = 1, check.names = FALSE)
otu_genus <- read.delim("otu_table_genus.txt", row.names = 1, check.names = FALSE)


otu_matrix <- as.matrix(otu_phylum)
otu_table_ps <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# sample_data
metadata_ps <- sample_data(metadata)

# Criar phyloseq
ps_phylum <- phyloseq(otu_table_ps, metadata_ps)
