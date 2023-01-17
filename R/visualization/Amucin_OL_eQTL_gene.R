### Purpose: plot correlation between co-mapped eQTL gene w/ Amucin MAG
### Created: 2021-08-17

library(qtl2)

options(stringsAsFactors = FALSE)
source("/Users/rootqz/R/QZ_functions.R")

phe_MAG <- readRDS("~/Desktop/DO/metagenomics/2020/data/DO.MAGs.GenomeCov.libSizeNorm.Amuciniphila.rds")
phe_taxa <- readRDS("~/Desktop/DO/metagenomics/2020/data/DO.summed.TPM.taxa.rds")
phe_rna <- readRDS("~/Desktop/DO/RNA-seq/analysis/eQTL/data/gbrs.gene.tpm.filter.rds")

phe_taxa_cor <- phe_taxa[which(rownames(phe_taxa) %in% rownames(phe_rna)),]
phe_rna_cor <- phe_rna[which(rownames(phe_rna) %in% rownames(phe_taxa)),]

# Atf3, ENSMUSG00000026628
phe_taxa_cor %>%
    cbind(phe_rna_cor) %>%
    ggplot(aes(x = log10(g_Akkermansia), y = ENSMUSG00000026628)) +
    geom_point() +
    theme_classic() +
    geom_smooth(method='lm')

# Tifa, ENSMUSG00000046688
phe_taxa_cor %>%
    cbind(phe_rna_cor) %>%
    ggplot(aes(x = log10(g_Akkermansia), y = ENSMUSG00000046688)) +
    geom_point() +
    theme_classic() +
    geom_smooth(method='lm')

# Jmjd8, ENSMUSG00000025736
phe_taxa_cor %>%
    cbind(phe_rna_cor) %>%
    ggplot(aes(x = log10(g_Akkermansia), y = ENSMUSG00000025736)) +
    geom_point() +
    theme_classic() +
    geom_smooth(method='lm')

phe_rna_cor %>%
    ggplot(aes(x = ENSMUSG00000026628, y = ENSMUSG00000046688)) +
    geom_point() +
    theme_classic() +
    geom_smooth(method='lm')
