### Purpose: generate figures for DO paper
###          1. metagenomic PCA
###          2. lipidomic PCA
###          3. RNA PCA

# load required libraries
options(stringsAsFactors = FALSE)
source("/Users/rootqz/R/QZ_functions.R")

############ Input data ############
# load meta table
meta_tbl <- read.table("~/Desktop/DO/data/metadata_500DO.tsv", sep = "\t", header = T)

# load microbial function traits matrix
mb.dir <- "/Users/rootqz/Desktop/DO/metagenomics/2020/"
KO_abun <- readRDS(paste0(mb.dir, "data/DO.summed.TPM.KEGG.orthology.rds"))
KO_abun_rankz <- readRDS(paste0(mb.dir, "data/DO.summed.TPM.KEGG.orthology.rankZ.rds"))

taxa_abun <- readRDS(paste0(mb.dir, "data/DO.summed.TPM.ncbi.genus.100genes.rds"))

### load mbQTL peak table
qtl.KO <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_kegg_orthology.rds"))
qtl.genus <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_ncbi_genus.rds"))
qtl.MAG <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_MAGs.rds"))

# filter genus taxa (at least have > 100 genes)
qtl.genus <- qtl.genus %>%
    filter(pheno %in% names(taxa_abun))

# annotations
KEGG_gene <- read.table(paste0(mb.dir, "data/KEGG_KO2gene2EC.tsv"), sep = "\t", header = T, quote = "")
KEGG_path <- read.table(paste0(mb.dir, "data/KEGG_KO2path.tsv"), sep = "\t", header = T, quote = "")

# trim KEGG_path
KEGG_path$path_definition <- gsub("\\s*\\[[^\\]+\\]","",KEGG_path$path_definition) %>%
    trimws()
KEGG_path <- unique(KEGG_path[,c("KO", "path_definition")]) %>%
    na.omit()

# 1.9 million metagenes annotation
anno.mb.gene <- readRDS(paste0(mb.dir, "data/DO_1.9M_NRGeneSet_KEGGanno_NCBIanno.rds"))
# metagenome-assembled genome (MAG) annotation
anno.mb.mag <- read.table(paste0(mb.dir, "data/MAG_bins_metaTable.tsv"), sep = "\t")
# WGS summary
summary.wgs <- read.table(paste0(mb.dir, "data/DO_WGS_reads_meta.tsv"), sep = "\t", header = T)

### load eQTL
rna.dir <- "/Users/rootqz/Desktop/DO/RNA-seq/analysis/eQTL/"
qtl.rna <- readRDS(paste0(rna.dir, "result/qtl_blupsEff_intstine_gbrs_diploid_gene_withGeneAnno.rds"))

### load lipid QTL
qtl.lipid <- read.table("/Users/rootqz/Desktop/DO/RNA-seq/analysis/eQTL/data/qtl.summary.cecum.vl.csv", sep = ",", header = T, row.names = 1)
###################

################### Figure S1 ##########################
##### gut microbiome
# pca
pr.out <- prcomp(KO_abun, scale. = T)
pr.var=pr.out$sdev ^2
pve_mb=pr.var/sum(pr.var)
pca_matrix <- pr.out$x
pca_1and2 <- as.data.frame(pca_matrix[,1:2])
pca_1and2$mouse <- rownames(pca_1and2)

pca_merge_mb <- pca_1and2 %>%
    left_join(meta_tbl[,c("mouse", "sex", "DOgen", "DOwave")], by = c("mouse" = "mouse"))

# sex
mb_sex <- ggplot(pca_merge_mb) + 
    geom_point(aes(x = PC1, y = PC2, color = sex)) + 
    theme_bw() +
    scale_color_manual(values = c("F" = "#13547a", "M" = "#80d0c7"), name = "Sex") +
    theme(legend.title = element_text(size = 14)) +
    xlab(paste0("PC1 (", round(pve_mb[1]*100, 2), "%)")) +
    ylab(paste0("PC2 (", round(pve_mb[2]*100, 2), "%)"))

# DO wave
mb_wave <- ggplot(pca_merge_mb) + 
    geom_point(aes(x = PC1, y = PC2, color = as.factor(DOwave))) + 
    theme_bw() +
    scale_color_manual(values = gg_color_hue(20)[c(1,5,15)], name = "DO wave") +
    theme(legend.title = element_text(size = 14)) +
    xlab(paste0("PC1 (", round(pve_mb[1]*100, 2), "%)")) +
    ylab(paste0("PC2 (", round(pve_mb[2]*100, 2), "%)"))

mb.plot <- grid.arrange(mb_sex, mb_wave, nrow=1, ncol=2)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/PCA_mb_function_trait.eps", height = 5, width = 12)
grid.arrange(mb_sex, mb_wave, nrow=1, ncol=2)
dev.off()


##### cecum lipid
# load cecum lipid matrix
lipid_abun <- readRDS("~/Desktop/DO/metabolomics/data_from_CoonLab/normalized/attie_cecum_lipids_normalized.rds")
lipid_abun_matrix <- lipid_abun[,-c(1:11)]

# pca
pr.out <- prcomp(lipid_abun_matrix, scale. = T)
pr.var <- pr.out$sdev ^2
pve_lipid <- pr.var/sum(pr.var)
pca_matrix <- pr.out$x
pca_1and2 <- as.data.frame(pca_matrix[,1:2])
pca_1and2$mouse <- rownames(pca_1and2)

pca_merge_lipid <- pca_1and2 %>%
    left_join(meta_tbl[,c("mouse", "sex", "DOgen", "DOwave")], by = c("mouse" = "mouse"))

# sex
lipid_sex <- ggplot(pca_merge_lipid) + 
    geom_point(aes(x = PC1, y = PC2, color = sex)) + 
    theme_bw() +
    scale_color_manual(values = c("F" = "#13547a", "M" = "#80d0c7"), name = "Sex") +
    theme(legend.title = element_text(size = 14)) +
    xlab(paste0("PC1 (", round(pve_lipid[1]*100, 2), "%)")) +
    ylab(paste0("PC2 (", round(pve_lipid[2]*100, 2), "%)"))

# DO wave
lipid_wave <- ggplot(pca_merge_lipid) + 
    geom_point(aes(x = PC1, y = PC2, color = as.factor(DOwave))) + 
    theme_bw() +
    scale_color_manual(values = gg_color_hue(20)[c(1,5,10,15)], name = "DO wave") +
    theme(legend.title = element_text(size = 14)) +
    xlab(paste0("PC1 (", round(pve_lipid[1]*100, 2), "%)")) +
    ylab(paste0("PC2 (", round(pve_lipid[2]*100, 2), "%)"))

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/PCA_lipid_trait.eps", height = 5, width = 12)
grid.arrange(lipid_sex, lipid_wave, nrow=1, ncol=2)
dev.off()

##### intestine transcript
### load DO small intestine RNA-Seq counts
rna.data.dir <- "/Users/rootqz/Desktop/DO/RNA-seq/analysis/eQTL/"
gene_tpm <- readRDS(paste0(rna.data.dir, "data/gbrs.gene.tpm.filter.rds"))

# pca
pr.out <- prcomp(gene_tpm, scale. = T)
pr.var <- pr.out$sdev ^2
pve_gene <-pr.var/sum(pr.var)
pca_matrix <- pr.out$x
pca_1and2 <- as.data.frame(pca_matrix[,1:2])
pca_1and2$mouse <- rownames(pca_1and2)

pca_merge_gene <- pca_1and2 %>%
    left_join(meta_tbl[,c("mouse", "sex", "DOgen", "DOwave")], by = c("mouse" = "mouse"))

# sex
gene_sex <- ggplot(pca_merge_gene) + 
    geom_point(aes(x = PC1, y = PC2, color = sex)) + 
    theme_bw() +
    scale_color_manual(values = c("F" = "#13547a", "M" = "#80d0c7"), name = "Sex") +
    theme(legend.title = element_text(size = 14)) +
    xlab(paste0("PC1 (", round(pve_gene[1]*100, 2), "%)")) +
    ylab(paste0("PC2 (", round(pve_gene[2]*100, 2), "%)"))

# DO wave
gene_wave <- ggplot(pca_merge_gene) + 
    geom_point(aes(x = PC1, y = PC2, color = as.factor(DOwave))) + 
    theme_bw() +
    scale_color_manual(values = gg_color_hue(20)[c(1,5,15,20)], name = "DO wave") +
    theme(legend.title = element_text(size = 14)) +
    xlab(paste0("PC1 (", round(pve_gene[1]*100, 2), "%)")) +
    ylab(paste0("PC2 (", round(pve_gene[2]*100, 2), "%)"))

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/PCA_intestine_transcript_trait.eps", height = 5, width = 12)
grid.arrange(gene_sex, gene_wave, nrow=1, ncol=2)
dev.off()
