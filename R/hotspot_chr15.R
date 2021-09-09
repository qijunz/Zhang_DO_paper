### Purpose: generate figures for DO paper
###          1. 
###          2. 
###          3. 

# load required libraries
library(dplyr)
library(ggplot2)
library(pheatmap)

options(stringsAsFactors = FALSE)
source("/Users/rootqz/R/QZ_functions.R")

mb.dir <- "/Users/rootqz/Desktop/DO/metagenomics/2020/"

# Load QTL table
qtl.KO <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_kegg_orthology.rds"))
qtl.taxa <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD5_taxa.rds"))

# load KEGG annotation
KEGG_gene <- read.table(paste0(mb.dir, "data/KEGG_KO2gene2EC.tsv"), sep = "\t", header = T, quote = "")
KEGG_path <- read.table(paste0(mb.dir, "data/KEGG_KO2path.tsv"), sep = "\t", header = T, quote = "")

# trim KEGG_path
KEGG_path$path_definition <- gsub("\\s*\\[[^\\]+\\]","",KEGG_path$path_definition) %>%
    trimws()
KEGG_path <- unique(KEGG_path[,c("KO", "path_definition")]) %>%
    na.omit()

### function to plot allele effect
#
# input:
#   ch:    the position of chromosome
#   start: start position of interested region
#   end:   end postion of interested region
# 
# return the tbl with phenotype name, A-H allele effects, heatmap label, phenotype type
#
qtl_selection <- function(ch, start, end) {
    meta.ko <- qtl.KO %>%
        mutate(label = paste0(pheno, ": ", gene_name, ", ", gene_definition)) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "KO") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    meta.taxa <- qtl.taxa %>%
        mutate(label = pheno) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "taxa") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    return(rbind(meta.ko, meta.taxa))
}

# function of allele effects heatmap of interested QTL region 
qtl_heatmap <- function(qtl_tbl, ch, start, end) {
    plot.df <- qtl_tbl %>%
        column_to_rownames("label") %>%
        select(A, B, C, D, E, F, G, H)
    
    anno.row <- qtl_tbl[,"group",drop=F]
    names(anno.row) <- "Group"
    rownames(anno.row) <- rownames(plot.df)
    
    p <- pheatmap(plot.df,
                  cluster_cols = F,
                  color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(400), 
                  breaks = seq(-2, 2, by = 0.01),
                  border_color = NA,
                  annotation_row = anno.row,
                  scale = "row",
                  clustering_method = "ward.D2",
                  annotation_colors = list(Group = c(KO = "grey40", 
                                                     taxa = "grey70"))[1],
                  cellwidth = 20,
                  cellheight = 5,
                  fontsize_row = 5,
                  main = paste0("Allele Effects: Chr", ch, ": ", start, "-", end, "Mbp"))
    
    cuttree <- length(table(cutree(p$tree_row,h = 5)))
    
    p_curtree <- pheatmap(plot.df,
                          cluster_cols = F,
                          color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(400), 
                          breaks = seq(-2, 2, by = 0.01),
                          border_color = NA,
                          annotation_row = anno.row,
                          scale = "row",
                          clustering_method = "ward.D2",
                          annotation_colors = list(Group = c(KO = "grey40", 
                                                             taxa = "grey70"))[1],
                          cellwidth = 20,
                          cellheight = 5,
                          fontsize_row = 5,
                          cutree_rows = cuttree, 
                          main = paste0("Allele Effects: Chr", ch, ": ", start, "-", end, "Mbp"))
    
    
    return(p_curtree)
}

# chr15
chr <- "15"
start <- 61
end <- 65
qtl_plot <- qtl_selection(chr,start,end)
qtl_plot$cor <- NA
for (i in 1:nrow(qtl_plot)) {
    this_cor <- cor.test(as.numeric(qtl_plot[i,5:12]), as.numeric(qtl_plot[143,5:12]), method = "pearson")
    qtl_plot$cor[i] <- this_cor$estimate
}

plot_out <- qtl_plot %>%
    filter(cor > 0.7 | cor < -0.7) %>%
    qtl_heatmap(chr,start,end)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/hotspot_chr15.eps", height = 20, width = 9)
plot_out
dev.off()

# chr15, make it readable by higher cutoff
chr <- "15"
start <- 61
end <- 65
qtl_plot <- qtl_selection(chr,start,end)
qtl_plot$cor <- NA
for (i in 1:nrow(qtl_plot)) {
    this_cor <- cor.test(as.numeric(qtl_plot[i,5:12]), as.numeric(qtl_plot[143,5:12]), method = "pearson")
    qtl_plot$cor[i] <- this_cor$estimate
}
qtl_plot <- qtl_plot[which(qtl_plot$group == "taxa" | qtl_plot$group == "KO" & qtl_plot$lod > 7.2),]
plot_out <- qtl_plot %>%
    filter(cor > 0.7 | cor < -0.7) %>%
    qtl_heatmap(chr,start,end)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/hotspot_chr15_cut.eps", height = 13, width = 9)
plot_out
dev.off()

### 16s in founders
taxa <- read.table("~/Desktop/DO/omics/integration/16S/founder/taxa_phylum.txt", sep = "\t", header = T, check.names = F)
taxa <- taxa %>%
    t() %>%
    as.data.frame()
names(taxa) <- taxa[1,]
taxa <- taxa[-1,,drop=F] %>% data.matrix() %>% as.data.frame()

taxa$mouse <- rownames(taxa)
taxa$strain <- NA

# rename
taxa$strain[which(startsWith(taxa$mouse, "A"))] <- "AJ"
taxa$strain[which(startsWith(taxa$mouse, "B6"))] <- "B6"
taxa$strain[which(startsWith(taxa$mouse, "129"))] <- "129"
taxa$strain[which(startsWith(taxa$mouse, "NOD"))] <- "NOD"
taxa$strain[which(startsWith(taxa$mouse, "NZO"))] <- "NZO"
taxa$strain[which(startsWith(taxa$mouse, "Cast"))] <- "CAST"
taxa$strain[which(startsWith(taxa$mouse, "Pwk"))] <- "PWK"
taxa$strain[which(startsWith(taxa$mouse, "WSB"))] <- "WSB"

taxa$strain <- factor(taxa$strain, levels = names(CCcolors))

p1 <- ggplot(taxa) +
    geom_boxplot(aes(x = strain, y = `k__Bacteria;p__Bacteroidetes`, fill = strain, color = strain)) +
    theme_classic() +
    scale_fill_manual(values = CCcolors) +
    scale_color_manual(values = CCcolors) +
    xlab("DO founder strains") +
    theme(axis.title = element_text(size = 12))

p2 <- ggplot(taxa) +
    geom_boxplot(aes(x = strain, y = `k__Bacteria;p__Firmicutes`, fill = strain, color = strain)) +
    theme_classic() +
    scale_fill_manual(values = CCcolors) +
    scale_color_manual(values = CCcolors) +
    xlab("DO founder strains") +
    theme(axis.title = element_text(size = 12))

taxa$B2F <- taxa$`k__Bacteria;p__Bacteroidetes`/taxa$`k__Bacteria;p__Firmicutes`
p3 <- ggplot(taxa) +
    geom_boxplot(aes(x = strain, y = B2F, fill = strain, color = strain)) +
    theme_classic() +
    scale_fill_manual(values = CCcolors) +
    scale_color_manual(values = CCcolors) +
    xlab("DO founder strains") +
    ylab("Bacteroidetes/Firmicutes ratio") +
    theme(axis.title = element_text(size = 12))

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/16s_taxa_FB.eps", height = 9, width = 9)
grid.arrange(p1,p2,p3,  nrow=3, ncol=1)
dev.off()

### spore genes in MAGs
mag.anno <- read.table(paste0(mb.dir, "data/MAG_bins_metaTable.tsv"), sep = "\t")
ko_profile <- read.table("omics/integration/binning/result/KEGG_profiles_present_in_HQbins.tsv", sep = "\t", header = T)

KEGG_gene$label <- paste0(KEGG_gene$KO, ": ", KEGG_gene$gene_definition)

ko_spore <- qtl_plot$pheno[c(grep("spore", qtl_plot$label), grep("sporulation", qtl_plot$label))]
ko_spore_all <- KEGG_gene$KO[c(grep("spore", KEGG_gene$gene_definition), grep("sporulation", KEGG_gene$gene_definition))]

ko_spore_label <- KEGG_gene[which(KEGG_gene$KO %in% ko_spore), "label"]

anno_mag <- mag.anno %>% filter(bin %in% rownames(ko_profile))
rownames(anno_mag) <- anno_mag$bin
anno_mag <- anno_mag[,c("phylum"),drop=F]
anno_mag <- anno_mag[order(anno_mag$phylum),,drop=F]

## spores KO
plot_negative <- ko_profile[rownames(anno_mag),which(names(ko_profile) %in% ko_spore)]
names(plot_negative) <- ko_spore_label

## LPS KO
ko_chr15 <- qtl_plot %>% filter(cor < -0.9, group == "KO")
# ko_chr15 <- qtl_plot %>% filter(cor > 0.9, group == "KO")
ko_chr15 <- ko_chr15$pheno
ko_chr15 <- ko_chr15[which(ko_chr15 %in% names(ko_profile))]

plot_positive <- ko_profile[rownames(anno_mag), which(names(ko_profile) %in% ko_chr15)]
names(plot_positive) <- KEGG_gene$label[match(ko_chr15, KEGG_gene$KO)]

plot_combine <- cbind(plot_negative, plot_positive)

anno_ko <- as.data.frame(matrix(nrow = ncol(plot_combine), ncol = 1))
names(anno_ko) <- "group"
rownames(anno_ko) <- c(names(plot_negative), names(plot_positive))
anno_ko$group[which(rownames(anno_ko) %in% names(plot_positive))] <- "positive"
anno_ko$group[which(rownames(anno_ko) %in% names(plot_negative))] <- "negative"

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/MAG_KO_chr15.eps", height = 9, width = 10)
pheatmap(plot_combine,
         show_rownames = F,
         annotation_row = anno_mag,
         annotation_col = anno_ko,
         clustering_method = "ward.D2",
         annotation_colors = list(phylum = c(Actinobacteria = pal_rickandmorty()(5)[1], 
                                             Bacteroidetes = pal_rickandmorty()(5)[5], 
                                             Firmicutes = pal_rickandmorty()(5)[3], 
                                             Tenericutes = pal_rickandmorty()(5)[4], 
                                             Verrucomicrobia = pal_rickandmorty()(5)[2],
                                             unclassified = "grey70"),
                                  group = c(positive = "#f5a962",
                                            negative = "#39a6a3")),
         color = c("black", "yellow"),
         breaks = c(0, 0.5, 1))
dev.off()

# enrichment, fisher's exact test
ko_chr15_posi <- qtl_plot %>% filter(cor < -0.9, group == "KO")
ko_chr15_posi <- ko_chr15_posi$pheno
fisher_chr15_posi <- fisher_exact_test(names(ko_profile), ko_chr15_posi)

ko_chr15_nega <- qtl_plot %>% filter(cor > 0.9, group == "KO")
ko_chr15_nega <- ko_chr15_nega$pheno
fisher_chr15_nega <- fisher_exact_test(names(ko_profile), ko_chr15_nega)

### allele effects plot

blups_F <- qtl.taxa %>%
    filter(pheno == "p_Firmicutes", chr == "15") %>%
    allele_eff_plot()

blups_B <- qtl.taxa %>%
    filter(pheno == "p_Bacteroidetes", chr == "15") %>%
    allele_eff_plot()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/allele_effects_taxa_Firmicutes.eps", height = 2.5, width = 6)
blups_F
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/allele_effects_taxa_Bacteroidetes.eps", height = 2.5, width = 6)
blups_B
dev.off()

## snps association
# load dataset containing all QTL input files:
load(file = "/Users/rootqz/Desktop/DO/data/DO_QTL_allInput_without_Phe_df_20181205.RData")
markers <- readRDS("/Users/rootqz/Desktop/DO/data/DO_69k_markers.rds")

# WHERE:
# pheno_clin        - clinical phenotypes (500 x 170)
# pheno_clin_dict   - dictionary of clinical phenotypes
# K                 - list of "loco" kinship matrices (500 x 500)
# pmap              - physical map for the 69k grid
# probs             - genotype probabilities
# markers           - data frame for 69k grid with marker, chr, pos, cM, bp

KO_abun <- readRDS(paste0(mb.dir, "data/DO.summed.TPM.KEGG.orthology.rds"))
KO_abun_rankz <- readRDS(paste0(mb.dir, "data/DO.summed.TPM.KEGG.orthology.rankZ.rds"))

# load snp and gene 
# SNP database
snp_db <- "/Users/rootqz/Desktop/DO/data/cc_variants.sqlite"
# mouse genes database
gene_db <- "/Users/rootqz/Desktop/DO/data/mouse_genes_mgi.sqlite"

# set covar
covar <- cbind(pheno_clin[,"mouse",drop=F],
               sex = (pheno_clin$sex == "M")*1,
               dietday = pheno_clin$diet_days,
               DOwave2 = (pheno_clin$DOwave == 2)*1,
               DOwave4 = (pheno_clin$DOwave == 4)*1)

covar$mouse <- NULL

# use K06385 as representative of trait
# use PC1 of KO w/ LOD > 7.8 at cht15:63 hotspot as representative of trait
chr <- "15"
kos <- qtl.KO %>% filter(chr == "15", lod > 7.8, peak_mbp > 62, peak_mbp < 65)
kos <- kos$pheno

pca_df <- KO_abun[,kos]
pca_out <- as.data.frame(prcomp(pca_df)$x)
phe_pca <- pca_out[,1,drop=F]
phe_pca$PC1 <- rankZ(phe_pca$PC1)
phe <- phe_pca

# QTL scan1 from list:
qtl_scan1 <- scan1(probs, phe, K, covar, cores = 0)

# query snps
query_variants <- create_variant_query_func(snp_db)
# find SNPs, keeping ALL OF THEM:
out_snps <- scan1snps(probs[,chr], pmap, phe, 
                      K[chr], addcovar=covar,
                      query_func=query_variants, chr=chr, start=61, end=65,
                      keep_all_snps=TRUE)

# top SNPs
topSNPs <- top_snps(out_snps$lod, out_snps$snpinfo, drop = 1.5, show_all_snps = FALSE)

# snp association
snpinfo <- index_snps(pmap, out_snps$snpinfo)                # create index
snpprobs <- genoprob_to_snpprob(probs, snpinfo)              # convert to SNP probs
scan_snppr <- scan1(snpprobs, phe, kinship = K[[chr]], addcovar = covar)

# query gene
query_genes <- create_gene_query_func(gene_db)
genes <- query_genes(chr = chr, 61, 65)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/KO_QTL_chr15_hs_snps.eps", height = 4, width = 10)
plot(out_snps$lod, out_snps$snpinfo, col= "grey30", xlab = "", ylim = c(1,9))
dev.off()

### 2021-07-29, allele effects heatmap for only taxa trait
chr <- "15"
start <- 61
end <- 65
qtl_plot <- qtl_selection(chr,start,end)
qtl_plot$cor <- NA
for (i in 1:nrow(qtl_plot)) {
    this_cor <- cor.test(as.numeric(qtl_plot[i,5:12]), as.numeric(qtl_plot[143,5:12]), method = "pearson")
    qtl_plot$cor[i] <- this_cor$estimate
}

plot_out <- qtl_plot %>%
    filter(cor > 0.7 | cor < -0.7) %>%
    filter(group == "taxa") %>%
    qtl_heatmap(chr,start,end)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/hotspot_chr15_taxa.eps", height = 4, width = 6)
plot_out
dev.off()


