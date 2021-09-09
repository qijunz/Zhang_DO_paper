### Purpose: generate figures for DO paper
###          1. comapping exams for lipids
###          2. comapping exams for hotsopts
###          3. 


# load required libraries
library(dplyr)
library(tibble)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

options(stringsAsFactors = FALSE)

# set path
mb.dir <- "/Users/rootqz/Desktop/DO/metagenomics/2020/"

### Load QTL table
qtl.KO <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_kegg_orthology.rds"))
qtl.taxa <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD5_taxa.rds"))
qtl.MAG <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_MAGs.rds"))
qtl.MAG.akk <- readRDS("~/Desktop/DO/omics/integration/project/QTL_MAG/A_mucin_221mouse/qtl_blup_tbl_minLOD5_withMarker_221mice.rds")

qtl.lipid <- readRDS("~/Desktop/DO/RNA-seq/analysis/eQTL/data/qtl.summary.cecum.vl.rds")
qtl.lipid.ol <- readRDS("~/Desktop/DO/metabolomics/OLs/cecum_ols_221mice_2020-07-16/qtl_blup_tbl_minLOD5_withMarker_221mice.rds")
qtl.lipid.ol$pheno <- gsub("UNK", "OL", qtl.lipid.ol$pheno)

qtl.rna <- readRDS("~/Desktop/DO/RNA-seq/analysis/eQTL/result/qtl_blupsEff_intstine_gbrs_diploid_gene_withGeneAnno.rds")

# load KEGG annotation
KEGG_gene <- read.table(paste0(mb.dir, "data/KEGG_KO2gene2EC.tsv"), sep = "\t", header = T, quote = "")
KEGG_path <- read.table(paste0(mb.dir, "data/KEGG_KO2path.tsv"), sep = "\t", header = T, quote = "")

# trim KEGG_path
KEGG_path$path_definition <- gsub("\\s*\\[[^\\]+\\]","",KEGG_path$path_definition) %>%
    trimws()
KEGG_path <- unique(KEGG_path[,c("KO", "path_definition")]) %>%
    na.omit()

# function to plot allele effect
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
    
    lipid <- qtl.lipid %>%
        mutate(label = identifier) %>%
        mutate(pheno = identifier, chr = qtl.chr, peak_mbp = qtl.pos, lod = qtl.lod) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "lipid") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    rna <- qtl.rna %>%
        mutate(label = paste0("eQTL: ", gene_name)) %>%
        mutate(chr = qtl_chr) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "eQTL") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    return(rbind(meta.ko, lipid, rna))
}

# function of allele effects heatmap of interested QTL region 
qtl_heatmap <- function(qtl_tbl, ch, start, end) {
    plot.df <- qtl_tbl %>%
        remove_rownames() %>%
        column_to_rownames(var="label") %>%
        select(A, B, C, D, E, F, G, H)
    
    names(plot.df) <- c("AJ", "B6", "129", "NOD",
                        "NZO", "CAST", "PWK", "WSB")
    
    anno.row <- qtl_tbl[,"group",drop=F]
    rownames(anno.row) <- rownames(plot.df)
    
    p <- pheatmap(plot.df,
                  cluster_cols = F,
                  color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(400), 
                  breaks = seq(-2, 2, by = 0.01),
                  border_color = NA,
                  scale = "row",
                  clustering_method = "ward.D2",
                  annotation_row = anno.row,
                  annotation_colors = list(group = c(lipid = "#63C29C", 
                                                     eQTL="#FFCB04", 
                                                     KO = "grey40", 
                                                     taxa = "grey70"))[1],
                  main = paste0("Allele Effects: Chr", ch, ": ", start, "-", end, "Mbp"))
    
    cuttree <- length(table(cutree(p$tree_row,h = 5)))
    
    p_curtree <- pheatmap(plot.df,
                          cluster_cols = F,
                          color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(400), 
                          breaks = seq(-2, 2, by = 0.01),
                          border_color = NA,
                          scale = "row",
                          clustering_method = "ward.D2",
                          annotation_row = anno.row,
                          annotation_colors = list(group = c(lipid = "#63C29C", 
                                                             eQTL="#FFCB04", 
                                                             KO = "grey40", 
                                                             taxa = "grey70", 
                                                             MAG="grey90"))[1],
                          cellwidth = 25,
                          cellheight = 10,
                          fontsize_row = 9,
                          cutree_rows = cuttree, 
                          main = paste0("Allele Effects: Chr", ch, ": ", start, "-", end, "Mbp"))
    
    return(p_curtree)
}

# chr4 hs
f5a <- qtl_selection("4",48,52)
f5a.filter <- f5a[c(1,5,19,21),]
f5a.fig <- qtl_heatmap(f5a.filter,"4",48,52)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/comapping_chr4_50mbp.eps", height = 4, width = 7)
f5a.fig
dev.off()

# chr17 hs
f5d <- qtl_selection("17",30,34)
f5d.filter <- f5d[c(21:43,57,58),]
f5d.fig <- qtl_heatmap(f5d.filter,"17",30,34)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/comapping_chr17_32mbp.eps", height = 8, width = 7)
f5d.fig
dev.off()


# chr15
hs_chr15 <- qtl_selection("15",61,65)
hs_chr15_filter <- hs_chr15[which(hs_chr15$group != "KO"),]
hs_chr15.fig <- qtl_heatmap(hs_chr15,"15",61,65)

# chr8
hs_chr8 <- qtl_selection("8",10.5,14.5)
hs_chr8_filter <- hs_chr8[c(which(hs_chr8$group == "KO"), 30,31,41),]
hs_chr8.fig <- qtl_heatmap(hs_chr8_filter,"8",10.5,14.5)

### 2021-06-15, Akk/OL co-mapped eQTL gene
qtl_selection <- function(ch, start, end) {
    meta.ko <- qtl.KO %>%
        mutate(label = paste0(pheno, ": ", gene_name, ", ", gene_definition)) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "KO") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    meta.mag <- qtl.MAG %>%
        mutate(label = paste0(pheno, ": ", genus)) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "MAG") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    meta.mag.akk <- qtl.MAG.akk %>%
        mutate(label = paste0(pheno, ": A.muciniphila")) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "MAG.Akk") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    meta.taxa <- qtl.taxa %>%
        mutate(label = pheno) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "taxa") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    lipid <- qtl.lipid %>%
        mutate(label = identifier) %>%
        mutate(pheno = identifier, chr = qtl.chr, peak_mbp = qtl.pos, lod = qtl.lod) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "lipid") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    lipid.ol <- qtl.lipid.ol %>%
        mutate(label = pheno) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "lipid") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    lipid <- rbind(lipid, lipid.ol)
    
    rna <- qtl.rna %>%
        mutate(label = paste0("eQTL: ", gene_name)) %>%
        mutate(chr = qtl_chr) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "eQTL") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    return(rbind(meta.ko, meta.mag, meta.mag.akk, meta.taxa, lipid, rna))
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
                  annotation_colors = list(Group = c(lipid = "#63C29C", 
                                                     eQTL= "#FFCB04", 
                                                     KO = "grey40", 
                                                     taxa = "grey70", 
                                                     MAG = "grey90",
                                                     MAG.Akk = "#3d84b8"))[1],
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
                          annotation_colors = list(Group = c(lipid = "#63C29C", 
                                                             eQTL = "#FFCB04", 
                                                             KO = "grey40", 
                                                             taxa = "grey70", 
                                                             MAG ="grey90",
                                                             MAG.Akk = "#3d84b8"))[1],
                          cellwidth = 20,
                          cellheight = 5,
                          fontsize_row = 5,
                          cutree_rows = cuttree, 
                          main = paste0("Allele Effects: Chr", ch, ": ", start, "-", end, "Mbp"))
    
    
    return(p_curtree)
}

### chr1
chr <- "1"
start <- 90
end <- 95
qtl_plot <- qtl_selection(chr,start,end)
qtl_plot <- qtl_plot[-c(42:44),]
qtl_plot$cor <- NA
for (i in 1:nrow(qtl_plot)) {
    this_cor <- cor.test(as.numeric(qtl_plot[i,5:12]), as.numeric(qtl_plot[28,5:12]), method = "pearson")
    qtl_plot$cor[i] <- this_cor$estimate
}
plot_out_chr1 <- qtl_plot %>%
    filter(cor > 0.7 | cor < -0.7) %>%
    qtl_heatmap(chr,start,end)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Akk_comapping_chr1.eps", height = 5, width = 7)
plot_out_chr1
dev.off()

### chr2
chr <- "2"
start <- 77
end <- 81
qtl_plot <- qtl_selection(chr,start,end)
qtl_plot$cor <- NA
for (i in 1:nrow(qtl_plot)) {
    this_cor <- cor.test(as.numeric(qtl_plot[i,5:12]), as.numeric(qtl_plot[66,5:12]), method = "pearson")
    qtl_plot$cor[i] <- this_cor$estimate
}
plot_out_chr2 <- qtl_plot %>%
    filter(group != "MAG.Akk" | group == "MAG.Akk" & lod > 7) %>%
    filter(cor > 0.6 | cor < -0.6) %>%
    qtl_heatmap(chr,start,end)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Akk_comapping_chr2.eps", height = 5, width = 12)
plot_out_chr2
dev.off()

### chr7
chr <- "7"
start <- 126
end <- 131
qtl_plot <- qtl_selection(chr,start,end)
qtl_plot$cor <- NA
for (i in 1:nrow(qtl_plot)) {
    this_cor <- cor.test(as.numeric(qtl_plot[i,5:12]), as.numeric(qtl_plot[33,5:12]), method = "pearson")
    qtl_plot$cor[i] <- this_cor$estimate
}
plot_out_chr7 <- qtl_plot %>%
    filter(group != "MAG.Akk" | group == "MAG.Akk" & lod > 5.5) %>%
    filter(cor > 0.6 | cor < -0.6) %>%
    qtl_heatmap(chr,start,end)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Akk_comapping_chr7.eps", height = 3, width = 8)
plot_out_chr7
dev.off()

### chr12
chr <- "12"
start <- 55
end <- 63
qtl_plot <- qtl_selection(chr,start,end)
qtl_plot <- qtl_plot[-c(44:45),]
qtl_plot$cor <- NA
for (i in 1:nrow(qtl_plot)) {
    this_cor <- cor.test(as.numeric(qtl_plot[i,5:12]), as.numeric(qtl_plot[45,5:12]), method = "pearson")
    qtl_plot$cor[i] <- this_cor$estimate
}
plot_out_chr12 <- qtl_plot %>%
    filter(group != "MAG.Akk" | group == "MAG.Akk" & lod > 5.5) %>%
    filter(cor > 0.6 | cor < -0.6) %>%
    qtl_heatmap(chr,start,end)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Akk_comapping_chr12.eps", height = 3, width = 8)
plot_out_chr12
dev.off()

### chr15
chr <- "15"
start <- 75
end <- 79
qtl_plot <- qtl_selection(chr,start,end)
qtl_plot$cor <- NA
for (i in 1:nrow(qtl_plot)) {
    this_cor <- cor.test(as.numeric(qtl_plot[i,5:12]), as.numeric(qtl_plot[48,5:12]), method = "pearson")
    qtl_plot$cor[i] <- this_cor$estimate
}
plot_out_chr15 <- qtl_plot %>%
    filter(group != "MAG.Akk" | group == "MAG.Akk" & lod > 7) %>%
    filter(cor > 0.6 | cor < -0.6) %>%
    qtl_heatmap(chr,start,end)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Akk_comapping_chr15.eps", height = 6, width = 10)
plot_out_chr15
dev.off()

### 2021-08-18, co-mapping pglyrp1 gene eQTL
chr <- "8"
start <- 10.5
end <- 14.5
qtl_plot <- qtl_selection(chr,start,end)
qtl_plot$cor <- NA
for (i in 1:nrow(qtl_plot)) {
    this_cor <- cor.test(as.numeric(qtl_plot[i,5:12]), as.numeric(qtl_plot[41,5:12]), method = "pearson")
    qtl_plot$cor[i] <- this_cor$estimate
}
plot_out_pglyrp1 <- qtl_plot %>%
    filter(group == "KO" | pheno == "ENSMUSG00000030413") %>%
    qtl_heatmap(chr,start,end)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/pglyrp1_comap.eps", height = 6, width = 10)
plot_out_pglyrp1
dev.off()
