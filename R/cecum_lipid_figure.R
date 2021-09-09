# Purpose: modify Vanessa's figures for paper

# load required library
library("ggplot2")
library("dplyr")
library("plyr")
library("pheatmap")
library("RColorBrewer")
library("grid")
library("gridExtra")

options(stringsAsFactors = F)

# define colors
coon_blue <- "#2CA7DF"
coon_grey <- "#566977"
coon_purp <- "#955CA5"
coon_red <- "#EC6B63"
coon_turq <- "#63C29C"
coon_yel <- "#FFCB04"

# load QTL peaks
meta.qtl.dir <- "/Users/rootqz/Desktop/DO/metagenomics/2020/result/"

qtl.KO <- readRDS(paste0(meta.qtl.dir, "qtl_blup_minLOD6_kegg_orthology.rds"))
qtl.genus <- readRDS(paste0(meta.qtl.dir, "qtl_blup_minLOD6_ncbi_genus.rds"))
qtl.MAG <- readRDS(paste0(meta.qtl.dir, "qtl_blup_minLOD6_MAGs.rds"))
qtl.lipid <- read.table("/Users/rootqz/Desktop/DO/RNA-seq/analysis/eQTL/data/qtl.summary.cecum.vl.csv", 
                        sep = ",", header = T, row.names = 1)
qtl.rna <- readRDS("/Users/rootqz/Desktop/DO/RNA-seq/analysis/eQTL/result/qtl_blupsEff_intstine_gbrs_diploid_gene_withGeneAnno.rds")


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
    
    meta.mag <- qtl.MAG %>%
        mutate(label = paste0(pheno, ": ", genus)) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "MAG") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    meta.genus <- qtl.genus %>%
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
    
    rna <- qtl.rna %>%
        mutate(label = paste0("eQTL: ", gene_name)) %>%
        mutate(chr = qtl_chr) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "eQTL") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    return(rbind(meta.ko, meta.mag, meta.genus, lipid, rna))
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
                  annotation_colors = list(Group = c(lipid = coon_turq, 
                                                     eQTL=coon_yel, 
                                                     KO = "grey40", 
                                                     taxa = "grey70", 
                                                     MAG="grey90"))[1],
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
                          annotation_colors = list(Group = c(lipid = coon_turq, 
                                                             eQTL=coon_yel, 
                                                             KO = "grey40", 
                                                             taxa = "grey70", 
                                                             MAG="grey90"))[1],
                          cellwidth = 20,
                          cellheight = 5,
                          fontsize_row = 5,
                          cutree_rows = cuttree, 
                          main = paste0("Allele Effects: Chr", ch, ": ", start, "-", end, "Mbp"))
    
    
    return(p_curtree)
}

# test function
chr <- "17"
start <- 30
end <- 34

test <- qtl_selection(chr,start,end)
# test.plot <- test[c(1,6,22),]
# rownames(test.plot) <- NULL
out.fig <- qtl_heatmap(test,chr,start,end)
out.fig

pdf(paste0("~/Desktop/ReyLab/paper/DO_metagenomic/figures/old_2020-08/fig/chr", chr, "-", start, "-", end, ".pdf"), 
    width = 10, height = 8)
out.fig
dev.off()


