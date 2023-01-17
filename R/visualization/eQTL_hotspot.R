### Purpose: generate figures for DO paper
###          1. eQTL hotspots

# load required libraries
source("/Users/rootqz/R/QZ_functions.R")
options(stringsAsFactors = FALSE)

### load eQTL
rna.dir <- "/Users/rootqz/Desktop/DO/RNA-seq/analysis/eQTL/"
qtl.rna <- readRDS(paste0(rna.dir, "result/qtl_blupsEff_intstine_gbrs_diploid_gene_withGeneAnno.rds"))

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
    rna <- qtl.rna %>%
        mutate(label = paste0("eQTL: ", gene_name)) %>%
        mutate(chr = qtl_chr) %>%
        select(pheno, chr, peak_mbp, lod, A, B, C, D, E, F, G, H, label) %>%
        mutate(group = "eQTL") %>%
        filter(chr == ch, peak_mbp > start, peak_mbp < end)
    
    return(rbind(rna))
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


chr18 <- qtl_selection("18",83,87)
chr18.fig <- qtl_heatmap(chr18,"18",83,87)

pdf(file = "~/Downloads/chr18.pdf", width = 7, height = 30)
chr18.fig
dev.off()

chr14 <- qtl_selection("14",105,109)
chr14.fig <- qtl_heatmap(chr14,"14",105,109)

pdf(file = "~/Downloads/chr14.pdf", width = 7, height = 30)
chr14.fig
dev.off()

chr10 <- qtl_selection("10",25,29)
chr10.fig <- qtl_heatmap(chr10,"10",25,29)

pdf(file = "~/Downloads/chr10.pdf", width = 7, height = 30)
chr10.fig
dev.off()

chr8 <- qtl_selection("8",17,23)
chr8.fig <- qtl_heatmap(chr8,"8",17,23)

pdf(file = "~/Downloads/chr8.pdf", width = 7, height = 30)
chr8.fig
dev.off()

chr4 <- qtl_selection("4",20,23)
chr4.fig <- qtl_heatmap(chr4,"4",20,23)

pdf(file = "~/Downloads/chr4.pdf", width = 7, height = 20)
chr4.fig
dev.off()

##################
load(file = "/Users/rootqz/Desktop/DO/data/DO_QTL_allInput_without_Phe_df_20181205.RData")
markers <- readRDS("/Users/rootqz/Desktop/DO/data/DO_69k_markers.rds")

# WHERE:
# pheno_clin        - clinical phenotypes (500 x 170)
# pheno_clin_dict   - dictionary of clinical phenotypes
# K                 - list of "loco" kinship matrices (500 x 500)
# pmap              - physical map for the 69k grid
# probs             - genotype probabilities
# markers           - data frame for 69k grid with marker, chr, pos, cM, bp

gene.abun <- readRDS(paste0(rna.dir, "data/gbrs.gene.tpm.filter.rds"))
gene.abun.trans <- readRDS(paste0(rna.dir, "data/gbrs.gene.tpm.filter.rankz.rds"))

meta.rna <- readRDS(paste0(rna.dir, "data/metadata_rna.rds"))
# load ensembl gene annotation, release 84
ens.gene <- readRDS(paste0(rna.dir, "data/mouseGenome/Mus_musculus.GRCm38.84.gene.rds"))

# set up covar
covar <- cbind(meta.rna[,"Sample",drop=F],
               sex = (meta.rna$sex == "M")*1,
               DOwave2 = (meta.rna$DOwave == 2)*1,
               DOwave4 = (meta.rna$DOwave == 4)*1,
               DOwave5 = (meta.rna$DOwave == 5)*1,
               DUI1 = (meta.rna$RNAPlate == 1)*1,
               DUI2 = (meta.rna$RNAPlate == 2)*1,
               DUI3 = (meta.rna$RNAPlate == 3)*1)
rownames(covar) <- covar$Sample
covar$Sample <- NULL

# chr18 hotspot genes
qtl.chr18 <- qtl.rna %>%
    filter(lod > 7.2) %>%
    filter(qtl_chr == "18", peak_mbp > 82, peak_mbp < 87)

pca.chr18 <- gene.abun.trans[,qtl.chr18$pheno] %>%
    as.matrix() %>%
    prcomp()

pca <- pca.chr18$x %>%
    as.data.frame()

qtl.scan <- scan1(genoprobs = probs, pheno = pca[,1,drop=F], kinship = K, addcovar = covar, cores = 0)
plot_scan1(qtl.scan, pmap)

gene.med.list <- qtl.chr18 %>%
    filter(cis)

covar.med <- covar %>%
    mutate(mouse = rownames(covar)) %>%
    left_join(gene.abun.trans[,gene.med.list$pheno[5],drop=F] %>% 
                  mutate(mouse = rownames(gene.abun.trans)), by = "mouse") %>%
    column_to_rownames("mouse")

qtl.scan.med <- scan1(genoprobs = probs, pheno = pca[,1,drop=F], kinship = K, addcovar = covar.med, cores = 0)
plot_scan1(qtl.scan.med, pmap)

# chr14 hotspot genes
qtl.chr14 <- qtl.rna %>%
    filter(lod > 7.2) %>%
    filter(qtl_chr == "14", peak_mbp > 105, peak_mbp < 109)

pca.chr18 <- gene.abun.trans[,qtl.chr18$pheno] %>%
    as.matrix() %>%
    prcomp()

pca <- pca.chr18$x %>%
    as.data.frame()

qtl.scan <- scan1(genoprobs = probs, pheno = pca[,1,drop=F], kinship = K, addcovar = covar, cores = 0)
plot_scan1(qtl.scan, pmap)

gene.med.list <- qtl.chr18 %>%
    filter(cis)

covar.med <- covar %>%
    mutate(mouse = rownames(covar)) %>%
    left_join(gene.abun.trans[,gene.med.list$pheno[5],drop=F] %>% 
                  mutate(mouse = rownames(gene.abun.trans)), by = "mouse") %>%
    column_to_rownames("mouse")

qtl.scan.med <- scan1(genoprobs = probs, pheno = pca[,1,drop=F], kinship = K, addcovar = covar.med, cores = 0)
plot_scan1(qtl.scan.med, pmap)

##################



