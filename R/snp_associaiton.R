### Purpose: SNP association for metagenomic QTLs

# load required libraries
library(qtl2)
library(qtl2convert)
library(qtl2ggplot)
library(qtl2pattern)
library(dplyr)
library(tidyverse)

source("/Users/rootqz/R/QZ_functions.R")

options(stringsAsFactors = FALSE)

# setup working dir
setwd("/Users/rootqz/Desktop/DO/metagenomics/2020/")

# load dataset containing all files:
load(file = "/Users/rootqz/Desktop/DO/data/DO_QTL_allInput_without_Phe_df_20181205.RData")
markers <- readRDS("/Users/rootqz/Desktop/DO/data/DO_69k_markers.rds")

# WHERE:
# pheno_clin        - clinical phenotypes (500 x 170)
# pheno_clin_dict   - dictionary of clinical phenotypes
# K                 - list of "loco" kinship matrices (500 x 500)
# pmap              - physical map for the 69k grid
# probs             - genotype probabilities
# markers           - data frame for 69k grid with marker, chr, pos, cM, bp

# load snp and gene 
# SNP database
snp_db <- "/Users/rootqz/Desktop/DO/data/cc_variants.sqlite"
# mouse genes database
gene_db <- "/Users/rootqz/Desktop/DO/data/mouse_genes_mgi.sqlite"

KO_abun <- readRDS("data/DO.summed.TPM.KEGG.orthology.rds")
KO_abun_rankz <- readRDS("data/DO.summed.TPM.KEGG.orthology.rankZ.rds")

genus_abun <- readRDS("data/DO.summed.TPM.ncbi.genus.100genes.rds")

# load KO QTL peak
qtl.KO.peak <- readRDS("result/qtl_blup_minLOD6_kegg_orthology.rds")
qtl.KO.scan <- readRDS("result/qtl_scan1_out_metagenomic_kegg_orthology.rds")

# set covariance
covar <- cbind(pheno_clin[,"mouse",drop=F],
               sex = (pheno_clin$sex == "M")*1,
               dietday = pheno_clin$diet_days,
               DOwave2 = (pheno_clin$DOwave == 2)*1,
               DOwave4 = (pheno_clin$DOwave == 4)*1)

covar$mouse <- NULL

# filter QTL with LOD > 9
qtl.KO.peak.filter <- qtl.KO.peak %>% filter(lod > 9)

# generate out table for top SNPs and patterns
tbl.out.snp <- data.frame(matrix(nrow = 0, ncol = 21))
names(tbl.out.snp) <- c("snp_id", "chr", "pos", "alleles", "sdp", "ensembl_gene", "consequence",
                     "A_J", "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ",
                     "type", "index", "interval", "on_map", "pheno", "lod")

tbl.out.pattern <- data.frame(matrix(nrow = 0, ncol = 8))
names(tbl.out.pattern) <- c("pattern", "max_lod", "max_pos", "pheno", "pct", "min_lod", "sdp", "max_snp")

for (i in 1:nrow(qtl.KO.peak.filter)) {
    peak <- qtl.KO.peak.filter[i,,drop=F]
    phename <- peak$pheno
    chr <- peak$chr
    
    phe <- KO_abun_rankz[,phename,drop=F]
    phe <- na.omit(phe)
    
    # QTL scan1 from list:
    qtl_scan1 <- qtl.KO.scan[,phename,drop=F]
    
    # query snps
    query_variants <- create_variant_query_func(snp_db)
    # find SNPs, keeping ALL OF THEM:
    out_snps <- scan1snps(probs[,chr], pmap, phe, 
                          K[chr], addcovar=covar,
                          query_func=query_variants, chr=chr, start=peak$ci_lo, end=peak$ci_hi,
                          keep_all_snps=TRUE)
    
    # query gene
    query_genes <- create_gene_query_func(gene_db)
    genes <- query_genes(chr = chr, peak$ci_lo, peak$ci_hi)
    
    # top SNPs
    topSNPs <- top_snps(out_snps$lod, out_snps$snpinfo, drop = 1.5, show_all_snps = FALSE)
    
    # snp association
    snpinfo <- index_snps(pmap, out_snps$snpinfo)                # create index
    snpprobs <- genoprob_to_snpprob(probs, snpinfo)              # convert to SNP probs
    scan_snppr <- scan1(snpprobs, phe, kinship = K[[chr]], addcovar = covar)
    
    # get top SNPs:
    top_snps_tbl <- top_snps_all(scan_snppr, snpinfo)
    
    # concat into output table
    tbl.out.snp <<- rbind(tbl.out.snp, top_snps_tbl)
    
    # return patterns, snpinfo, snpprobs
    patterns <- arrange(summary(top_snps_tbl), desc(max_lod))
    
    # concat into output table
    tbl.out.pattern <<- rbind(tbl.out.pattern, patterns)
}


