### Purpose: generate figures for DO paper
###          1. candidate genes for metagenomic hotspots

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

# load ens gene
ens_gene <- readRDS("~/Desktop/DO/RNA-seq/analysis/eQTL/data/mouseGenome/Mus_musculus.GRCm38.84.gene.rds")

# KO QTL hotspots
qtl.KO.hs <- qtl.KO %>% filter(lod > 9)

### snps association
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

for (i in 1:nrow(qtl.KO.hs)) {
    KO <- qtl.KO.hs$pheno[i]
    chr <- qtl.KO.hs$chr[i]
    start <- qtl.KO.hs$ci_lo[i]-2
    end <- qtl.KO.hs$ci_hi[i]+2
    
    phe <- KO_abun_rankz[,KO,drop=F]
    
    # QTL scan1 from list:
    qtl_scan1 <- scan1(probs, phe, K, covar, cores = 0)
    
    # query snps
    query_variants <- create_variant_query_func(snp_db)
    # find SNPs, keeping ALL OF THEM:
    out_snps <- scan1snps(probs[,chr], pmap, phe, 
                          K[chr], addcovar=covar,
                          query_func=query_variants, chr=chr, start=start, end=end,
                          keep_all_snps=TRUE)
    
    # top SNPs
    topSNPs <- top_snps(out_snps$lod, out_snps$snpinfo, drop = 1.5, show_all_snps = FALSE)
    topSNPs$phe <- KO
    
    if (i == 1) {
        topSNPs_out <- topSNPs
    } else {
        topSNPs_out <<- rbind(topSNPs_out, topSNPs)
    }

}

topSNPs_out_merge <- topSNPs_out %>%
    left_join(KEGG_gene, by = c("phe" = "KO")) %>%
    left_join(ens_gene[,c("gene_id", "gene_name", "gene_biotype")], by = c("ensembl_gene" = "gene_id"))

