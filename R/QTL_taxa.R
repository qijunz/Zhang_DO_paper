### Purpose: QTL mapping for metagenomic taxa
#          1. 5 taxa level, including phylum, order, class, family, genus

# load required libraries
library(dplyr)

options(stringsAsFactors = FALSE)
source("/Users/rootqz/R/QZ_functions.R")

# set path
setwd("/Users/rootqz/Desktop/DO/")

mb.dir <- "/Users/rootqz/Desktop/DO/metagenomics/2020/"

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

taxa_abun <- readRDS(paste0(mb.dir, "data/DO.summed.TPM.taxa.rds"))
taxa_abun_rankz <- readRDS(paste0(mb.dir, "data/DO.summed.TPM.taxa.rankZ.rds"))
taxa_info <- readRDS(paste0(mb.dir, "data/DO_1.9M_NRGeneSet_taxaInfo.rds"))

### QTL scan
phe.df <- taxa_abun_rankz

# generate list of gene names
phenotypes <- colnames(phe.df)

# set covar
covar <- cbind(pheno_clin[,"mouse",drop=F],
               sex = (pheno_clin$sex == "M")*1,
               dietday = pheno_clin$diet_days,
               DOwave2 = (pheno_clin$DOwave == 2)*1,
               DOwave4 = (pheno_clin$DOwave == 4)*1)

covar$mouse <- NULL

# create output table:
qtl.scan.df <- as.data.frame(matrix(ncol = 7, nrow = 1)) # contains scan 1 output
colnames(qtl.scan.df)<- c("lodindex", "lodcolumn", "chr", "pos", "lod", "ci_lo", "ci_hi")

# qtl scan1 output and other list
qtl.scan.list <- list()

##### Loop for QTL scan
# track progress
tally <- 0
total <- length(phenotypes)

for (phename in phenotypes) {
    tally = sum(tally, 1)
    cat ("On phenotype ", tally, " of", total, "\n")
    
    phe <- phe.df[,phename, drop=FALSE]
    phe <- na.omit(phe)
    
    # QTL scan1:
    qtl.scan1 <- scan1(probs, phe, K, covar, cores = 0)
    
    wgs.peaks <- find_peaks(qtl.scan1, pmap, threshold = 5, prob=0.95)
    qtl.scan.df <<- rbind(qtl.scan.df, wgs.peaks)
    
    ##### Saving list ######
    qtl.scan.list[[phename]] <- qtl.scan1
}

# combine output into a single scan1 object with cbind()
qtl.scan <- do.call("cbind", qtl.scan.list)

qtl.scan.df<- na.omit(qtl.scan.df)
qtl.scan.df$lodindex <- NULL                      
colnames(qtl.scan.df)<- c("pheno", "chr", "peak_mbp", "lod", "ci_lo", "ci_hi")

# add marker info to each QTL
qtl.peak <- qtl.scan.df
for (i in 1:nrow(qtl.peak)) {
    qtl.peak$marker[i] <- rownames(markers)[which(markers$chr == qtl.peak$chr[i] & markers$pos == qtl.peak$peak_mbp[i])]
}

qtl.peak <- merge(qtl.peak, taxa_info, by.x = "pheno", by.y = "new", all.x = T)

saveRDS(qtl.scan, file = paste0(mb.dir, "result/qtl_scan1_out_metagenomic_taxa.rds"))
saveRDS(qtl.peak, file = paste0(mb.dir, "result/qtl_scan1_peaks_minLOD5_metagenomic_taxa.rds"))

### BLUPs QTL genetic effect
# a new table adding blup to each QTL
qtl.blup <- data.frame(matrix(nrow = nrow(qtl.peak), ncol = 13))
names(qtl.blup) <- c("A", "B", "C", "D", "E", "F", "G", "H", "intercept", "sex", "dietday", "DOwave2", "DOwave4")
qtl.blup <- cbind(qtl.peak, qtl.blup)

# loop for BLUPs estimation
total <- nrow(qtl.peak)
for (i in 1:nrow(qtl.peak)) {
    # track progress
    cat ("On peak ", i, " of", total, "\n")
    
    peak <- qtl.peak[i,,drop=F]
    phename <- peak$pheno
    chr = peak$chr
    
    phe <- taxa_abun_rankz[,phename, drop=FALSE]
    phe <- na.omit(phe)
    
    # QTL scan1 from list:
    this.scan <- qtl.scan[,phename,drop=F]
    
    ######### QTL effects blup scan  #############
    print(paste("effect scan:", phename, "chr", chr))
    this.blup <- scan1blup(probs[,chr], phe, K[[chr]], covar, cores = 0)
    
    # add blup to QTL table
    qtl.blup[i,11:23] <- this.blup[which(rownames(this.blup) == qtl.peak$marker[i]),,drop=F]
    
}

saveRDS(qtl.blup, file = paste0(mb.dir, "result/qtl_blup_minLOD5_taxa.rds"))
# write.table(qtl.blup, file = paste0(mb.dir, "result/qtl_blup_minLOD5_taxa.tsv"), sep = "\t", quote = F, row.names = F)

### p_Bacteroidetes/p_Firmicutes ratio QTL
phe <- as.data.frame(rankZ(taxa_abun$p_Bacteroidetes/taxa_abun$p_Firmicutes))
rownames(phe) <- rownames(taxa_abun)
names(phe) <- "B2F"

chr = "15"
qtl.scan1 <- scan1(probs, phe, K, covar, cores = 0)
peaks <- find_peaks(qtl.scan1, pmap, threshold = 5, prob=0.95)
this.blup <- scan1blup(probs[,chr], phe, K[[chr]], covar, cores = 0)
b2f_out <- this.blup[which(rownames(this.blup) == "15_64221007"),,drop=F] %>% allele_eff_plot()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/allele_effects_taxa_B2F.eps", height = 2.5, width = 6)
b2f_out
dev.off()

