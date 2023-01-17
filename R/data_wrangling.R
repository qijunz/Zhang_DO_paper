# metagenomic data wrangling

library(dplyr)
library(tibble)
library(data.table)

options(stringsAsFactors = FALSE)

data.dir <- "/Users/rootqz/Desktop/DO/metagenomics/2020/data/"

# WGS reads metadata
meta.tbl <- read.table("~/Desktop/DO/metagenomics/2019/results/stats/assembly_stats.tsv", sep = "\t")

# add sample mix-up correction
# mixup DO samples
# label   should be	
# 53	  54	  simple mixup
# 54	  53	  simple mixup
# 83	  88	  simple mixup
# 85	  83	  simple mixup
# 88	  85	  simple mixup
# 360	  370	  simple mixup
# 370	  360	  simple mixup

meta.tbl$realLabel <- meta.tbl$mouse
meta.tbl$realLabel[which(meta.tbl$mouse == "DO053")] <- "DO054"
meta.tbl$realLabel[which(meta.tbl$mouse == "DO054")] <- "DO053"
meta.tbl$realLabel[which(meta.tbl$mouse == "DO083")] <- "DO088"
meta.tbl$realLabel[which(meta.tbl$mouse == "DO085")] <- "DO083"
meta.tbl$realLabel[which(meta.tbl$mouse == "DO088")] <- "DO085"
meta.tbl$realLabel[which(meta.tbl$mouse == "DO360")] <- "DO370"
meta.tbl$realLabel[which(meta.tbl$mouse == "DO370")] <- "DO360"

# write.table(meta.tbl, file = paste0(data.dir, "DO_WGS_reads_meta.tsv"), sep = "\t", quote = F, row.names = F)
meta.tbl <- read.table(paste0(data.dir, "DO_WGS_reads_meta.tsv"), sep = "\t", header = T)

# load RSEM TPM table, gene normalization using gene length (not effective length)
tpm <- fread(paste0(data.dir, "DO.1.9M.rsem.bowtie2.tpm.tsv")) %>% as.data.frame()
names(tpm)[1] <- "gene"

# load 1.9M gene annotation
gene.anno <- readRDS(paste0(data.dir, "DO_1.9M_NRGeneSet_KEGGanno_NCBIanno.rds"))

# only keep Bacteria gene
gene.anno <- gene.anno %>%
    filter(ncbi_superkingdom == "Bacteria")

# sample mix-up correction
if (all.equal(names(tpm)[-1], meta.tbl$mouse)) {
    names(tpm)[-1] <- meta.tbl$realLabel
}

# low quanlity sample removal
omic.meta <- read.table(paste0(data.dir, "DO_500mice_omics_info.tsv"), sep = "\t")
sample.list <- omic.meta$mouse[which(omic.meta$metagenomic)]
tpm <- tpm[,c("gene", sample.list),drop=F]

# rankZ function
rankZ = function(x) {
    x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
    return(qnorm(x))
}

### sum gene to KO level
tpm.ko <- tpm %>%
    left_join(gene.anno[,c("gene", "KO"),drop=F], by=c("gene" = "gene")) %>%
    filter(KO != "") %>%
    select(-gene) %>%
    group_by(KO) %>%
    filter(n() >= 10) %>%
    summarise_all(funs(sum)) %>%
    column_to_rownames(var="KO") %>%
    t() %>%
    as.data.frame()

# write.table(tpm.ko, file = paste0(data.dir, "DO.summed.TPM.KEGG.orthology.tsv"), sep = "\t", quote = F)
saveRDS(tpm.ko, file = paste0(data.dir, "DO.summed.TPM.KEGG.orthology.rds"))

# rankZ transformation
tpm.ko.rankz <- tpm.ko
for (i in 1:ncol(tpm.ko)){
    tpm.ko.rankz[,i] <- rankZ(tpm.ko[,i])
}
saveRDS(tpm.ko.rankz, file = paste0(data.dir, "DO.summed.TPM.KEGG.orthology.rankZ.rds"))


### sum gene to genus 
tpm.genus <- tpm %>%
    left_join(gene.anno[,c("gene", "ncbi_genus"),drop=F], by=c("gene" = "gene")) %>%
    filter(ncbi_genus != "") %>%
    select(-gene) %>%
    group_by(ncbi_genus) %>%
    filter(n() >= 100) %>%
    summarise_all(funs(sum)) %>%
    column_to_rownames(var="ncbi_genus") %>%
    t() %>%
    as.data.frame()

# write.table(tpm.genus, file = paste0(data.dir, "DO.summed.TPM.ncbi.genus.tsv"), sep = "\t", quote = F)
saveRDS(tpm.genus, file = paste0(data.dir, "DO.summed.TPM.ncbi.genus.rds"))
saveRDS(tpm.genus, file = paste0(data.dir, "DO.summed.TPM.ncbi.genus.100genes.rds"))

# rankZ transformation
tpm.genus.rankz <- tpm.genus
for (i in 1:ncol(tpm.genus)){
    tpm.genus.rankz[,i] <- rankZ(tpm.genus[,i])
}
saveRDS(tpm.genus.rankz, file = paste0(data.dir, "DO.summed.TPM.ncbi.genus.rankZ.rds"))


### MAGs matrix
mag.matrix <- readRDS("/Users/rootqz/Desktop/DO/integration/binning/result/ReconstGenomes_mashRep_CovQuanti_libSizeNorm_keepAllDO_mixupCorrected_noTrans.rds")

mag.matrix <- mag.matrix %>%
    t() %>%
    as.data.frame() %>%
    select(sample.list) %>%
    t() %>%
    as.data.frame()

saveRDS(mag.matrix, file = paste0(data.dir, "DO.MAGs.GenomeCov.libSizeNorm.rds"))

# rankZ transformation
mag.matrix.rankz <- mag.matrix
for (i in 1:ncol(mag.matrix)){
    mag.matrix.rankz[,i] <- rankZ(mag.matrix[,i])
}
saveRDS(mag.matrix.rankz, file = paste0(data.dir, "DO.MAGs.GenomeCov.libSizeNorm.rankZ.rds"))


### Akkermansia gene
akk.matrix <- tpm %>%
    left_join(gene.anno[,c("gene", "ncbi_genus"),drop=F], by=c("gene" = "gene")) %>%
    filter(ncbi_genus == "Akkermansia") %>%
    select(-ncbi_genus) %>%
    unique() %>%
    remove_rownames() %>%        # why need this??
    column_to_rownames(var="gene") %>%
    t() %>%
    as.data.frame()

# write.table(tpm.genus, file = paste0(data.dir, "DO.summed.TPM.ncbi.genus.tsv"), sep = "\t", quote = F)
saveRDS(akk.matrix, file = paste0(data.dir, "DO.Akkermansia.gene.TPM.rds"))

# rankZ transformation
akk.matrix.rankz <- akk.matrix
for (i in 1:ncol(akk.matrix)){
    akk.matrix.rankz[,i] <- rankZ(akk.matrix[,i])
}
saveRDS(akk.matrix.rankz, file = paste0(data.dir, "DO.Akkermansia.gene.TPM.rankZ.rds"))


############ 2021-05-13, taxa for other levels ############
# sum gene to genus level
tpm.g <- tpm %>%
    left_join(gene.anno[,c("gene", "ncbi_genus"),drop=F], by=c("gene" = "gene")) %>%
    filter(ncbi_genus != "") %>%
    select(-gene) %>%
    group_by(ncbi_genus) %>%
    filter(n() >= 100) %>%
    summarise_all(funs(sum)) %>%
    column_to_rownames(var="ncbi_genus") %>%
    t() %>%
    as.data.frame()

name_g <- as.data.frame(matrix(ncol = 2, nrow = ncol(tpm.g)))
names(name_g) <- c("old", "new")
name_g$old <- names(tpm.g)
name_g$new <- paste0("g_", name_g$old)
names(tpm.g) <- name_g$new
freq.g <- table(gene.anno$ncbi_genus) %>% as.data.frame()
name_g <- merge(name_g, freq.g, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_g$level <- "genus"

# sum gene to family level
tpm.f <- tpm %>%
    left_join(gene.anno[,c("gene", "ncbi_family"),drop=F], by=c("gene" = "gene")) %>%
    filter(ncbi_family != "") %>%
    select(-gene) %>%
    group_by(ncbi_family) %>%
    filter(n() >= 1000) %>%
    summarise_all(funs(sum)) %>%
    column_to_rownames(var="ncbi_family") %>%
    t() %>%
    as.data.frame()

name_f <- as.data.frame(matrix(ncol = 2, nrow = ncol(tpm.f)))
names(name_f) <- c("old", "new")
name_f$old <- names(tpm.f)
name_f$new <- paste0("f_", name_f$old)
names(tpm.f) <- name_f$new
freq.f <- table(gene.anno$ncbi_family) %>% as.data.frame()
name_f <- merge(name_f, freq.f, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_f$level <- "family"

# sum gene to order level
tpm.o <- tpm %>%
    left_join(gene.anno[,c("gene", "ncbi_order"),drop=F], by=c("gene" = "gene")) %>%
    filter(ncbi_order != "") %>%
    select(-gene) %>%
    group_by(ncbi_order) %>%
    filter(n() >= 1000) %>%
    summarise_all(funs(sum)) %>%
    column_to_rownames(var="ncbi_order") %>%
    t() %>%
    as.data.frame()

name_o <- as.data.frame(matrix(ncol = 2, nrow = ncol(tpm.o)))
names(name_o) <- c("old", "new")
name_o$old <- names(tpm.o)
name_o$new <- paste0("o_", name_o$old)
names(tpm.o) <- name_o$new
freq.o <- table(gene.anno$ncbi_order) %>% as.data.frame()
name_o <- merge(name_o, freq.o, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_o$level <- "order"

# sum gene to class level
tpm.c <- tpm %>%
    left_join(gene.anno[,c("gene", "ncbi_class"),drop=F], by=c("gene" = "gene")) %>%
    filter(ncbi_class != "") %>%
    select(-gene) %>%
    group_by(ncbi_class) %>%
    filter(n() >= 1000) %>%
    summarise_all(funs(sum)) %>%
    column_to_rownames(var="ncbi_class") %>%
    t() %>%
    as.data.frame()

name_c <- as.data.frame(matrix(ncol = 2, nrow = ncol(tpm.c)))
names(name_c) <- c("old", "new")
name_c$old <- names(tpm.c)
name_c$new <- paste0("c_", name_c$old)
names(tpm.c) <- name_c$new
freq.c <- table(gene.anno$ncbi_class) %>% as.data.frame()
name_c <- merge(name_c, freq.c, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_c$level <- "class"

# sum gene to phylum level
tpm.p <- tpm %>%
    left_join(gene.anno[,c("gene", "ncbi_phylum"),drop=F], by=c("gene" = "gene")) %>%
    filter(ncbi_phylum != "") %>%
    select(-gene) %>%
    group_by(ncbi_phylum) %>%
    filter(n() >= 1000) %>%
    summarise_all(funs(sum)) %>%
    column_to_rownames(var="ncbi_phylum") %>%
    t() %>%
    as.data.frame()

name_p <- as.data.frame(matrix(ncol = 2, nrow = ncol(tpm.p)))
names(name_p) <- c("old", "new")
name_p$old <- names(tpm.p)
name_p$new <- paste0("p_", name_p$old)
names(tpm.p) <- name_p$new
freq.p <- table(gene.anno$ncbi_phylum) %>% as.data.frame()
name_p <- merge(name_p, freq.p, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_p$level <- "phylum"

name_merge <- name_g %>%
    rbind(name_f) %>%
    rbind(name_o) %>%
    rbind(name_c) %>%
    rbind(name_p)

tpm_taxa <- tpm.g %>%
    cbind(tpm.f) %>%
    cbind(tpm.o) %>%
    cbind(tpm.c) %>%
    cbind(tpm.p)

tpm_taxa_rankZ <- tpm_taxa %>%
    mutate_all(rankZ)

saveRDS(name_merge, file = paste0(data.dir, "DO_1.9M_NRGeneSet_taxaInfo.rds"))

saveRDS(tpm_taxa, file = paste0(data.dir, "DO.summed.TPM.taxa.rds"))
saveRDS(tpm_taxa_rankZ, file = paste0(data.dir, "DO.summed.TPM.taxa.rankZ.rds"))

