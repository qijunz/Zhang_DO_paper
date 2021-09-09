### Purpose: script for supplementary tables/notes
### Created: 2021-08-18

library(qtl2)
library(pheatmap)

options(stringsAsFactors = FALSE)
source("/Users/rootqz/R/QZ_functions.R")

wgs_stats <- read.table("~/Desktop/DO/metagenomics/2020/data/WGS_assembly_stats.tsv", sep = "\t", header = T)

### Sup Tab for QTL
qtl_KO <- readRDS("~/Desktop/DO/metagenomics/2020/result/qtl_blup_minLOD6_kegg_orthology.rds")
qtl_taxa <- readRDS("~/Desktop/DO/metagenomics/2020/result/qtl_blup_minLOD5_taxa.rds")
qtl_lipid <- readRDS("~/Desktop/DO/RNA-seq/analysis/eQTL/data/qtl.summary.cecum.vl.rds")
qtl_rna <- readRDS("~/Desktop/DO/RNA-seq/analysis/eQTL/result/qtl_blupsEff_intstine_gbrs_diploid_gene_withGeneAnno.rds")

# write out csv table
write.table(qtl_KO, "~/Desktop/ReyLab/paper/DO_metagenomic/table/Supplementary_Table1_DO_GutMicrobiome_Function_QTL.tsv", 
            sep = "\t", quote = F, row.names = F)
write.table(qtl_taxa, "~/Desktop/ReyLab/paper/DO_metagenomic/table/Supplementary_Table2_DO_GutMicrobiome_Taxon_QTL.tsv", 
            sep = "\t", quote = F, row.names = F)
write.table(qtl_lipid, "~/Desktop/ReyLab/paper/DO_metagenomic/table/Supplementary_Table3_DO_CecumLipid_QTL.tsv", 
            sep = "\t", quote = F, row.names = F)
write.table(qtl_rna, "~/Desktop/ReyLab/paper/DO_metagenomic/table/Supplementary_Table4_DO_SmallIntestine_eQTL.tsv", 
            sep = "\t", quote = F, row.names = F)

### Sup Tabl for metagenomic annotation
# 1.9 million metagenes annotation
anno_mb_gene <- readRDS("~/Desktop/DO/metagenomics/2020/data/DO_1.9M_NRGeneSet_KEGGanno_NCBIanno.rds")
# filter out low abundance
phe_gene <- readRDS("~/Desktop/DO/omics/integration/data/DO_1.9M_metagene_with_GhostKOALA_anno_prokaryotes_only_SCnorm_min10reads_0.1sample_rankZ_20190405.rds")
anno_mb_gene_filter <- anno_mb_gene[which(anno_mb_gene$gene %in% names(phe_gene)),]

# metagenome-assembled genome (MAG) annotation
anno_mb_mag <- read.table("~/Desktop/DO/metagenomics/2020/data/MAG_bins_metaTable.tsv", sep = "\t")

write.table(anno_mb_gene_filter, "~/Desktop/ReyLab/paper/DO_metagenomic/table/Supplementary_Table5_DO_GutMicrobiome_Metagene_Annotation.tsv", 
            sep = "\t", quote = F, row.names = F)
write.table(anno_mb_mag, "~/Desktop/ReyLab/paper/DO_metagenomic/table/Supplementary_Table6_DO_GutMicrobiome_MAG_Annotation.tsv", 
            sep = "\t", quote = F, row.names = F)



