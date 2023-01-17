### Purpose: generate figures for DO paper
###          1. MAG reconstruction

# load required libraries
options(stringsAsFactors = FALSE)
source("/Users/rootqz/R/QZ_functions.R")

mb.dir <- "/Users/rootqz/Desktop/DO/metagenomics/2020/"

# load metagenome-assembled genome (MAG) annotation
mag.anno <- read.table(paste0(mb.dir, "data/MAG_bins_metaTable.tsv"), sep = "\t")

# plot MAG completeness vs. contamination

mag.plot <- mag.anno %>%
    filter(number_contigs < 1000, Contamination < 50) %>%
    ggplot() +
    geom_point(aes(x = Completeness, y = Contamination, color = number_contigs)) +
    theme_bw() +
    scale_color_gradient(limits = c(0, 1000), name = "Number of Contigs") +
    scale_x_continuous(limits = c(0,100)) +
    scale_y_continuous(limits = c(0,50)) +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12))

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/MAG_checkm_all.eps", height = 5, width = 7)
mag.plot
dev.off()

mag.hq.plot <- mag.anno %>%
    filter(Contamination < 5, Completeness > 50) %>%
    ggplot() +
    geom_point(aes(x = Completeness, y = Contamination, color = number_contigs)) +
    theme_bw() +
    scale_color_gradient(limits = c(0, 1000), name = "Number of Contigs") +
    scale_x_continuous(limits = c(50,100)) +
    scale_y_continuous(limits = c(0,5)) +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 12))

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/MAG_checkm_MQ.eps", height = 5, width = 7)
mag.hq.plot
dev.off()

### mash for Akkermansia MAGs pair-wise ANI
akk_mash <- read.table("~/Desktop/DO/omics/integration/binning/grouping_bins/mash-OSX64-v2.2/mash_out_Akk_HQ_bins.tsv", sep = "\t", header = T)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/MAG_Akk_ANI.eps", height = 9, width = 7)
pheatmap(akk_mash/1000,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         show_rownames = FALSE)
dev.off()


### Akkermansia abundance, metagenomic vs. 16s
taxa <- readRDS("~/Desktop/DO/data/JHK/Rey_DO_Data/MicrobiotaComposition/pheno_16s_taxon.RDS")
taxa <- as.data.frame(taxa)

taxa_16s <- taxa[which(rownames(taxa) %in% rownames(taxa_abun)),"g_Akkermansia",drop=F]
names(taxa_16s) <- "g_Akkermansia_16s"
taxa_meta <- taxa_abun[,"g_Akkermansia",drop=F]
names(taxa_meta) <- "g_Akkermansia_metagenomic"
taxa_meta$g_Akkermansia_metagenomic <- log10(taxa_meta$g_Akkermansia_metagenomic)

taxa_combined <- cbind(taxa_16s, taxa_meta)

hist_plot <- ggplot(taxa_combined) + 
    geom_point(aes(x = g_Akkermansia_metagenomic, y = g_Akkermansia_16s)) + 
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    xlab("metagenomic abundance of g_Akkermansia (log10)") +
    ylab("16S rRNA abundance of g_Akkermansia")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Akk_metagenomic_vs_16s.eps", height = 6, width = 6)
hist_plot
dev.off()

