### Purpose: generate figures for DO paper
###          1. metagenomic stats

# load required libraries
options(stringsAsFactors = FALSE)
source("/Users/rootqz/R/QZ_functions.R")

assembly_stats <- read.table("~/Desktop/DO/metagenomics/2019-old/results/stats/assembly_stats.tsv", sep = "\t", header = T)

### plot assembly ratio vs. #PE reads
f1 <- ggplot(assembly_stats) +
    geom_point(aes(x = PE_reads/1000000, y=assembled_ratio)) +
    theme_bw() +
    xlab("Number of PE reads (Million per DO mice)") +
    ylab("de novo assembly ratio (%)") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15))

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/metagenomic_assebmly_ratio.eps", height = 5, width = 6.5)
f1
dev.off()

### plot allignment rate comparison
# load the csv file for 2.6M GeneSet
bowtie2_2.6M_stats <- read.csv("~/Desktop/DO/metagenomics/2019-old/scripts/scripts-bowtie2-2.6M_geneset_Nat/bowtie2_2.6M_genes_stats.csv", stringsAsFactors = F, check.names = F)
# rename mouse
bowtie2_2.6M_stats[,1] <- gsub("bowtie2_2.6M_genes_", "", bowtie2_2.6M_stats[,1])
rownames(bowtie2_2.6M_stats) <- bowtie2_2.6M_stats[,1]

# load the csv file for 1.9M GeneSet
bowtie2_1.9M_stats <- read.csv("~/Desktop/DO/metagenomics/2019-old/scripts/scripts-bowtie2-1.9M_geneset/bowtie2_1.9M_genes_stats.csv", stringsAsFactors = F, check.names = F)
# rename mouse
bowtie2_1.9M_stats[,1] <- gsub("bowtie2_1.9M_genes_", "", bowtie2_1.9M_stats[,1])
rownames(bowtie2_1.9M_stats) <- bowtie2_1.9M_stats[,1]

overall_1.9 <- bowtie2_1.9M_stats[,c("mouse", "overall_rate")]
names(overall_1.9)[2] <- c("overall_rate_1.9")
overall_2.6 <- bowtie2_2.6M_stats[,c("mouse", "overall_rate")]
names(overall_2.6)[2] <- c("overall_rate_2.6")

overall_plot_df <- merge(overall_1.9, overall_2.6, by = "mouse")
overall_plot_df <- overall_plot_df[which(overall_plot_df$overall_rate_1.9 != 0),]
overall_plot_df <- overall_plot_df[order(overall_plot_df$overall_rate_1.9),]
overall_plot_df$mouse <- factor(overall_plot_df$mouse, level = overall_plot_df$mouse)

plot_df <- melt(overall_plot_df, id = "mouse")

f2 <- ggplot(data = plot_df) + 
    theme_classic() + 
    geom_boxplot(aes(x= variable, y = value)) + 
    ylab("Overall alignment rate by Bowtie2 (%)") +
    theme(axis.title.x = element_blank(),
          axis.title = element_text(size=15),
          axis.text = element_text(size=12)) +
    scale_x_discrete(labels = c("Assembled gene\ncatalog", "Reference gene\ncatalog"))

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/metagenomic_align_comparison_boxplot.eps", height = 6, width = 4)
f2
dev.off()

