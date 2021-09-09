### Purpose: generate figures for DO paper
###          1. Akk 16s in founders
###          2. 
###          3. 


# load required libraries
library(dplyr)
library(tibble)
library(qtl2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

source("/Users/rootqz/R/QZ_functions.R")

taxa <- read.table("~/Desktop/DO/omics/integration/16S/founder/taxa_genus.txt", sep = "\t", header = T, check.names = F)
taxa_akk <- taxa %>% 
    filter(`OTU ID` == "k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Verrucomicrobiaceae;g__Akkermansia") %>%
    t() %>%
    as.data.frame()
names(taxa_akk) <- "Akkermansia_16s"
taxa_akk <- taxa_akk[-1,,drop=F]
taxa_akk$mouse <- rownames(taxa_akk)
taxa_akk$strain <- NA

taxa_akk$Akkermansia_16s <- log10(as.numeric(taxa_akk$Akkermansia_16s))

# rename
taxa_akk$strain[which(startsWith(taxa_akk$mouse, "A"))] <- "AJ"
taxa_akk$strain[which(startsWith(taxa_akk$mouse, "B6"))] <- "B6"
taxa_akk$strain[which(startsWith(taxa_akk$mouse, "129"))] <- "129"
taxa_akk$strain[which(startsWith(taxa_akk$mouse, "NOD"))] <- "NOD"
taxa_akk$strain[which(startsWith(taxa_akk$mouse, "NZO"))] <- "NZO"
taxa_akk$strain[which(startsWith(taxa_akk$mouse, "Cast"))] <- "CAST"
taxa_akk$strain[which(startsWith(taxa_akk$mouse, "Pwk"))] <- "PWK"
taxa_akk$strain[which(startsWith(taxa_akk$mouse, "WSB"))] <- "WSB"

taxa_akk$strain <- factor(taxa_akk$strain, levels = names(CCcolors))

out <- ggplot(taxa_akk) +
    geom_boxplot(aes(x = strain, y = Akkermansia_16s, fill = strain, color = strain)) +
    theme_classic() +
    scale_fill_manual(values = CCcolors) +
    scale_color_manual(values = CCcolors) +
    xlab("DO founder strains") +
    ylab("Relative abundance of A.muciniphila (log10)") +
    theme(axis.title = element_text(size = 12))

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/16S_akk.eps", height = 4, width = 6)
out
dev.off()

