### Purpose: generate figures for DO paper
###          1. mass spectrometry for OL detection
###          2. 
###          3. 

# load required libraries
library(dplyr)
library(tibble)
library(pheatmap)
library(ggplot2)
library(ggsci)
library(viridis)
library(gridExtra)
library(RColorBrewer)

options(stringsAsFactors = F)
source("/Users/rootqz/R/QZ_functions.R")

### load input data
data.dir <- "~/Desktop/DO/experiment/LC-MS_lipidomics_2020-09-25/analysis/"

master <- read.table(paste0(data.dir, "OLs_detection_2020-09-16.csv"), sep = ",", header = T)

# ms data for phosphate growth condition
ms.po4 <- read.table(paste0(data.dir, "Final_Results_phosphate_study.csv"), sep = ",", header = T, check.names = F)
# rename column name
names(ms.po4)[8:24] <- c(master$Sample_ID[1:15], "blank1", "blank2")
names(ms.po4)

# ms data for AmEVs 
ms.ev <- read.table(paste0(data.dir, "Final_Results_AmEV.csv"), sep = ",", header = T, check.names = F)
# rename column name
names(ms.ev)[8:26] <- master$Sample_ID[16:34]
names(ms.ev)

ms.ol <- read.table(paste0(data.dir, "Final_Results_OL_standard.csv"), sep = ",", header = T, check.names = F)
names(ms.ol)

### heatmap plotting
# plot heatmap of cells in PO4 conditions
po4.plot <- ms.po4 %>%
    filter(`Lipid Class` == "OL") %>%
    select(c(master$Sample_ID[1:15], "blank1", "blank2")) %>%
    log10()

# rename with lipid features
tmp <- ms.po4 %>%
    filter(`Lipid Class` == "OL") %>%
    select(c(1,2,3,5))

names(tmp) <- c("a", "b", "c", "d")
rownames(po4.plot) <- paste0(tmp$d, "_", tmp$a, "_", tmp$b, "_", tmp$c)

po4.plot.AM <- po4.plot[,c(1:3,6:8,11:13)]
names(po4.plot.AM) <- c("Amucinphia#1_20uM", "Amucinphia#2_20uM","Amucinphia#3_20uM",
                        "Amucinphia#1_200uM", "Amucinphia#2_200uM","Amucinphia#3_200uM",
                        "Amucinphia#1_2000uM", "Amucinphia#2_2000uM","Amucinphia#3_2000uM")

my_col <- data.frame(phosphate = rep(c("20uM", "200uM", "2000uM"), each=3))
rownames(my_col) <- names(po4.plot.AM)

my_color <- list(Lipid.class = c("OL" = pal_jco()(5)[1],
                                 "PE" = pal_jco()(5)[2],
                                 "PG" = pal_jco()(5)[4],
                                 "CL" = pal_jco()(5)[5],
                                 "UNK" = "grey75"),
                 phosphate = c("20uM" = viridis(10)[10],
                               "200uM" = viridis(10)[7],
                               "2000uM" = viridis(10)[4])
)

# heatmap for OL, Amucin in different phosphate condition
po4.out <- pheatmap(po4.plot.AM,
                    cluster_cols = T,
                    annotation_col = my_col,
                    annotation_colors = my_color,
                    cutree_cols = 3,
                    cellwidth = 25,
                    cellheight = 10,
                    color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(400),
                    border_color = NA)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/heatmap_AmPO4_OL.eps", height = 9, width = 8)
po4.out
dev.off()

# heatmap for all lipids, Amucin in different phosphate condition
po4.all <- ms.po4 %>%
    select(c(1:6, master$Sample_ID[c(1:3,6:8,11:13)]))
po4.all$mean <- rowMeans(subset(po4.all, select = c(7:15)))
po4.all.order <- po4.all[order(po4.all$mean, decreasing = T),]
po4.all.order.top <- po4.all.order[1:100,]

po4.all.order.top$`Lipid Class`[which(po4.all.order.top$`Lipid Class` == "")] <- "UNK"
po4.all.order.top$`Lipid Class`[which(po4.all.order.top$`Lipid Class` == "PE-NMe")] <- "PE"
po4.all.order.top$`Lipid Class`[which(po4.all.order.top$`Lipid Class` == "PE-NMe2")] <- "PE"
po4.all.order.top$`Lipid Class`[which(po4.all.order.top$`Lipid Class` == "MLCL")] <- "CL"


po4.plot.all <- log10(po4.all.order.top[,c(master$Sample_ID[c(1:3,6:8,11:13)])])
rownames(po4.plot.all) <- paste0(po4.all.order.top$`Lipid Class`, "_", po4.all.order.top$`Retention Time (min)`, "_", po4.all.order.top$`Quant Ion`)
names(po4.plot.all) <- c("Amucinphia#1_20uM", "Amucinphia#2_20uM","Amucinphia#3_20uM",
                        "Amucinphia#1_200uM", "Amucinphia#2_200uM","Amucinphia#3_200uM",
                        "Amucinphia#1_2000uM", "Amucinphia#2_2000uM","Amucinphia#3_2000uM")


my_row <- data.frame(`Lipid class` = po4.all.order.top$`Lipid Class`)
rownames(my_row) <- rownames(po4.plot.all)

my_col <- data.frame(phosphate = rep(c("20uM", "200uM", "2000uM"), each=3))
rownames(my_col) <- names(po4.plot.all)

my_color <- list(Lipid.class = c("OL" = pal_jco()(5)[1],
                                   "PE" = pal_jco()(5)[2],
                                   "PG" = pal_jco()(5)[4],
                                   "CL" = pal_jco()(5)[5],
                                   "UNK" = "grey75"),
                 phosphate = c("20uM" = viridis(10)[10],
                               "200uM" = viridis(10)[7],
                               "2000uM" = viridis(10)[4])
                 )

po4.out.all <- pheatmap(po4.plot.all,
                        cluster_cols = T,
                        annotation_col = my_col,
                        annotation_row = my_row,
                        annotation_colors = my_color,
                        cutree_cols = 3,
                        cellwidth = 25,
                        cellheight = 4,
                        color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(400),
                        border_color = NA,
                        show_rownames = F)
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/heatmap_AmPO4_top100_lipids.eps", height = 9, width = 7)
po4.out.all
dev.off()

### heatmap for AmEVs
ev.all <- ms.ev %>%
    select(c(1:6, master$Sample_ID[c(16:18,22:24,28:30)]))
ev.all$mean <- rowMeans(subset(ev.all, select = c(7:15)))
ev.all.order <- ev.all[order(ev.all$mean, decreasing = T),]
ev.all.order.top <- ev.all.order[1:500,]

ev.all.order.top$`Lipid Class`[which(ev.all.order.top$`Lipid Class` == "")] <- "UNK"
ev.all.order.top$`Lipid Class`[which(ev.all.order.top$`Lipid Class` == "PE-NMe")] <- "PE"
ev.all.order.top$`Lipid Class`[which(ev.all.order.top$`Lipid Class` == "PE-NMe2")] <- "PE"
ev.all.order.top$`Lipid Class`[which(ev.all.order.top$`Lipid Class` == "MLCL")] <- "CL"


ev.plot.all <- log10(ev.all.order.top[,c(master$Sample_ID[c(16:18,22:24,28:30)])])
rownames(ev.plot.all) <- paste0(ev.all.order.top$`Lipid Class`, "_", ev.all.order.top$`Retention Time (min)`, "_", ev.all.order.top$`Quant Ion`)
names(ev.plot.all) <- rep(c("AmEV#1", "AmEV#2","AmEV#3"), each=3)

my_row <- data.frame(`Lipid class` = ev.all.order.top$`Lipid Class`)
rownames(my_row) <- rownames(ev.plot.all)

my_color <- list(Lipid.class = c("OL" = pal_jco()(5)[1],
                                 "PE" = pal_jco()(5)[2],
                                 "PG" = pal_jco()(5)[4],
                                 "CL" = pal_jco()(5)[5],
                                 "PC" = pal_jco()(10)[6],
                                 "TG" = pal_jco()(10)[7],
                                 "SP" = pal_jco()(10)[9],
                                 "UNK" = "grey75"))

ev.out.all <- pheatmap(ev.plot.all,
                        cluster_cols = F,
                        annotation_row = my_row,
                        annotation_colors = my_color,
                        cellwidth = 25,
                        cellheight = 1,
                        color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(400),
                        border_color = NA,
                        show_rownames = F)
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/heatmap_AmEV_top500_lipids.eps", height = 9, width = 7)
ev.out.all
dev.off()

# synthetic OL standard

### colonized mouse cecum content
ms.mouse <- read.table("~/Desktop/DO/experiment/LC-MS_lipidomics/Lipidomics_AllSignals.csv", sep = ",", header = T)

