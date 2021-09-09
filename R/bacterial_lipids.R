### Purpose: correlation between KO/taxa vs. cecum lipids, co-mapping example

# load required libraries
library(ggplot2)
library(gridExtra)

source("/Users/rootqz/R/QZ_functions.R")
options(stringsAsFactors = FALSE)

# load meta table
meta_tbl <- read.table("~/Desktop/DO/data/metadata_500DO.tsv", sep = "\t", header = T)

# set path
mb.dir <- "/Users/rootqz/Desktop/DO/metagenomics/2020/"

### Load QTL table
qtl.KO <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_kegg_orthology.rds"))
qtl.taxa <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD5_taxa.rds"))
qtl.MAG <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_MAGs.rds"))
qtl.lipid <- readRDS("~/Desktop/DO/RNA-seq/analysis/eQTL/data/qtl.summary.cecum.vl.rds")

phe_KO <- readRDS("~/Desktop/DO/metagenomics/2020/data/DO.summed.TPM.KEGG.orthology.rds")
phe_taxa <- readRDS("~/Desktop/DO/metagenomics/2020/data/DO.summed.TPM.taxa.rds")
phe_MAG <- readRDS("~/Desktop/DO/metagenomics/2020/data/DO.MAGs.GenomeCov.libSizeNorm.rds")
phe_lipid_ori <- readRDS("~/Desktop/DO/metabolomics/data_raw/attie_cecum_lipids_normalized.rds")

mag.anno <- read.table(paste0(mb.dir, "data/MAG_bins_metaTable.tsv"), sep = "\t")

lipids <- names(phe_lipid_ori)[12:3382]
lipids <- lipids[which(!startsWith(lipids, "UNK"))]
phe_lipid <- phe_lipid_ori[which(rownames(phe_lipid_ori) %in% rownames(phe_MAG)),lipids]

# MAG taxa annotation
taxa_anno_all <- mag.anno[which(mag.anno$bin %in% names(phe_MAG)), c("bin", "phylum", "class", "order", "family", "genus", "species")]
taxa_anno <- taxa_anno_all[, c("bin", "phylum")]
rownames(taxa_anno) <- taxa_anno$bin
taxa_anno$bin <- NULL

# lipid annotation at category level
lipid_anno_all <- qtl.lipid[which(qtl.lipid$identifier %in% names(phe_lipid)),c("identifier", "class", "category")] %>% unique()
lipid_anno <- lipid_anno_all[,c("identifier", "category")] %>% unique()
rownames(lipid_anno) <- lipid_anno$identifier
lipid_anno$identifier <- NULL

# lipid annotation at class level
lipid_anno_class <- lipid_anno_all[,c("identifier", "class")] %>% unique()
rownames(lipid_anno_class) <- lipid_anno_class$identifier
lipid_anno_class$identifier <- NULL

# correlation between MAG vs. lipid
phe_MAG_cor <- phe_MAG[which(rownames(phe_MAG) %in% rownames(phe_lipid)),]
phe_lipid_cor <- phe_lipid[,rownames(lipid_anno)]
cor_lipid_MAG <- cor(phe_MAG_cor, phe_lipid_cor, method = "spearman", use = "complete.obs")

saveRDS(cor_lipid_MAG, file = "~/Desktop/DO/omics/2021/data/cor_MAG_cecumLipid.rds")

#####
# heatmap for all
pdf("~/Desktop/ReyLab/paper/DO_metagenomic/figure/cor_MAGvsCL_lipid.pdf", height = 8, width = 15)
pheatmap(cor_lipid_MAG,
         show_rownames = F,
         show_colnames = F,
         clustering_method = "ward.D2",
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(120), 
         breaks = seq(-0.6, 0.6, by = 0.01),
         border_color = NA,
         annotation_row = taxa_anno,
         annotation_col = lipid_anno,
         annotation_colors = list(category = c(Fatty.Acyl = pal_startrek()(5)[1], 
                                               Glycerolipid = pal_startrek()(5)[2], 
                                               Sphingolipid = pal_startrek()(5)[3], 
                                               Phospholipid = pal_startrek()(5)[4], 
                                               UNK = "grey70"),
                                  phylum = c(Actinobacteria = pal_rickandmorty()(5)[1], 
                                             Bacteroidetes = pal_rickandmorty()(5)[5], 
                                             Firmicutes = pal_rickandmorty()(5)[3], 
                                             Tenericutes = pal_rickandmorty()(5)[4], 
                                             Verrucomicrobia = pal_rickandmorty()(5)[2],
                                             unclassified = "grey70")),
         cutree_rows = 5,
         cutree_cols = 6)
dev.off()

# Sphingolipid
lipids_Sphingolipid <- lipids[which(lipids %in% qtl.lipid$identifier[which(qtl.lipid$category == "Sphingolipid")])]
phe_lipid_Sphingolipid_cor <- phe_lipid[,lipids_Sphingolipid]

cor_lipid_MAG_Sphingolipid <- cor(phe_MAG_cor, phe_lipid_Sphingolipid_cor, method = "spearman", use = "complete.obs")

lipid_anno_class_Sphingolipid <- lipid_anno_class[which(rownames(lipid_anno_class) %in% lipids_Sphingolipid),,drop=F]

pheatmap(cor_lipid_MAG_Sphingolipid,
         show_rownames = F,
         show_colnames = F,
         clustering_method = "ward.D2",
         annotation_row = taxa_anno,
         annotation_col = lipid_anno_class_Sphingolipid,
         annotation_colors = list(phylum = c(Actinobacteria = pal_rickandmorty()(5)[1], 
                                             Bacteroidetes = pal_rickandmorty()(5)[5], 
                                             Firmicutes = pal_rickandmorty()(5)[3], 
                                             Tenericutes = pal_rickandmorty()(5)[4], 
                                             Verrucomicrobia = pal_rickandmorty()(5)[2],
                                             unclassified = "grey70")),
         cutree_rows = 5,
         cutree_cols = 5)

# Phospholipid
lipids_Phospholipid <- lipids[which(lipids %in% qtl.lipid$identifier[which(qtl.lipid$category == "Phospholipid")])]
phe_lipid_Phospholipid_cor <- phe_lipid[,lipids_Phospholipid]

cor_lipid_MAG_Phospholipid <- cor(phe_MAG_cor, phe_lipid_Phospholipid_cor, method = "spearman", use = "complete.obs")

lipid_anno_class_Phospholipid <- lipid_anno_class[which(rownames(lipid_anno_class) %in% lipids_Phospholipid),,drop=F]

pheatmap(cor_lipid_MAG_Phospholipid,
         show_rownames = F,
         show_colnames = F,
         clustering_method = "ward.D2",
         annotation_row = taxa_anno,
         annotation_col = lipid_anno_class_Phospholipid,
         annotation_colors = list(phylum = c(Actinobacteria = pal_rickandmorty()(5)[1], 
                                             Bacteroidetes = pal_rickandmorty()(5)[5], 
                                             Firmicutes = pal_rickandmorty()(5)[3], 
                                             Tenericutes = pal_rickandmorty()(5)[4], 
                                             Verrucomicrobia = pal_rickandmorty()(5)[2],
                                             unclassified = "grey70")),
         cutree_rows = 5,
         cutree_cols = 5)
#####

###
heatmap_all <- pheatmap(cor_lipid_MAG,
                        show_rownames = F,
                        show_colnames = F,
                        clustering_method = "ward.D2",
                        color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(120), 
                        breaks = seq(-0.6, 0.6, by = 0.01),
                        border_color = NA,
                        annotation_row = taxa_anno,
                        annotation_col = lipid_anno,
                        annotation_colors = list(category = c(Fatty.Acyl = pal_startrek()(5)[1], 
                                                              Glycerolipid = pal_startrek()(5)[2], 
                                                              Sphingolipid = pal_startrek()(5)[3], 
                                                              Phospholipid = pal_startrek()(5)[4], 
                                                              UNK = "grey70"),
                                                 phylum = c(Actinobacteria = pal_rickandmorty()(5)[1], 
                                                            Bacteroidetes = pal_rickandmorty()(5)[5], 
                                                            Firmicutes = pal_rickandmorty()(5)[3], 
                                                            Tenericutes = pal_rickandmorty()(5)[4], 
                                                            Verrucomicrobia = pal_rickandmorty()(5)[2],
                                                            unclassified = "grey70")),
                        cutree_rows = 5,
                        cutree_cols = 6)


tree_row <- cutree(heatmap_all$tree_row, k = 5)
tree_col <- cutree(heatmap_all$tree_col, k = 6)

pheatmap(cor_lipid_MAG[names(tree_row)[which(tree_row == 1)],names(tree_col)[which(tree_col == 6)]],
         show_rownames = F,
         show_colnames = F,
         clustering_method = "ward.D2",
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(120), 
         breaks = seq(-0.6, 0.6, by = 0.01),
         border_color = NA,
         annotation_row = taxa_anno,
         annotation_col = lipid_anno,
         annotation_colors = list(category = c(Fatty.Acyl = pal_startrek()(5)[1], 
                                               Glycerolipid = pal_startrek()(5)[2], 
                                               Sphingolipid = pal_startrek()(5)[3], 
                                               Phospholipid = pal_startrek()(5)[4], 
                                               UNK = "grey70"),
                                  phylum = c(Actinobacteria = pal_rickandmorty()(5)[1], 
                                             Bacteroidetes = pal_rickandmorty()(5)[5], 
                                             Firmicutes = pal_rickandmorty()(5)[3], 
                                             Tenericutes = pal_rickandmorty()(5)[4], 
                                             Verrucomicrobia = pal_rickandmorty()(5)[2],
                                             unclassified = "grey70"))
         )

# merge lipid annotation with cluster information
lipid_anno_cluster <- as.data.frame(tree_col)
names(lipid_anno_cluster) <- "cluster"
lipid_anno_cluster$cluster <- as.character(lipid_anno_cluster$cluster)
lipid_anno_cluster$identifier <- rownames(lipid_anno_cluster)
lipid_anno_cluster <- lipid_anno_all %>% 
    left_join(lipid_anno_cluster, by = c("identifier" = "identifier"))

# merge MAG annotation with cluster information
taxa_anno_cluster <- as.data.frame(tree_row)
names(taxa_anno_cluster) <- "cluster"
taxa_anno_cluster$cluster <- as.character(taxa_anno_cluster$cluster)
taxa_anno_cluster$bin <- rownames(taxa_anno_cluster)
taxa_anno_cluster <- taxa_anno_all %>% 
    left_join(taxa_anno_cluster, by = c("bin" = "bin"))

### Enrichment for lipid class
rownames(lipid_anno_all) <- lipid_anno_all$identifier
lipid_all_class <- unique(lipid_anno_cluster[,2:3])
fisher_exact_test <- function(allL, myL){
    
    allTbl <- lipid_anno_all[allL,]
    myTbl <- lipid_anno_all[myL,]
    
    allTbl_freq <- as.data.frame(table(lipid_anno_all$class))
    names(allTbl_freq) <- c("class", "allFreq")
    myTbl_freq <- as.data.frame(table(myTbl$class))
    names(myTbl_freq) <- c("class", "myFreq")
    
    tbl_freq <- merge(allTbl_freq, myTbl_freq, by = "class", all = T)
    tbl_freq[is.na(tbl_freq)] <- 0
    
    tbl_out <- cbind(tbl_freq, as.data.frame(matrix(ncol = 5, nrow = nrow(tbl_freq))))
    names(tbl_out)[4:8] <- c("pathInt", "pathNoInt", "backInt", "backNoInt", "pval")
    
    for (i in 1:nrow(tbl_out)){
        tbl_out$pathInt[i] <- tbl_out$myFreq[i]
        tbl_out$pathNoInt[i] <- tbl_out$allFreq[i] - tbl_out$myFreq[i]
        tbl_out$backInt[i] <- sum(tbl_out$myFreq) - tbl_out$myFreq[i]
        tbl_out$backNoInt[i] <- sum(tbl_out$allFreq) - sum(tbl_out$myFreq) - (tbl_out$allFreq[i] - tbl_out$myFreq[i])
        
        fisher_tbl <- matrix(as.integer(tbl_out[i,4:7]), byrow = TRUE, 2,2)
        fisher_result <- fisher.test(fisher_tbl, alternative = "greater", conf.int = TRUE, conf.level = 0.95)
        
        tbl_out$pval[i] <- fisher_result$p.value
    }
    
    tbl_out <- tbl_out[which(tbl_out$myFreq > 1),]
    tbl_out$bonferroni <- p.adjust(tbl_out$pval, method = "bonferroni")
    tbl_out$bh <- p.adjust(tbl_out$pval, method = "BH")
    tbl_out <- tbl_out[order(tbl_out$bh),]
    tbl_out <- tbl_out %>% left_join(lipid_all_class, by = c("class" = "class"))
    
    return(tbl_out)
}

enrichment_plot <- function(tbl_out) {
    tbl_out$label <- paste0(tbl_out$category, ": ", tbl_out$class)
    tbl_out <- tbl_out %>% filter(bh < 0.1)
    tbl_out$label <- factor(tbl_out$label, levels = rev(tbl_out$label))
    
    plot_out <- ggplot(tbl_out) +
        geom_bar(aes(y = -log10(bh), x = label), stat="identity", fill = "grey50") +
        theme_classic() +
        scale_y_continuous(expand = expansion(add=c(0,0))) +
        coord_flip() +
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size = 20)) +
        ylab("-log10(1/p value after BH correction)") +
        xlab("Cecum Lipid Class")
}

allL <- lipid_anno_cluster$identifier
fisher_lipid <- list()
for (i in 1:6) {
    myL <- lipid_anno_cluster$identifier[which(lipid_anno_cluster$cluster == as.character(i))]
    tbl_out <- fisher_exact_test(allL, myL)
    fisher_lipid[[i]] <- tbl_out
    
    # plot bar
    tbl_out$label <- paste0(tbl_out$category, ": ", tbl_out$class)
    tbl_out <- tbl_out %>% filter(bh < 0.1)
    tbl_out$label <- factor(tbl_out$label, levels = rev(tbl_out$label))
    
    plot_out <- ggplot(tbl_out) +
        geom_bar(aes(y = -log10(bh), x = label), stat="identity", fill = "grey50") +
        theme_classic() +
        scale_y_continuous(limits = c(0,10), expand = expansion(add=c(0,0))) +
        coord_flip() +
        theme(axis.text = element_text(size = 15),
              axis.title = element_blank(),
              axis.text.y = element_blank()) +
        ylab("") +
        xlab("")
    
    setEPS()
    postscript(paste0("~/Desktop/ReyLab/paper/DO_metagenomic/figure/cor_MAGvsCL_lipid_fisher_enrichment_", i,".eps"), height = 3, width = 10)
    plot_out
    dev.off()
}

# put all together
for (i in 1:6) {
    myL <- lipid_anno_cluster$identifier[which(lipid_anno_cluster$cluster == as.character(i))]
    tbl_out <- fisher_exact_test(allL, myL)
    fisher_lipid[[i]] <- tbl_out
    
    # plot bar
    tbl_out$label <- paste0(tbl_out$category, ": ", tbl_out$class)
    tbl_out <- tbl_out %>% filter(bh < 0.1)
    tbl_out$label <- factor(tbl_out$label, levels = rev(tbl_out$label))
    tbl_out$group <- paste0("cluster", i)
    
    if (i == 1) {
        tbl_out_all <- tbl_out
    } else {
        tbl_out_all <<- rbind(tbl_out_all, tbl_out)
        }
    
}

tbl_out_all <- tbl_out_all %>%
    mutate(cluster = case_when(
        group == "cluster1" ~ "cluster4",
        group == "cluster2" ~ "cluster1",
        group == "cluster3" ~ "cluster5",
        group == "cluster4" ~ "cluster2",
        group == "cluster5" ~ "cluster6",
        group == "cluster6" ~ "cluster3"
    ))

scaleFUN <- function(x) sprintf("%.0f", x)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/cor_MAGvsCL_lipid_fisher_enrichment_all.eps", height = 10, width = 12)
ggplot(tbl_out_all) +
    geom_bar(aes(y = -log10(bh), x = label, fill = category), stat="identity") +
    scale_fill_manual(values = c(pal_startrek()(5)[1:4],  "grey70")) +
    theme_classic() +
    scale_y_continuous(trans = "log10",
                       breaks = c(0,1,2,3,4,5,10,20,80),
                       expand = expansion(add=c(0,0))
                       ) +
    coord_flip() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_blank()) +
    facet_grid(cluster ~ ., scales = "free_y", space ="free_y") +
    ylab("") +
    xlab("")
dev.off()

### 2021-08-22, get the table of # of positive/negative correlations in each lipid cluster vs. MAG cluster

# MAG true cluster label
#
# cluster 1 -> cluster 4
# cluster 2 -> cluster 1
# cluster 3 -> cluster 2
# cluster 4 -> cluster 5
# cluster 5 -> cluster 3

# Lipid true cluster label
#
# cluster 1 -> cluster 4
# cluster 2 -> cluster 1
# cluster 3 -> cluster 5
# cluster 4 -> cluster 2
# cluster 5 -> cluster 6
# cluster 6 -> cluster 3

for (i in 1:5) {
    for (j in 1:6) {
        dim <- dim(cor_lipid_MAG[names(tree_row)[which(tree_row == i)],names(tree_col)[which(tree_col == j)]])
        n_all <- dim[1] * dim[2]
        n_posi <- length(which(cor_lipid_MAG[names(tree_row)[which(tree_row == i)],names(tree_col)[which(tree_col == j)]] > 0))
        cat(paste0("taxa: ", i, "; lipid: ", j, " - ", n_posi/n_all, "\n"))
    }
}




