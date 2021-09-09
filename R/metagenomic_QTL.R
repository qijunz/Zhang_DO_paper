### Purpose: generate figures for DO paper
###          1. KO QTL Manhattan plot
###          2. taxa (genus level) QTL Manhattan plot
###          3. MAGs QTL Manhattan plot


# load required libraries
library("qtl2")
library("ggplot2")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("grid")
library("gridExtra")

mb.dir <- "/Users/rootqz/Desktop/DO/metagenomics/2020/"

KEGG_gene <- read.table(paste0(mb.dir, "data/KEGG_KO2gene2EC.tsv"), sep = "\t", header = T, quote = "")
KEGG_path <- read.table(paste0(mb.dir, "data/KEGG_KO2path.tsv"), sep = "\t", header = T, quote = "")

# trim KEGG_path
KEGG_path$path_definition <- gsub("\\s*\\[[^\\]+\\]","",KEGG_path$path_definition) %>%
    trimws()
KEGG_path <- unique(KEGG_path[,c("KO", "path_definition")]) %>%
    na.omit()

KO_abun <- readRDS(paste0(mb.dir, "data/DO.summed.TPM.KEGG.orthology.rds"))
KO_abun_rankz <- readRDS(paste0(mb.dir, "data/DO.summed.TPM.KEGG.orthology.rankZ.rds"))

# load QTL results
qtl.scan.KO <- readRDS(paste0(mb.dir, "result/qtl_scan1_out_metagenomic_kegg_orthology.rds"))
qtl.scan.taxa <- readRDS(paste0(mb.dir, "result/qtl_scan1_out_metagenomic_taxa.rds"))
qtl.scan.MAG <- readRDS(paste0(mb.dir, "result/qtl_scan1_out_MAGs.rds"))

qtl.KO <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_kegg_orthology.rds"))
qtl.taxa <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD5_taxa.rds"))
qtl.MAG <- readRDS(paste0(mb.dir, "result/qtl_blup_minLOD6_MAGs.rds"))

markers <- readRDS("/Users/rootqz/Desktop/DO/data/DO_69k_markers.rds")

pseudo.peak <- markers %>%
    group_by(chr) %>%
    mutate(max = max(bp)/1000000, min = min(bp)/1000000) %>%
    select(chr, max, min) %>%
    unique() %>%
    reshape(idvar = "chr", v.names = "peak_mbp", 
            direction = "long", varying = c("min", "max"), new.row.names = 1:40) %>%
    select(-time) %>%
    mutate(lod = 0) %>%
    as.data.frame()

color1 <- c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19")
color2 <- c("2", "4", "6", "8", "10", "12", "14", "16", "18", "X")

qtl.peak.KO.plot <- qtl.KO[,c("chr", "peak_mbp", "lod")]
qtl.peak.KO.plot <- rbind(qtl.peak.KO.plot, pseudo.peak)
qtl.peak.KO.plot$color <- NA

for (i in 1:nrow(qtl.KO)) {
    if (qtl.peak.KO.plot$chr[i] %in%  color1){
        qtl.peak.KO.plot$color[i] <- "color1"
    } else {
        qtl.peak.KO.plot$color[i] <- "color2"
    }
}

qtl.peak.KO.plot$chr <- factor(qtl.peak.KO.plot$chr, levels = c("1","2","3","4","5","6","7","8","9","10",
                                                                "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))

# calculate cumulative marker position
mh.plot <- qtl.peak.KO.plot %>%
    # compute chr size
    group_by(chr) %>%
    summarise(chr_len = max(peak_mbp)) %>%
    
    # calculate cumulative position of each chr
    mutate(total = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # add this infor to extra column
    left_join(qtl.peak.KO.plot, ., by=c("chr" = "chr")) %>%
    
    # add a cumulative position of each marker
    arrange(chr, peak_mbp) %>%
    mutate(pos_cum = peak_mbp+total)

# make x axis with chr display
axis.df <- mh.plot %>%
    group_by(chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

out.plot.KO <- ggplot(mh.plot, aes(x=pos_cum, y=lod)) +
    geom_point(aes(color=color), alpha=1, size=0.7) +
    
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(6, 12), expand = expansion(add=c(0,0))) +
    
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "white"),
          legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    scale_color_manual(values = c("color1" = "#13547a", "color2" = "#80d0c7")) +
    ylab("LOD") + xlab("Mouse Chromosome")


setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/Manhattan_plot_KO_QTL.eps", height = 4, width = 16)
out.plot.KO
dev.off()

### taxa QTL
qtl.peak.taxa.plot <- qtl.taxa[,c("chr", "peak_mbp", "lod")]
qtl.peak.taxa.plot <- rbind(qtl.peak.taxa.plot, pseudo.peak)
qtl.peak.taxa.plot$color <- NA

for (i in 1:nrow(qtl.taxa)) {
    if (qtl.peak.taxa.plot$chr[i] %in%  color1){
        qtl.peak.taxa.plot$color[i] <- "color1"
    } else {
        qtl.peak.taxa.plot$color[i] <- "color2"
    }
}

qtl.peak.taxa.plot$chr <- factor(qtl.peak.taxa.plot$chr, levels = c("1","2","3","4","5","6","7","8","9","10",
                                                                      "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))

# calculate cumulative marker position
mh.plot <- qtl.peak.taxa.plot %>%
    # compute chr size
    group_by(chr) %>%
    summarise(chr_len = max(peak_mbp)) %>%
    
    # calculate cumulative position of each chr
    mutate(total = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # add this infor to extra column
    left_join(qtl.peak.taxa.plot, ., by=c("chr" = "chr")) %>%
    
    # add a cumulative position of each marker
    arrange(chr, peak_mbp) %>%
    mutate(pos_cum = peak_mbp+total)

# make x axis with chr display
axis.df <- mh.plot %>%
    group_by(chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

out.plot.taxa <- ggplot(mh.plot, aes(x=pos_cum, y=lod)) +
    geom_point(aes(color=color), alpha=1, size=0.7) +
    
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(6, 12), expand = expansion(add=c(0,0))) +
    
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "white"),
          legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    scale_color_manual(values = c("color1" = "#FF7B89", "color2" = "#BA5082")) +
    ylab("LOD") + xlab("Mouse Chromosome")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Manhattan_plot_taxa_QTL.eps", height = 4, width = 16)
out.plot.taxa
dev.off()


# fisher's exact test
### fisher's exact test
fisher_exact_test <- function(allKO, myKO){
    
    allTbl <- KEGG_path[which(KEGG_path$KO %in% allKO),]
    myTbl <- KEGG_path[which(KEGG_path$KO %in% myKO),]
    
    allTbl_freq <- as.data.frame(table(allTbl$path_definition))
    names(allTbl_freq) <- c("path", "allFreq")
    myTbl_freq <- as.data.frame(table(myTbl$path_definition))
    names(myTbl_freq) <- c("path", "myFreq")
    
    tbl_freq <- merge(allTbl_freq, myTbl_freq, by = "path", all = T)
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
    
    return(tbl_out)
}

##
all.ko <- names(KO_abun)
my.ko <- qtl.blup.KO %>%
    filter(chr == 15, peak_mbp > 61, peak_mbp < 65, C < 0) %>%
    select(pheno)

tmp <- fisher_exact_test(all.ko, my.ko$pheno) %>%
    mutate(logPval = -log10(bh)) %>%
    filter(bh < 0.2)

fish_plot <- ggplot(tmp, aes(reorder(path, logPval), logPval)) +
    geom_bar(stat = "identity", fill= "#F2BC94", width = 0.9) +
    coord_flip() + theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.text.y = element_text(size = 12, colour = c("white")),
          strip.background = element_rect(fill = "grey70")) +
    scale_y_continuous(expand = expansion(add=c(0,0))) +
    labs(x = "", y= "-log10(1/BH-Pvalue)")


setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/KO_QTL_chr15_hs_fisher.eps", height = 5, width = 10)
fish_plot
dev.off()

# chr15 hoptspot zoom in
qtl.peak.chr15 <- qtl.blup.KO %>%
    filter(chr == 15, peak_mbp > 61, peak_mbp < 65)
qtl.peak.chr15$group <- NA
qtl.peak.chr15$group[which(grepl("sporulation", qtl.peak.chr15$gene_definition))] <- "spore"

ko.motility <- as.character(KEGG_path$KO[which(KEGG_path$path_definition == "Bacterial motility proteins")])
ko.cellGrowth <- as.character(KEGG_path$KO[which(KEGG_path$path_definition == "Cell growth")])

qtl.peak.chr15$group[which(qtl.peak.chr15$pheno %in% ko.motility)] <- "motility"
# qtl.peak.chr15$group[which(qtl.peak.chr15$pheno %in% ko.cellGrowth)] <- "cell_growth"

qtl.peak.chr15$group[which(is.na(qtl.peak.chr15$group))] <- "others"

qtl.peak.chr15$group <- factor(qtl.peak.chr15$group, levels = c("others", "spore", "motility"))
qtl.peak.chr15 <- qtl.peak.chr15 %>%
    arrange(group)

chr15_zoomin <- ggplot(qtl.peak.chr15) +
    theme_bw() +
    geom_point(aes(x=peak_mbp, y=lod, color = group)) +
    scale_y_continuous(limits = c(6, 11.3), expand = c(0,0)) +
    scale_x_continuous(limits = c(61, 65), expand = c(0,0)) +
    scale_color_manual(values = c("spore" = "#F8766D", 
                                  "motility" = "#00A9FF",
                                  "others" = "grey80"), name = "KEGG Orthology") +
    theme(legend.position="none",
          axis.text = element_text(size = 10)) +
    xlab("") +
    ylab("")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/KO_QTL_chr15_hs_zoomin.eps", height = 7, width = 7)
chr15_zoomin
dev.off()


# Chr15 hotspots ~65Mbp

spor.qtl <- qtl.blup.KO %>%
    filter(chr == "15", peak_mbp > 61, peak_mbp < 65) %>%
    filter(grepl("sporulation", gene_definition))

motility.qtl <- qtl.blup.KO[which(qtl.blup.KO$pheno %in% ko.motility),] %>%
    filter(chr == "15", peak_mbp > 61, peak_mbp < 65) 

combine.qtl <- rbind(spor.qtl, motility.qtl)

spor.ko <- spor.qtl$pheno

spor.gene <- gene.anno %>%
    filter(KO %in% spor.ko)

spor.heatmap <- spor.qtl %>%
    select("A", "B", "C", "D", "E", "F", "G", "H") 
rownames(spor.heatmap) <- paste0(spor.qtl$pheno, ": ", spor.qtl$gene_definition)
colnames(spor.heatmap) <- names(CCcolors)

spor.plot <- pheatmap(spor.heatmap,
                      cluster_cols = F,
                      color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(400), 
                      breaks = seq(-2, 2, by = 0.01),
                      border_color = NA,
                      scale = "row")

mot.heatmap <- motility.qtl %>%
    select("A", "B", "C", "D", "E", "F", "G", "H") 
rownames(mot.heatmap) <- paste0(motility.qtl$pheno, ": ", motility.qtl$gene_definition)
colnames(mot.heatmap) <- names(CCcolors)

mot.plot <- pheatmap(mot.heatmap,
                      cluster_cols = F,
                      color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(400), 
                      breaks = seq(-2, 2, by = 0.01),
                      border_color = NA,
                      scale = "row")


setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/KO_QTL_chr15_hs_sporeKOs_heatmap.eps", height = 3, width = 10)
print(spor.plot)
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/KO_QTL_chr15_hs_motilityKOs_heatmap.eps", height = 3, width = 11)
print(mot.plot)
dev.off()



######## gene
# mouse genes database
gene_db <- "~/Desktop/ReyLab_Project/Data/CCdb/mouse_genes_mgi.sqlite"
### query gene
query_genes <- create_gene_query_func(gene_db)
genes <- query_genes(chr = 15, 61, 65)
genes <- genes[which(!startsWith(genes$mgiName, "predicted")),]

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/KO_QTL_chr15_hs_genes.eps", height = 3, width = 9)
plot_genes(genes = genes, xlab="", xlim=c(61,65))
dev.off()

