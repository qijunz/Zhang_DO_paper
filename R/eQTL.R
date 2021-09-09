### Purpose: generate figures for DO paper
###          1. eQTL hotspots
###          2. 
###          3. 


# load required libraries
library(dplyr)
library(tibble)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)

source("/Users/rootqz/R/QZ_functions.R")

options(stringsAsFactors = FALSE)

### load eQTL
rna.dir <- "/Users/rootqz/Desktop/DO/RNA-seq/analysis/eQTL/"
qtl.rna <- readRDS(paste0(rna.dir, "result/qtl_blupsEff_intstine_gbrs_diploid_gene_withGeneAnno.rds"))

qtl.plot <- qtl.rna %>%
    # remove gene on chrY
    filter(gene_chr != "Y") %>%
    mutate(gene_chr = as.character(gene_chr)) %>%
    mutate(qtl_chr = factor(qtl_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X"))) %>%
    mutate(gene_chr = factor(gene_chr, levels = (c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")))) %>%
    na.omit()

# plot eQTL, QTL marker as x-axis and gene position as y-axis
# calculate cumulative marker position
qtl.plot.x <- qtl.plot %>%
    # compute chr size
    group_by(qtl_chr) %>%
    summarise(chr_len_qtl = max(peak_mbp)) %>%

    # calculate cumulative position of each chr
    mutate(total_qtl = cumsum(chr_len_qtl) - chr_len_qtl) %>%
    select(-chr_len_qtl) %>%
    
    # add this infor to extra column
    left_join(qtl.plot, ., by=c("qtl_chr" = "qtl_chr")) %>%
    
    # add a cumulative position of each marker
    arrange(qtl_chr, peak_mbp) %>%
    mutate(pos_cum_qtl = peak_mbp+total_qtl)

qtl.plot.x.y <- qtl.plot.x %>%
    # for gene pos
    group_by(gene_chr) %>%
    summarise(chr_len_gene = max(gene_pos)) %>%
    mutate(total_gene = cumsum(chr_len_gene) - chr_len_gene) %>%
    select(-chr_len_gene) %>%
    left_join(qtl.plot.x, ., by=c("gene_chr" = "gene_chr")) %>%
    arrange(gene_chr, gene_pos) %>%
    mutate(pos_cum_gene = gene_pos+total_gene)

# make x axis with chr display
axis.x <- qtl.plot.x.y %>%
    group_by(qtl_chr) %>%
    summarise(center=(max(pos_cum_qtl) + min(pos_cum_qtl))/2)

axis.y <- qtl.plot.x.y %>%
    group_by(gene_chr) %>%
    summarise(center=(max(pos_cum_gene) + min(pos_cum_gene))/2)

index_cis <- which(qtl.plot.x.y$cis)
index_trans <- which(!qtl.plot.x.y$cis)
qtl.plot.x.y$cis[index_cis] <- "cis"
qtl.plot.x.y$cis[index_trans] <- "trans"

eqtl.map <- ggplot(qtl.plot.x.y) +
    geom_point(aes(x=pos_cum_qtl, y=pos_cum_gene, color = cis), 
               alpha=1, size=0.7) +
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.x$center, expand = c(0,1)) +
    scale_y_continuous(labels = c(1:19, "X"), breaks = axis.y$center, expand = c(0,1)) +
    theme_classic() +
    theme(legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    scale_color_manual(values = c("cis" = "#13547a", "trans" = "grey50")) +
    ylab("Gene Position") + 
    xlab("eQTL Position")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/eQTL_map.eps", height = 12, width = 12)
eqtl.map
dev.off()

############
# cis-eQTL density, 4462 cis-eQTL totally
breaks <- matrix(c(seq(0, 200, 4), seq(1, 201, 4), seq(2, 202, 4), seq(3, 203, 4)), ncol = 4)
tmp <- as.list(1:ncol(breaks)) 
for(i in 1:ncol(breaks)) {
    tmp[[i]] <- qtl.rna %>%
        filter(lod >= 7.18 & cis == TRUE) %>%
        arrange(qtl_chr, peak_mbp) %>%
        group_by(qtl_chr) %>%
        mutate(win = cut(peak_mbp, breaks = breaks[,i])) %>%
        group_by(qtl_chr, win) %>% 
        summarize(cnt = n()) %>%
        separate(win, into = c("other", "prox", "dist")) %>%
        mutate(prox = as.numeric(prox), 
               dist = as.numeric(dist), 
               mid = 0.5 * (prox + dist)) %>%
        select(qtl_chr, mid, cnt)
}

cis <- bind_rows(tmp[[1]], tmp[[2]], tmp[[3]], tmp[[4]]) %>%
    mutate(qtl_chr = factor(qtl_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X"))) %>%
    na.omit()
rm(tmp)

cis.cum <- cis %>%
    group_by(qtl_chr) %>%
    summarise(chr_len = max(mid)) %>%
    mutate(total = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(cis, ., by=c("qtl_chr" = "qtl_chr")) %>%
    arrange(qtl_chr, mid) %>%
    mutate(pos_cum = mid+total) %>%
    mutate(cnt = cnt -1)

cis.axis <- cis.cum %>%
    group_by(qtl_chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

cis.plot <- ggplot(cis.cum) +
    geom_line(aes(x = pos_cum, y = cnt), color = "#13547a") +
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.x$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(0, max(cis.cum$cnt + 2)), expand = expansion(add=c(0,0))) +
#   geom_hline(aes(yintercept = 100), linetype = 2, color = "grey50") +
    theme_classic() +
    theme(axis.text = element_text(size = 12)) +
    labs(x = "", y = "Number of cis-eQTL")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/eQTL_cis_density.eps", height = 4, width = 12)
cis.plot
dev.off()

# trans-eQTL density, 10855 cis-eQTL totally
breaks <- matrix(c(seq(0, 200, 4), seq(1, 201, 4), seq(2, 202, 4), seq(3, 203, 4)), ncol = 4)
tmp <- as.list(1:ncol(breaks)) 
for(i in 1:ncol(breaks)) {
    tmp[[i]] <- qtl.rna %>%
        filter(lod >= 7.18 & cis != TRUE) %>%
        arrange(qtl_chr, peak_mbp) %>%
        group_by(qtl_chr) %>%
        mutate(win = cut(peak_mbp, breaks = breaks[,i])) %>%
        group_by(qtl_chr, win) %>% 
        summarize(cnt = n()) %>%
        separate(win, into = c("other", "prox", "dist")) %>%
        mutate(prox = as.numeric(prox), 
               dist = as.numeric(dist), 
               mid = 0.5 * (prox + dist)) %>%
        select(qtl_chr, mid, cnt)
}

trans <- bind_rows(tmp[[1]], tmp[[2]], tmp[[3]], tmp[[4]]) %>%
    mutate(qtl_chr = factor(qtl_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X"))) %>%
    na.omit()
rm(tmp)

trans.cum <- trans %>%
    group_by(qtl_chr) %>%
    summarise(chr_len = max(mid)) %>%
    mutate(total = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(trans, ., by=c("qtl_chr" = "qtl_chr")) %>%
    arrange(qtl_chr, mid) %>%
    mutate(pos_cum = mid+total) %>%
    mutate(cnt = cnt -1)

trans.axis <- trans.cum %>%
    group_by(qtl_chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

trans.plot <- ggplot(trans.cum) +
    geom_line(aes(x = pos_cum, y = cnt), color = "grey50") +
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.x$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(0, max(trans.cum$cnt + 2)), expand = expansion(add=c(0,0))) +
    #   geom_hline(aes(yintercept = 100), linetype = 2, color = "grey50") +
    theme_classic() +
    theme(axis.text = element_text(size = 12)) +
    labs(x = "", y = "Number of trans-eQTL")


setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/eQTL_trans_density.eps", height = 4, width = 12)
trans.plot
dev.off()






