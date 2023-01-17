### Purpose: generate figures for DO paper

options(stringsAsFactors = F)

# define colors
coon_blue <- "#2CA7DF"
coon_grey <- "#566977"
coon_purp <- "#955CA5"
coon_red <- "#EC6B63"
coon_turq <- "#63C29C"
coon_yel <- "#FFCB04"

qtl.lipid <- read.table("/Users/rootqz/Desktop/DO/RNA-seq/analysis/eQTL/data/qtl.summary.cecum.vl.csv", sep = ",", header = T, row.names = 1)
qtl.lipid$qtl.chr <- factor(qtl.lipid$qtl.chr, levels = c("1","2","3","4","5","6","7","8","9","10",
                                                          "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))


# load markers, create pseudo markers
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

pseudo.peak$color <- NA

qtl.lipid.g1 <- qtl.lipid %>% filter(id.status != "un-identified")
qtl.lipid.g2 <- qtl.lipid %>% filter(id.status == "un-identified")

### plot lipid QTL with identification
qtl.lipid.g1 <- qtl.lipid.g1 %>%
    mutate(color = case_when(category == "Fatty.Acyl" ~ "#EC6B63",
                             category == "Glycerolipid" ~ "#955CA5",
                             category == "Sphingolipid" ~ "#FFCB04",
                             category == "Phospholipid" ~ "#63C29C"))

qtl.lipid.g1.plot <- qtl.lipid.g1[,c("qtl.chr", "qtl.pos", "qtl.lod", "color")]
names(qtl.lipid.g1.plot) <- names(pseudo.peak)
qtl.lipid.g1.plot <- rbind(qtl.lipid.g1.plot, pseudo.peak)

# calculate cumulative marker position
p1 <- qtl.lipid.g1.plot %>%
    # compute chr size
    group_by(chr) %>%
    summarise(chr_len = max(peak_mbp)) %>%
    
    # calculate cumulative position of each chr
    mutate(total = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # add this infor to extra column
    left_join(qtl.lipid.g1.plot, ., by=c("chr" = "chr")) %>%
    
    # add a cumulative position of each marker
    arrange(chr, peak_mbp) %>%
    mutate(pos_cum = peak_mbp+total)

# make x axis with chr display
axis.df <- p1 %>%
    group_by(chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

p1_out <- ggplot(p1, aes(x=pos_cum, y=lod)) +
    geom_point(aes(color=color), alpha=1, size=0.7) +
    
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(6, 14), expand = expansion(add=c(0,0))) +
    
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "white"),
          legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    scale_color_manual(values = c("#EC6B63" = "#EC6B63", "#955CA5" = "#955CA5",
                                  "#FFCB04" = "#FFCB04", "#63C29C" = "#63C29C")) +
    ylab("LOD") + xlab("Mouse Chromosome")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Manhattan_plot_cecumLipid_QTL_p1.eps", height = 4, width = 16)
p1_out
dev.off()

### plot lipids QTL unknown
qtl.lipid.g2.plot <- qtl.lipid.g2[,c("qtl.chr", "qtl.pos", "qtl.lod")]
qtl.lipid.g2.plot$color <- NA
names(qtl.lipid.g2.plot) <- names(pseudo.peak)
qtl.lipid.g2.plot <- rbind(qtl.lipid.g2.plot, pseudo.peak)

color1 <- c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19")
color2 <- c("2", "4", "6", "8", "10", "12", "14", "16", "18", "X")

for (i in 1:nrow(qtl.lipid.g2.plot)) {
    if (qtl.lipid.g2.plot$chr[i] %in%  color1){
        qtl.lipid.g2.plot$color[i] <- "color1"
    } else {
        qtl.lipid.g2.plot$color[i] <- "color2"
    }
}

# calculate cumulative marker position
p2 <- qtl.lipid.g2.plot %>%
    # compute chr size
    group_by(chr) %>%
    summarise(chr_len = max(peak_mbp)) %>%
    
    # calculate cumulative position of each chr
    mutate(total = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # add this infor to extra column
    left_join(qtl.lipid.g2.plot, ., by=c("chr" = "chr")) %>%
    
    # add a cumulative position of each marker
    arrange(chr, peak_mbp) %>%
    mutate(pos_cum = peak_mbp+total)

# make x axis with chr display
axis.df <- p2 %>%
    group_by(chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

scaleFUN <- function(x) sprintf("%.0f", x)

p2_out <- ggplot(p2, aes(x=pos_cum, y=lod)) +
    geom_point(aes(color=color), alpha=1, size=0.7) +
    
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(6, 50), expand = expansion(add=c(0,0))) +
    
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "white"),
          legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    scale_color_manual(values = c("color1" = "grey", "color2" = "#566977")) +
    scale_y_continuous(labels = scaleFUN, trans = "log10", 
                       breaks = seq(0, 50, by = 10), 
                       limits = c(6, 50), 
                       expand = expand_scale(add=c(0,0))) +
    ylab("LOD") + xlab("Mouse Chromosome")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Manhattan_plot_cecumLipid_QTL_p2.eps", height = 4, width = 16)
p2_out
dev.off()

