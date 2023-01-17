### Purpose: plot figures of A.muciniphila and OL QTL

# load required libraries
options(stringsAsFactors = FALSE)
source("/Users/rootqz/R/QZ_functions.R")

markers <- readRDS("/Users/rootqz/Desktop/DO/metagenomics_QTL/data/data_raw/DO_69k_markers.rds")

# load mediation result table
data.dir <- "/Users/rootqz/Desktop/DO/omics/CHTC_running/results_mediation/"

cecal_l <- readRDS(paste0(data.dir, "mediation_metagene2cecumLipid_LOD_drop_tbl.rds"))
cecal_l$lod.drop.prec <- cecal_l$lod.drop/cecal_l$lod.org*100

# plot LOD_drop vs. LOD_ori
# plot original lod vs. add metagenes as cov lod
cecal_l$group <- cecal_l$genus
cecal_l$group[which(is.na(cecal_l$genus))] <- "others"
cecal_l$group[which(cecal_l$genus != "Akkermansia")] <- "others"
cecal_l$lod.drop.prec <- cecal_l$lod.drop/cecal_l$lod.org*100

cecal_l_plot <- cecal_l %>% filter(lod.org > 7.2)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/Mediation_QTL_CecalLipid.eps", height = 7, width = 9)
ggplot(data = cecal_l_plot, aes(x=lod.drop.prec, y=lod.org, color=group)) + 
    theme_bw() +
    geom_point(alpha=1, size=0.1) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 15)) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    ylab("Orginal cecal lipid QTL LOD score") +
    xlab("LOD drop percentage (%)")  +
    ggtitle("Cecal Lipids mediation QTL") +
    scale_x_continuous(limits = c(-75, 75), expand = c(0,0)) +
    scale_color_manual(values = c("Akkermansia" = "#F8766D", "others" = "grey70"), name = "Metagene annotation")
dev.off()

#### Manhatton plot for A.muciniphila MAGs and OL
# load phenotype matrix
mag_phe <- readRDS("~/Desktop/DO/omics/integration/project/QTL_MAG/data/A_mucin_MAGs_Cov_libSizeNorm_rankZ_264mice.rds")
mag_qtl_out <- readRDS("~/Desktop/DO/omics/integration/project/QTL_MAG/A_mucin_221mouse/qtl_scan1_out_221mice.rds")
mag_qtl_peak <- readRDS("~/Desktop/DO/omics/integration/project/QTL_MAG/A_mucin_221mouse/qtl_blup_tbl_minLOD5_withMarker_221mice.rds")

# load cecal lipid matrix
cecal_lipid <- readRDS("~/Desktop/DO/metabolomics/data_raw/attie_cecum_lipids_normalized.rds")

ol_qtl_out <- readRDS("~/Desktop/DO/metabolomics/OLs/cecum_ols_221mice_2020-07-16/qtl_scan1_out_221mice.rds")
ol_qtl_peak <- readRDS("~/Desktop/DO/metabolomics/OLs/cecum_ols_221mice_2020-07-16/qtl_blup_tbl_minLOD5_withMarker_221mice.rds")


mh_plot <- function(phe, phe_qtl_out, xmax, lodcut) {
    # inputs are:
    #    1. phenotype name
    #    2. the lod score data frame in 69k marker position
    #    3. x-axis maximal limit
    #    4. cutoff of lod for highlight
    
    phe_qtl <- as.data.frame(phe_qtl_out[,phe,drop=F])
    names(phe_qtl) <- "lod"
    phe_qtl$marker <- rownames(phe_qtl)
    phe_qtl_plot <- phe_qtl %>% left_join(markers, by = c("marker" = "marker"))
    
    color1 <- c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19")
    color2 <- c("2", "4", "6", "8", "10", "12", "14", "16", "18", "X")
    
    qtl.plot <- phe_qtl_plot[,c("chr", "pos", "lod")]
    qtl.plot$color <- NA
    
    for (i in 1:nrow(qtl.plot)) {
        if (qtl.plot$chr[i] %in%  color1){
            qtl.plot$color[i] <- "color1"
        } else {
            qtl.plot$color[i] <- "color2"
        }
    }
    
    qtl.plot$color[which(qtl.plot$lod > lodcut)] <- "color3" 
    
    qtl.plot$chr <- factor(qtl.plot$chr, levels = c("1","2","3","4","5","6","7","8","9","10",
                                                    "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))
    
    # calculate cumulative marker position
    mh.plot <- qtl.plot %>%
        # compute chr size
        group_by(chr) %>%
        summarise(chr_len = max(pos)) %>%
        
        # calculate cumulative position of each chr
        mutate(total = cumsum(chr_len) - chr_len) %>%
        select(-chr_len) %>%
        
        # add this infor to extra column
        left_join(qtl.plot, ., by=c("chr" = "chr")) %>%
        
        # add a cumulative position of each marker
        arrange(chr, pos) %>%
        mutate(pos_cum = pos+total)
    
    # make x axis with chr display
    axis.df <- mh.plot %>%
        group_by(chr) %>%
        summarise(center=(max(pos_cum) + min(pos_cum))/2)
    
    out.plot <- ggplot(mh.plot, aes(x=pos_cum, y=lod)) +
        geom_point(aes(color=color), alpha=1, size=0.5) +
        
        # custom x,y axis
        scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
        scale_y_continuous(limits = c(0,xmax), expand = expansion(add=c(0,0))) +
        
        theme_classic() +
        theme(legend.position="none") +
        scale_color_manual(values = c("color1" = "#13547a", "color2" = "#80d0c7", "color3" = "#f05454")) +
        ylab("LOD") + xlab("Mouse Chromosome") +
        ggtitle(phe)
    
    return(out.plot)
}

# plot OL
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/UNK_8.373_597.52032_plus_mhplot.eps", height = 2.5, width = 12)
mh_plot("UNK_8.373_597.52032_plus", ol_qtl_out,10,6)
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/UNK_8.771_611.539_plus_mhplot.eps", height = 2.5, width = 12)
mh_plot("UNK_8.771_611.539_plus", ol_qtl_out,10,6)
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/UNK_9.183_625.55115_plus_mhplot.eps", height = 2.5, width = 12)
mh_plot("UNK_9.183_625.55115_plus", ol_qtl_out,10,6)
dev.off()

# plot Akkermansia MAGs
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/chr1_DO097_DBSCAN_round3_1_mhplot.eps", height = 2.5, width = 12)
mh_plot("DO097_DBSCAN_round3_1", mag_qtl_out,8,5.5)
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/chr7_DO048_DBSCAN_round1_3_mhplot.eps", height = 2.5, width = 12)
mh_plot("DO048_DBSCAN_round1_3", mag_qtl_out,8,5.4)
dev.off()


### function to plot founder allele effects, give 8 values
#   Inputs are 8 founder allele effects vector with A-H order (can have in their names)
CCcolor <- CCcolors
allele_eff_plot <- function(eff) {
    plot.df <- as.data.frame(matrix(nrow = 8, ncol = 0))
    plot.df$eff <- as.numeric(eff[,c("A", "B", "C", "D", "E", "F", "G", "H")])
    plot.df$strain <- factor(names(CCcolors), level=names(CCcolors))
    ylim <- max(abs(plot.df$eff)) + 0.2
    p <- ggplot(plot.df) + 
        geom_point(aes(x=strain, y=eff, color=strain), size=4) +
        geom_hline(yintercept=0) +
        theme_classic() +
        theme(axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 11)) +
        scale_color_manual(values =  CCcolor) +
        scale_y_continuous(limits = c(-ylim, ylim)) +
        xlab("") +
        ylab("")
    
    return(p)
}

akk_qtl_out <- readRDS("~/Desktop/DO/integration/project/QTL_MAG/A_mucin_221mouse/qtl_scan1_out_221mice.rds")
akk_qtl_peak <- readRDS("~/Desktop/DO/integration/project/QTL_MAG/A_mucin_221mouse/qtl_blup_tbl_minLOD5_withMarker_221mice.rds")

chr1_akk_phe <- "DO097_DBSCAN_round3_1"
chr2_akk_phe <- "DO050_DBSCAN_round1_3"
chr7_akk_phe <- "DO048_DBSCAN_round1_3"
chr12_akk_phe <- "DO025_DBSCAN_round1_2"
chr15_akk_phe <- "DO374_DBSCAN_round2_3"


chr1_ol_phe <- "UNK_8.771_611.539_plus"
chr2_ol_phe <- "UNK_8.373_597.52032_plus"
chr7_ol_phe <- "UNK_8.771_611.539_plus"
chr12_ol_phe <- "UNK_8.373_597.52032_plus"
chr15_ol_phe <- "UNK_8.373_597.52032_plus"

# Akk MGAs 
p_akkchr1 <- akk_qtl_peak %>% filter(pheno == chr1_akk_phe,  chr == 1) %>% allele_eff_plot()
p_akkchr2 <- akk_qtl_peak %>% filter(pheno == chr2_akk_phe,  chr == 2) %>% allele_eff_plot()
p_akkchr7 <- akk_qtl_peak %>% filter(pheno == chr7_akk_phe,  chr == 7) %>% allele_eff_plot()
p_akkchr12 <- akk_qtl_peak %>% filter(pheno == chr12_akk_phe,  chr == 12) %>% allele_eff_plot()
p_akkchr15 <- akk_qtl_peak %>% filter(pheno == chr15_akk_phe,  chr == 15) %>% allele_eff_plot()

# cecal OL
p_olchr1 <- ol_qtl_peak %>% filter(pheno == chr1_ol_phe,  chr == 1) %>% allele_eff_plot()
p_olchr2 <- ol_qtl_peak %>% filter(pheno == chr2_ol_phe,  chr == 2) %>% allele_eff_plot()
p_olchr7 <- ol_qtl_peak %>% filter(pheno == chr7_ol_phe,  chr == 7) %>% allele_eff_plot()
p_olchr12 <- ol_qtl_peak %>% filter(pheno == chr12_ol_phe,  chr == 12) %>% allele_eff_plot()
p_olchr15 <- ol_qtl_peak %>% filter(pheno == chr15_ol_phe,  chr == 15) %>% allele_eff_plot()


# put together
G_p_akkchr1 <- ggplotGrob(p_akkchr1)
G_p_akkchr2 <- ggplotGrob(p_akkchr2)
G_p_akkchr7 <- ggplotGrob(p_akkchr7)
G_p_akkchr12 <- ggplotGrob(p_akkchr12)
G_p_akkchr15 <- ggplotGrob(p_akkchr15)

G_p_olchr1 <- ggplotGrob(p_olchr1)
G_p_olchr2 <- ggplotGrob(p_olchr2)
G_p_olchr7 <- ggplotGrob(p_olchr7)
G_p_olchr12 <- ggplotGrob(p_olchr12)
G_p_olchr15 <- ggplotGrob(p_olchr15)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/Amucin_OL_QTL_effects.eps", height = 10, width = 7)
grid::grid.newpage()
grid::grid.draw(cbind(rbind(G_p_akkchr1, G_p_akkchr2, G_p_akkchr7, G_p_akkchr12, G_p_akkchr15),
                      rbind(G_p_olchr1, G_p_olchr2, G_p_olchr7, G_p_olchr12, G_p_olchr15)))
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figures/Amucin_OL_QTL_effects_legend.eps", height = 2, width = 5)
p_akkchr1
dev.off()
