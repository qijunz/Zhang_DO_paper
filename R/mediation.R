### Purpose: generate figures for DO paper
###          1. mediation: Genetics -> cecum lipid -> A.mucin
###          2. mediation: Genetics -> bacterial ORFs -> cecum lipid, filtering to reduce image size
###          3. 


# load required libraries
library("qtl2")
library("ggplot2")
library("dplyr")
library("grid")
library("gridExtra")

### mediation: Genetics -> cecum lipid -> A.mucin
mediation_tbl <- readRDS("~/Desktop/DO/QTL/result/mediation/qtl_mediation_OL2AkkGenes.rds")
mediation_tbl <- mediation_tbl %>%
    filter(lod < 100) %>%
    na.omit()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Mediation_QTL_OL2Akk.eps", height = 7, width = 7)
ggplot(data = mediation_tbl, aes(x=lod_drop_perc, y=lod)) + 
    theme_bw() +
    geom_point(alpha=1, size=0.1, color="#F8766D") +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 15)) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    ylab("Orginal A.muciniphila gene QTL LOD score") +
    xlab("LOD drop percentage (%)")  +
    scale_x_continuous(limits = c(-75, 75), expand = c(0,0)) 
dev.off()

# mediation: genetic -> lipid -> microbiome
med_cl2mg <- readRDS("~/Desktop/DO/omics/CHTC_running/mediation_cecumLipid2metagene/mediation_cecumLipid2metagene_LOD_drop.rds")

# UNK_8.373_597.52032_plus
# UNK_8.771_611.539_plus
# UNK_9.183_625.55115_plus

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Mediation_ORF2OL_8.373_597.52032_plus.eps", height = 5, width = 7)
med_cl2mg %>% 
    filter(covar_med == "UNK_8.373_597.52032_plus") %>%
    mutate(color = case_when(genus == "Akkermansia" ~ "Akkermansia",
                             TRUE ~ "no-Akkermansia")) %>%
    ggplot() + 
    geom_point(aes(x = lod_change_percentage, y = lod, color=color), alpha=1, size=0.4) +
    theme_bw() +
    scale_color_manual(values = c("#F8766D", "grey50")) +
    scale_x_continuous(limits = c(-1, 1), expand = c(0,0)) +
    xlab("LOD drop of OL_8.373_597.52032_plus") +
    ylab("")
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Mediation_ORF2OL_8.771_611.539_plus.eps", height = 5, width = 7)
med_cl2mg %>% 
    filter(covar_med == "UNK_8.771_611.539_plus") %>%
    mutate(color = case_when(genus == "Akkermansia" ~ "Akkermansia",
                             TRUE ~ "no-Akkermansia")) %>%
    ggplot() + 
    geom_point(aes(x = lod_change_percentage, y = lod, color=color), alpha=1, size=0.4) +
    theme_bw() +
    scale_color_manual(values = c("#F8766D", "grey50")) +
    scale_x_continuous(limits = c(-1, 1), expand = c(0,0)) +
    xlab("LOD drop of OL_8.771_611.539_plus") +
    ylab("")
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Mediation_ORF2OL_9.183_625.55115_plus.eps", height = 5, width = 7)
med_cl2mg %>% 
    filter(covar_med == "UNK_9.183_625.55115_plus") %>%
    mutate(color = case_when(genus == "Akkermansia" ~ "Akkermansia",
                             TRUE ~ "no-Akkermansia")) %>%
    ggplot() + 
    geom_point(aes(x = lod_change_percentage, y = lod, color=color), alpha=1, size=0.4) +
    theme_bw() +
    scale_color_manual(values = c("#F8766D", "grey50")) +
    scale_x_continuous(limits = c(-1, 1), expand = c(0,0)) +
    xlab("LOD drop of OL_9.183_625.55115_plus") +
    ylab("")
dev.off()

### mediation: Genetics -> bacterial ORFs -> cecum lipid, filtering to reduce image size
mediation_tbl_B2L <- readRDS("~/Desktop/DO/omics/CHTC_running/results_mediation/mediation_metagene2cecumLipid_LOD_drop_tbl.rds")
mediation_tbl_B2L <- mediation_tbl_B2L %>% 
    filter(!is.na(genus)) %>%
    filter(lod.org > 7)

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/Mediation_B2L.eps", height = 5, width = 7)
mediation_tbl_B2L %>% 
    mutate(color = case_when(genus == "Akkermansia" ~ "Akkermansia",
                             TRUE ~ "no-Akkermansia")) %>%
    ggplot() + 
    geom_point(aes(x = lod.drop/lod.org, y = lod.org, color=color), alpha=1, size=0.4) +
    theme_bw() +
    scale_color_manual(values = c("#F8766D", "grey50")) +
    scale_x_continuous(limits = c(-1, 1), expand = c(0,0)) +
    xlab("") +
    ylab("")
dev.off()


