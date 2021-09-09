### Purpose: generate figures for DO paper
###          1. in vitro, bacteria culture lipid data
###          2. in vivo, cecum lipids from mouse colonized with bacteria
###          3. WLS human fecal lipid data

# load required libraries
library("ggplot2")
library("dplyr")
library("tidyr")
library("reshape2")
library("pheatmap")
library("RColorBrewer")
library("grid")
library("gridExtra")

options(stringsAsFactors = F)

# load lipids results from targeted MS
ms.vitro <- read.table("~/Desktop/DO/experiment/LC-MS_lipidomics/metadata.csv", sep = ",", header = T)
# load lipids results from discovery MS
ms.vitro.all <- read.table("~/Desktop/DO/experiment/LC-MS_lipidomics/Lipidomics_AllSignals.csv", sep = ",", header = T, check.names = F)

ms.vitro.ol <- ms.vitro.all %>%
    filter(Lipid_Class == "OL")

# log transformation for OL signal
ms.vitro.ol.logTrans <- ms.vitro.ol
ms.vitro.ol.logTrans[,c(8:36)] <- log10(ms.vitro.ol[,c(8:36)])

ol.tbl <- ms.vitro.ol.logTrans[c(6,8,9),c(5,8:36)] %>%
    tidyr::gather(old, value, -Identification)

rename.tbl <- data.frame(old=unique(ol.tbl$ol))
rename.tbl$new <- c("Blank_bead", "Blank_nobead", rep("Amucin", times=4), 
                    rep("Btheta", times=3), rep("Ecoli", times=3), 
                    rep("B6_Amucin", times=4), rep("B6_AE", times=3), 
                    rep("B6_Btheta", times=4), rep("B6_Ecoli", times=3), rep("B6_GF", times=3))
rename.tbl$group <- c("cecum", rep("culture",times=11), rep("cecum",times=17))

ol.tbl <- merge(ol.tbl, rename.tbl, by = "old")

# subtract the baseline (control)
ol.tbl$value_ctrl[which(ol.tbl$group == "culture")] <- ol.tbl$value[which(ol.tbl$group == "culture")] - mean(ol.tbl$value[which(ol.tbl$group == "culture" & ol.tbl$new == "Blank_nobead")])
ol.tbl$value_ctrl[which(ol.tbl$group == "cecum")] <- ol.tbl$value[which(ol.tbl$group == "cecum")] - mean(ol.tbl$value[which(ol.tbl$group == "cecum" & ol.tbl$new == "Blank_bead")])

# get mean and sd
ol.tbl.stats <- ol.tbl %>%
    group_by(Identification, new) %>%
    summarise(mean=mean(value), sd=sd(value))

ol.tbl.stats <- merge(ol.tbl, ol.tbl.stats, by = c("Identification", "new")) %>%
    select(Identification, new, group, mean, sd) %>%
    unique()

ol.tbl.stats$new <- factor(ol.tbl.stats$new, levels=c("Blank_bead", "B6_Amucin", "B6_Btheta", "B6_Ecoli", "B6_AE", "B6_GF", "Blank_nobead", "Amucin", "Btheta", "Ecoli"))

# subtract the baseline (control)
for (ol in c("OL a-30:0", "OL a-31:0", "OL a-32:0")) {
    ol.tbl.stats$mean[which(ol.tbl.stats$group == "culture" & ol.tbl.stats$Identification == ol)] <- 
        ol.tbl.stats$mean[which(ol.tbl.stats$group == "culture"& ol.tbl.stats$Identification == ol)] - 
        ol.tbl.stats$mean[which(ol.tbl.stats$new == "Blank_nobead" & ol.tbl.stats$Identification == ol)]
    
    ol.tbl.stats$mean[which(ol.tbl.stats$group == "cecum" & ol.tbl.stats$Identification == ol)] <- 
        ol.tbl.stats$mean[which(ol.tbl.stats$group == "cecum"& ol.tbl.stats$Identification == ol)] - 
        ol.tbl.stats$mean[which(ol.tbl.stats$new == "Blank_bead" & ol.tbl.stats$Identification == ol)]
}

# boxplot using mean and sd to show error bar
ggplot(ol.tbl.stats) + 
    geom_boxplot(aes(new, mean), fill="grey70", width = 0.8) + 
    theme_bw() +
    geom_errorbar(aes(new, ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
    facet_grid(Identification ~ group, scales = "free", space = "free") +
    theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1),
          axis.title = element_text(size = 15),
          strip.text = element_text(size = 12)) +
    labs(x = "", y= "log10 LC-MS signal")


ol.tbl.plot <- ol.tbl %>%
    left_join(ol.tbl.stats)

# in vitro culture
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/OLs_boxplot_culture.eps", height = 9, width = 4)
ol.tbl.plot %>%
    filter(group == "culture", new != "Blank_nobead") %>%
    ggplot() + 
    geom_boxplot(aes(new, value_ctrl), fill="grey70", width = 0.8) + 
    theme_bw() +
    scale_y_continuous(limits = c(0,6)) +
    facet_grid(Identification ~ group, scales = "free_x", space = "free") + 
        geom_errorbar(aes(new, ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
        theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1),
              axis.title = element_text(size = 15),
              strip.text = element_text(size = 12)) +
        labs(x = "", y= "log10 LC-MS signal")
dev.off()


# in vivo colonized mice cecum
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/OLs_boxplot_cecum.eps", height = 9, width = 4)
ol.tbl.plot %>%
    filter(group == "cecum", new != "Blank_bead") %>%
    ggplot() + 
    geom_boxplot(aes(new, value_ctrl), fill="grey70", width = 0.8) + 
    theme_bw() +
    scale_y_continuous(limits = c(0,6)) +
    facet_grid(Identification ~ group, scales = "free_x", space = "free") + 
    geom_errorbar(aes(new, ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
    theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1),
          axis.title = element_text(size = 15),
          strip.text = element_text(size = 12)) +
    labs(x = "", y= "log10 LC-MS signal")
dev.off()

############### WLS Human fecal sample lipid data ############### 
ms.human <- read.table("~/Desktop/DO/experiment/LC-MS_lipidomics_2020-06-20/Final_Results_humanfecal.csv", sep = ",", header = T)
wls.meta <- read.table("~/Desktop/DO/experiment/LC-MS_lipidomics_2020-06-20/WLS_HumanFeca_sample_forOL_detection.csv", sep = ",", header = T)

# rename sample id
wls.meta <- wls.meta %>%
    separate(col = sample.id, sep = "_", remove = F, into = c(NA, "ms.number"))
wls.meta$ms.id <- paste0("WLS", wls.meta$ms.number)

# OL a-30:0_8.45_597.52051(+)
# OL a-31:0_8.785_611.53571(+)
# OL a-32:0_9.224_625.55109(+)

wls.plot <- ms.human %>%
    filter(Lipid.Class == "OL") %>%
    filter(Identifier %in% c("OL a-30:0_8.45_597.52051(+)", "OL a-31:0_8.785_611.53571(+)", "OL a-32:0_9.224_625.55109(+)")) %>%
    select(c("Identifier", wls.meta$ms.id)) %>% 
    melt()

wls.meta$sample.id <- gsub("_", "", wls.meta$sample.id)

names(wls.plot) <- c("OL", "WLS", "MS")
wls.plot <- wls.plot %>%
    left_join(wls.meta[,c("sample.id", "akk.otu")], by = c("WLS" = "sample.id"))

wls.plot$MS <- log10(wls.plot$MS)

cor.test(wls.plot$MS[which(wls.plot$OL == "OL a-30:0_8.45_597.52051(+)")], 
         wls.plot$akk.otu[which(wls.plot$OL == "OL a-30:0_8.45_597.52051(+)")])

cor.test(wls.plot$MS[which(wls.plot$OL == "OL a-31:0_8.785_611.53571(+)")], 
         wls.plot$akk.otu[which(wls.plot$OL == "OL a-31:0_8.785_611.53571(+)")])

cor.test(wls.plot$MS[which(wls.plot$OL == "OL a-32:0_9.224_625.55109(+)")], 
         wls.plot$akk.otu[which(wls.plot$OL == "OL a-32:0_9.224_625.55109(+)")])

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/OLs_WLS_human_OL30.eps", height = 3, width = 4)
wls.plot %>%
    filter(OL == "OL a-30:0_8.45_597.52051(+)") %>%
    ggplot(aes(x = akk.otu, y = MS)) +
    geom_point() +
    theme_classic() +
    geom_smooth(method='lm') +
    xlab("Gut A. muciniphila abundance") +
    ylab("LC-MS/MS signal log10")
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/OLs_WLS_human_OL31.eps", height = 3, width = 4)
wls.plot %>%
    filter(OL == "OL a-31:0_8.785_611.53571(+)") %>%
    ggplot(aes(x = akk.otu, y = MS)) +
    geom_point() +
    theme_classic() +
    geom_smooth(method='lm') +
    xlab("Gut A. muciniphila abundance") +
    ylab("LC-MS/MS signal log10")
dev.off()

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/OLs_WLS_human_OL32.eps", height = 3, width = 4)
wls.plot %>%
    filter(OL == "OL a-32:0_9.224_625.55109(+)") %>%
    ggplot(aes(x = akk.otu, y = MS)) +
    geom_point() +
    theme_classic() +
    geom_smooth(method='lm') +
    xlab("Gut A. muciniphila abundance") +
    ylab("LC-MS/MS signal log10")
dev.off()

### statistical test
# Differences between groups were evaluated using unpaired two-tailed Welch’s t-test
ols <- c("OL a-30:0", "OL a-31:0", "OL a-32:0")
for (i in 1:3) {
    t.test(ol.tbl$value[which(ol.tbl$Identification == ols[i] & ol.tbl$group == "culture" & ol.tbl$new == "Amucin")],
           ol.tbl$value[which(ol.tbl$Identification == ols[i] & ol.tbl$group == "culture" & ol.tbl$new == "Btheta")],
           alternative = "two.sided", paired = F)
}



