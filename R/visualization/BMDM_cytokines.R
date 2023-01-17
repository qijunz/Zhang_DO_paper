### Purpose: generate figures for DO paper

# load required libraries
options(stringsAsFactors = FALSE)
setwd("~/Desktop/DO/experiment/BMDM/")

tnf_b6 <- read.table("ELISA/data/B6_TNF-alpha_20210527.tsv", sep = "\t", header = T)
tnf_129 <- read.table("ELISA/data/129_TNF-alpha_20210528.tsv", sep = "\t", header = T)
il6_b6 <- read.table("ELISA/data/B6_IL6_20210527.tsv", sep = "\t", header = T)
il6_129 <- read.table("ELISA/data/129_IL6_20210528.tsv", sep = "\t", header = T)

# add mouse strain
tnf_b6$mouse <- "B6"
il6_b6$mouse <- "B6"
tnf_129$mouse <- "129"
il6_129$mouse <- "129"

### co-culture group
baxplot <- function(df, mouse, mouse_color, cytokine) {
    plot_out <- df %>%
        filter(group == "OL+LPS") %>%
        mutate(conc = factor(ifelse(sample_ID %in% c("LPS-c1", "LPS-c2", "LPS-c3"), "0", OL_conc_ng_per_mL), 
                             levels = c("0", "100", "500", "1000", "5000", "10000"))) %>%
        ggplot() +
        geom_boxplot(aes(x = conc, y = concentration_pg_per_mL, fill = mouse)) +
        scale_fill_manual(values = mouse_color) + 
        theme_classic() + 
        theme(axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text = element_text(size = 10)) +
        xlab("OL conc. ng/ml (LPS = 10ng/mL)") +
        ylab(paste0(cytokine, " concentration (pg/mL)"))
    
    return(plot_out)
}

### plot B6
plot_tnf_b6 <- baxplot(tnf_b6, "B6", CCcolors[2], "TNF-alpha")
plot_il6_b6 <- baxplot(il6_b6, "B6", CCcolors[2], "IL6")

plot_tnf_129 <- baxplot(tnf_129, "129", CCcolors[3], "TNF-alpha")
plot_il6_129 <- baxplot(il6_129, "129", CCcolors[3], "IL6")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/BMDM_B6_129_coculture.eps", height = 6, width = 9)
grid.arrange(plot_tnf_b6, plot_tnf_129, plot_il6_b6, plot_il6_129,  nrow=2, ncol=2)
dev.off()

### LPS group ###
lps_baxplot <- function(df, mouse, mouse_color, cytokine) {
    plot_out <- df %>%
        filter(group == "LPS") %>%
        ggplot() +
        geom_boxplot(aes(x = as.factor(LPS_conc_ng_per_mL), y = concentration_pg_per_mL, fill = mouse)) +
        scale_fill_manual(values = mouse_color) + 
        #        scale_y_continuous(limits = c(0,41000)) +
        theme_classic() + 
        theme(axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text = element_text(size = 10)) +
        xlab("LPS conc. (ng/mL)") +
        ylab(paste0(cytokine, " concentration (pg/mL)"))
    
    return(plot_out)
}
lps_tnf_b6 <- lps_baxplot(tnf_b6, "B6", CCcolors[2], "TNF-alpha")
lps_il6_b6 <- lps_baxplot(il6_b6, "B6", CCcolors[2], "IL6")

lps_tnf_129 <- lps_baxplot(tnf_129, "129", CCcolors[3], "TNF-alpha")
lps_il6_129 <- lps_baxplot(il6_129, "129", CCcolors[3], "IL6")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/BMDM_B6_129_LPS.eps", height = 6, width = 9)
grid.arrange(lps_tnf_b6, lps_tnf_129, lps_il6_b6, lps_il6_129,  nrow=2, ncol=2)
dev.off()

### OL group
ol_baxplot <- function(df, mouse, mouse_color, cytokine) {
    plot_out <- df %>%
        filter(group == "OL") %>%
        ggplot() +
        geom_boxplot(aes(x = as.factor(OL_conc_ng_per_mL), y = concentration_pg_per_mL, fill = mouse)) +
        scale_fill_manual(values = mouse_color) + 
        theme_classic() + 
        theme(axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text = element_text(size = 10)) +
        xlab("OL conc. (ng/mL)") +
        ylab(paste0(cytokine, " concentration (pg/mL)"))
    
    return(plot_out)
}
ol_tnf_b6 <- ol_baxplot(tnf_b6, "B6", CCcolors[2], "TNF-alpha")
ol_il6_b6 <- ol_baxplot(il6_b6, "B6", CCcolors[2], "IL6")

ol_tnf_129 <- ol_baxplot(tnf_129, "129", CCcolors[3], "TNF-alpha")
ol_il6_129 <- ol_baxplot(il6_129, "129", CCcolors[3], "IL6")

setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/BMDM_B6_129_OL.eps", height = 6, width = 9)
grid.arrange(ol_tnf_b6, ol_tnf_129, ol_il6_b6, ol_il6_129,  nrow=2, ncol=2)
dev.off()

### cytokines from Quansys w/ extra samples
cyto_all <- read.table("Quansys/quansys_result.csv", sep = ",", header = T, check.names = F)
cyto_all$mouse <- factor(cyto_all$mouse, levels =  c("B6", "129"))
cyto_all$OL_conc_ug_ml <- factor(cyto_all$OL_conc_ug_ml, levels = c("no-LPS", "0", "0.5", "1"))

cytos <- names(cyto_all)[c(8:9,14:16,18:23)]
out <- list()

for (i in 1:length(cytos)) {
    cyto <- cytos[i]
    plot_df <- cyto_all %>%
        filter(OL_conc_ug_ml != 5) %>%
        select(c("OL_conc_ug_ml", "mouse", "new_ID", cyto))
    names(plot_df)[4] <- "cyto"
    
    out[[i]] <- plot_df %>%
        ggplot() +
        geom_boxplot(aes(x = as.factor(OL_conc_ug_ml), y = cyto, color = mouse)) +
        facet_grid(~ mouse, scales = "free_x", space ="free_x") +
        theme_bw() +
        scale_color_manual(values = c(CCcolors[2], CCcolors[3])) +
        ggtitle(cyto) +
        xlab("OL conc. ug/mL") +
        ylab(paste0(cyto, " conc. pg/mL")) +
        theme(axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              plot.title = element_text(size = 16, hjust = 0.5),
              strip.background = element_rect(fill = c(CCcolors[2], CCcolors[3])),
              strip.text = element_text(size=12, color="white"))
    
}


setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/BMDM_B6_129_quansys_wExtraSample.eps", height = 9, width = 16)
grid.arrange(out[[1]], out[[2]], out[[3]], out[[4]], 
             out[[5]], out[[6]], out[[7]], out[[8]], 
             out[[9]], out[[10]], out[[11]],
             nrow=3, ncol=4)
dev.off()

### 2021-07-28
cyto_all <- read.table("Quansys/quansys_result.csv", sep = ",", header = T, check.names = F)
cyto_all$mouse <- factor(cyto_all$mouse, levels =  c("B6", "129"))
cyto_all$OL_conc_ug_ml <- factor(cyto_all$OL_conc_ug_ml, levels = c("no-LPS", "0", "0.5", "1"))

cytos <- names(cyto_all)[c(9,14:16,18:23)]

for (i in 1:length(cytos)) {
    cyto <- cytos[i]
    plot_df <- cyto_all %>%
        filter(OL_conc_ug_ml != 5) %>%
        select(c("OL_conc_ug_ml", "mouse", "new_ID", cyto))
    names(plot_df)[4] <- "cyto"
    
    setEPS()
    postscript(paste0("~/Desktop/ReyLab/paper/DO_metagenomic/figure/BMDM_B6_129_quansys_", cyto,".eps"), height = 4, width = 4)
    plot_df %>%
        ggplot() +
        geom_boxplot(aes(x = as.factor(OL_conc_ug_ml), y = cyto, color = mouse)) +
        facet_grid(~ mouse, scales = "free_x", space ="free_x") +
        theme_classic() +
        scale_color_manual(values = c(CCcolors[2], CCcolors[3])) +
        ggtitle(cyto) +
        xlab("OL conc. ug/mL") +
        ylab(paste0(cyto, " conc. pg/mL")) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              legend.position="none",
              plot.title = element_text(size = 16, hjust = 0.5),
              strip.background = element_rect(fill = c(CCcolors[2], CCcolors[3])),
              strip.text = element_text(size=12, color="white"))
    dev.off()
    
}

# heatmap
cyto_all <- na.omit(cyto_all)
cyto_all$group <- paste0(cyto_all$mouse, "_", cyto_all$treatment, "_", cyto_all$OL_conc_ug_ml)

cyto_ave <- cyto_all %>%
    group_by(group) %>%
    mutate(`mIL-1b_mean` = mean(`mIL-1b`),
           `mIL-6_mean` = mean(`mIL-6`),
           `mIL-10_mean` = mean(`mIL-10`),
           `mIL-12_mean` = mean(`mIL-12`),
           `mMCP-1_mean` = mean(`mMCP-1`),
           `mTNFa_mean` = mean(`mTNFa`),
           `mMIP-1a_mean` = mean(`mMIP-1a`),
           `mGM-CSF_mean` = mean(`mGM-CSF`),
           `mRANTES_mean` = mean(`mRANTES`))

cyto_ave_plot <- cyto_ave[,c(24:33)] %>% unique() %>% as.data.frame()
rownames(cyto_ave_plot) <- cyto_ave_plot$group
cyto_ave_plot$group <- NULL

cyto_ave_plot <- as.data.frame(t(cyto_ave_plot))
cyto_ave_plot <- cyto_ave_plot[,c(1,3,2,4,5,7,6,8)]

pheatmap(cyto_ave_plot,
         cluster_cols = F,
         scale = "row")

### 2021-08-16
tnf_b6 <- read.table("ELISA/data/B6_TNF-alpha_20210527.tsv", sep = "\t", header = T)
tnf_129 <- read.table("ELISA/data/129_TNF-alpha_20210528.tsv", sep = "\t", header = T)
il6_b6 <- read.table("ELISA/data/B6_IL6_20210527.tsv", sep = "\t", header = T)
il6_129 <- read.table("ELISA/data/129_IL6_20210528.tsv", sep = "\t", header = T)

# add mouse strain
tnf_b6$mouse <- "B6"
il6_b6$mouse <- "B6"
tnf_129$mouse <- "129"
il6_129$mouse <- "129"

tnf <- rbind(tnf_b6, tnf_129)
il6 <- rbind(il6_b6, il6_129)

tnf$mouse <- factor(tnf$mouse, levels = c("B6", "129"))
il6$mouse <- factor(il6$mouse, levels = c("B6", "129"))

# LPS, tnf
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/BMDM_ELISA_tnf_LPS.eps", height = 5, width = 7.5)
tnf %>%
    filter(group == "LPS") %>%
    ggplot() +
    geom_boxplot(aes(x = as.factor(LPS_conc_ng_per_mL), y = concentration_pg_per_mL, fill = mouse)) +
    scale_fill_manual(values = c("#888888", "#F08080")) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 15),
          legend.position="none") +
    scale_y_continuous(limits = c(0,32000)) +
    facet_grid(~ mouse, scales = "free_y", space ="free_y") +
    xlab("") +
    ylab("")
dev.off()

# OL, tnf
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/BMDM_ELISA_tnf_OL.eps", height = 5, width = 7.5)
tnf %>%
    filter(group == "OL") %>%
    ggplot() +
    geom_boxplot(aes(x = as.factor(OL_conc_ng_per_mL), y = concentration_pg_per_mL, fill = mouse)) +
    scale_fill_manual(values = c("#888888", "#F08080")) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 15),
          legend.position="none") +
    scale_y_continuous(limits = c(0,32000)) +
    facet_grid(~ mouse, scales = "free_y", space ="free_y") +
    xlab("") +
    ylab("")
dev.off()

# LPS, il6
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/BMDM_ELISA_il6_LPS.eps", height = 5, width = 7.5)
il6 %>%
    filter(group == "LPS") %>%
    ggplot() +
    geom_boxplot(aes(x = as.factor(LPS_conc_ng_per_mL), y = concentration_pg_per_mL, fill = mouse)) +
    scale_fill_manual(values = c("#888888", "#F08080")) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 15),
          legend.position="none") +
    scale_y_continuous(limits = c(0,42000)) +
    facet_grid(~ mouse, scales = "free_y", space ="free_y") +
    xlab("") +
    ylab("")
dev.off()

# OL, il6
setEPS()
postscript("~/Desktop/ReyLab/paper/DO_metagenomic/figure/BMDM_ELISA_il6_OL.eps", height = 5, width = 7.5)
il6 %>%
    filter(group == "OL") %>%
    ggplot() +
    geom_boxplot(aes(x = as.factor(OL_conc_ng_per_mL), y = concentration_pg_per_mL, fill = mouse)) +
    scale_fill_manual(values = c("#888888", "#F08080")) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 15),
          legend.position="none") +
    scale_y_continuous(limits = c(0,42000)) +
    facet_grid(~ mouse, scales = "free_y", space ="free_y") +
    xlab("") +
    ylab("")
dev.off()
