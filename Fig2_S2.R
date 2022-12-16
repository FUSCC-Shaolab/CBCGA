rm(list = ls())

#library(dplyr)
#library(stringr)
#library(tidyr)
#library(data.table)
#library(survtype)
#library(ggplot2)
library(ggridges)
#library(escape)
library(maftools)
library(trackViewer)
library(export)
#library(BradleyTerry2)
library(tidyverse)


# Data preparation --------------------------------------------------------
load("/data/CBCGA.Extended_MergedData_V2.5_220722.Rdata")
CBCGAClin <- CBCGA_Cohort.Info
#rm(list = ls()[-c(20:21)])
rm(list = ls()[! ls() %in% c("CBCGAClin","CBCGA_WES_Somatic")])
load("/data/Fig2_S2_data.RData")

CBCGAClin <- janitor::clean_names(CBCGAClin)
colnames(CBCGAClin)[1] <- "Tumor_Sample_Barcode"
CBCGAClin <- subset(CBCGAClin, histological_type == "IDC")
CBCGAClin_WES <- subset(CBCGAClin, exome_sequencing_paired == "Yes")

TCGAClin$Tumor_Sample_Barcode <- str_sub(TCGAClin$sampleID, 1, 12)
TCGAClin_white <- subset(TCGAClin, Race.Category == "WHITE" & histological_type == "Infiltrating Ductal Carcinoma")
TCGAClin_WES <- subset(TCGAClin, Tumor_Sample_Barcode %in% TCGAMaf_df$Tumor_Sample_Barcode)
TCGAClin_WES_white <- subset(TCGAClin_white, Tumor_Sample_Barcode %in% TCGAMaf_df$Tumor_Sample_Barcode)

# Figure 2A ---------------------------------------------------------------
CBCGAvsTCGA_age <- data.frame(
  age = c(CBCGAClin_WES$age_at_surgery, CBCGAClin_WES$age_at_surgery, TCGAClin_WES_white$age_at_initial_pathologic_diagnosis, TCGAClin_WES_white$age_at_initial_pathologic_diagnosis),
  cohort = c(rep("CBCGA", nrow(CBCGAClin_WES)), rep("CBCGA", nrow(CBCGAClin_WES)), rep("TCGA", nrow(TCGAClin_WES_white)), rep("TCGA", nrow(TCGAClin_WES_white))),
  pam50 = c(rep("All", nrow(CBCGAClin_WES)), CBCGAClin_WES$intrinsic_subtype_pam50, rep("All", nrow(TCGAClin_WES_white)), TCGAClin_WES_white$PAM50Call_RNAseq)
)
CBCGAvsTCGA_age$pam50 <- ifelse(CBCGAvsTCGA_age$pam50 == "" | is.na(CBCGAvsTCGA_age$pam50) == T, "Unknown", CBCGAvsTCGA_age$pam50)
CBCGAvsTCGA_age$pam50 <- paste(CBCGAvsTCGA_age$pam50, CBCGAvsTCGA_age$cohort, sep = "_")
CBCGAvsTCGA_age <- subset(CBCGAvsTCGA_age, !pam50 %in% c("Unknown_CBCGA", "Unknown_TCGA"))
CBCGAvsTCGA_age$pam50 <- factor(CBCGAvsTCGA_age$pam50, levels = rev(c(
  "All_CBCGA", "All_TCGA", "LumA_CBCGA", "LumA_TCGA", "LumB_CBCGA", "LumB_TCGA",
  "Her2_CBCGA", "Her2_TCGA", "Basal_CBCGA", "Basal_TCGA", "Normal_CBCGA", "Normal_TCGA"
)))

p = ggplot(CBCGAvsTCGA_age, aes(x = age, y = pam50, fill = cohort)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01) +
  scale_fill_manual(values = c("CBCGA" = "#0077C9", "TCGA" = "#91B720")) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(
    x = guide_axis(title = "Age"),
    fill = guide_legend(title = "")
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(colour = "grey"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "bottom"
  )
graph2pdf(p,file = '/results/Fig2A.pdf', width = 4.2, height = 6.9)

aggregate(age ~ pam50, CBCGAvsTCGA_age, median)
pval <- NULL
for (i in c("All", "LumA", "LumB", "Her2", "Basal", "Normal")) {
  tmp_function <- paste("age", "pam50", sep = "~") %>% as.formula()
  tmp_df <- subset(CBCGAvsTCGA_age, pam50 %in% paste(i, c("CBCGA", "TCGA"), sep = "_"))
  tmp_wcox <- wilcox.test(tmp_function, tmp_df)
  tmp_pval <- tmp_wcox$p.value
  pval <- c(pval, tmp_pval)
}
fdr <- p.adjust(as.numeric(pval), method = "fdr")
fdr


# Figure 2B ---------------------------------------------------------------
genes <- c(
  "TP53", "PIK3CA", "GATA3", "MAP3K1", "KMT2C", "AKT1", "PTEN", "NF1",
  "ARID1A", "SF3B1", "FAT4", "LRP1B", "MAP2K4"
)
mycol <- c(
  "Splice_Site" = "#6b4694", "Missense_Mutation" = "#4d93cd", "Nonsense_Mutation" = "#e63a3a",
  "Frame_Shift_Ins" = "#fab71b", "Frame_Shift_Del" = "#fab71b", "In_Frame_Ins" = "#ecd71e",
  "In_Frame_Del" = "#ecd71e", "Nonstop_Mutation" = "#ab5b9e", "Translation_Start_Site" = "#018b38",
  "Multi_Hit" = "#cecece"
)
vc_nonSyn <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del",
  "Splice_Site", "Splice_Region", "Nonstop_Mutation", "Translation_Start_Site"
)


CBCGAMaf_df <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin$Tumor_Sample_Barcode)
CBCGAMaf <- read.maf(CBCGAMaf_df)
TCGAMaf_df_white <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white$Tumor_Sample_Barcode)
TCGAMaf_white <- read.maf(TCGAMaf_df_white, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA <- coBarplot(TCGAMaf_white, CBCGAMaf,
                         genes = rev(genes),
                         m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                         colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_TCGA <- CBCGAvsTCGA$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_CBCGA <- CBCGAvsTCGA$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_TCGA$Variant_Classification,
                                                  levels = rev(c(
                                                    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                    "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                  ))
)
CBCGAvsTCGA_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_CBCGA$Variant_Classification,
                                                   levels = rev(c(
                                                     "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                     "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                   ))
)
CBCGAvsTCGA_TCGA$name <- factor(CBCGAvsTCGA_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_CBCGA$name <- factor(CBCGAvsTCGA_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_TCGA <- CBCGAvsTCGA_TCGA[order(CBCGAvsTCGA_TCGA$name, CBCGAvsTCGA_TCGA$Variant_Classification), ]
CBCGAvsTCGA_CBCGA <- CBCGAvsTCGA_CBCGA[order(CBCGAvsTCGA_CBCGA$name, CBCGAvsTCGA_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_TCGA$value <- -CBCGAvsTCGA_TCGA$value

p2b <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
p2b
graph2pdf(p2b,file = '/results/Fig2B-2.pdf', width = 6, height = 4.3)

fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_white, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_white$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/All.csv', row.names = F)

# Figure 2C ---------------------------------------------------------------
# Prepare the data for Figure 2C, S2C-G
CBCGAClin_LumA <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "LumA")
CBCGAMaf_df_LumA <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_LumA$Tumor_Sample_Barcode)
CBCGAMaf_LumA <- read.maf(CBCGAMaf_df_LumA)
TCGAClin_white_LumA <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "LumA")
TCGAMaf_df_LumA <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_LumA$Tumor_Sample_Barcode)
TCGAMaf_LumA <- read.maf(TCGAMaf_df_LumA, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_LumA <- coBarplot(TCGAMaf_LumA, CBCGAMaf_LumA,
                              genes = rev(genes),
                              m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                              colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_LumA_TCGA <- CBCGAvsTCGA_LumA$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_LumA_CBCGA <- CBCGAvsTCGA_LumA$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_LumA_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_LumA_TCGA$Variant_Classification,
                                                       levels = rev(c(
                                                         "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                         "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                       ))
)
CBCGAvsTCGA_LumA_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_LumA_CBCGA$Variant_Classification,
                                                        levels = rev(c(
                                                          "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                          "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                        ))
)
CBCGAvsTCGA_LumA_TCGA$name <- factor(CBCGAvsTCGA_LumA_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_LumA_CBCGA$name <- factor(CBCGAvsTCGA_LumA_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_LumA_TCGA <- CBCGAvsTCGA_LumA_TCGA[order(CBCGAvsTCGA_LumA_TCGA$name, CBCGAvsTCGA_LumA_TCGA$Variant_Classification), ]
CBCGAvsTCGA_LumA_CBCGA <- CBCGAvsTCGA_LumA_CBCGA[order(CBCGAvsTCGA_LumA_CBCGA$name, CBCGAvsTCGA_LumA_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_LumA_TCGA$value <- -CBCGAvsTCGA_LumA_TCGA$value


CBCGAClin_LumB <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "LumB")
CBCGAMaf_df_LumB <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_LumB$Tumor_Sample_Barcode)
CBCGAMaf_LumB <- read.maf(CBCGAMaf_df_LumB)
TCGAClin_white_LumB <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "LumB")
TCGAMaf_df_LumB <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_LumB$Tumor_Sample_Barcode)
TCGAMaf_LumB <- read.maf(TCGAMaf_df_LumB, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_LumB <- coBarplot(TCGAMaf_LumB, CBCGAMaf_LumB,
                              genes = rev(genes),
                              m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                              colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_LumB_TCGA <- CBCGAvsTCGA_LumB$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_LumB_CBCGA <- CBCGAvsTCGA_LumB$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_LumB_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_LumB_TCGA$Variant_Classification,
                                                       levels = rev(c(
                                                         "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                         "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                       ))
)
CBCGAvsTCGA_LumB_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_LumB_CBCGA$Variant_Classification,
                                                        levels = rev(c(
                                                          "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                          "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                        ))
)
CBCGAvsTCGA_LumB_TCGA$name <- factor(CBCGAvsTCGA_LumB_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_LumB_CBCGA$name <- factor(CBCGAvsTCGA_LumB_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_LumB_TCGA <- CBCGAvsTCGA_LumB_TCGA[order(CBCGAvsTCGA_LumB_TCGA$name, CBCGAvsTCGA_LumB_TCGA$Variant_Classification), ]
CBCGAvsTCGA_LumB_CBCGA <- CBCGAvsTCGA_LumB_CBCGA[order(CBCGAvsTCGA_LumB_CBCGA$name, CBCGAvsTCGA_LumB_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_LumB_TCGA$value <- -CBCGAvsTCGA_LumB_TCGA$value


CBCGAClin_Her2 <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "Her2")
CBCGAMaf_df_Her2 <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_Her2$Tumor_Sample_Barcode)
CBCGAMaf_Her2 <- read.maf(CBCGAMaf_df_Her2)
TCGAClin_white_Her2 <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "Her2")
TCGAMaf_df_Her2 <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_Her2$Tumor_Sample_Barcode)
TCGAMaf_Her2 <- read.maf(TCGAMaf_df_Her2, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_Her2 <- coBarplot(TCGAMaf_Her2, CBCGAMaf_Her2,
                              genes = rev(genes),
                              m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                              colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_Her2_TCGA <- CBCGAvsTCGA_Her2$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_Her2_CBCGA <- CBCGAvsTCGA_Her2$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_Her2_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_Her2_TCGA$Variant_Classification,
                                                       levels = rev(c(
                                                         "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                         "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                       ))
)
CBCGAvsTCGA_Her2_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_Her2_CBCGA$Variant_Classification,
                                                        levels = rev(c(
                                                          "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                          "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                        ))
)
CBCGAvsTCGA_Her2_TCGA$name <- factor(CBCGAvsTCGA_Her2_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_Her2_CBCGA$name <- factor(CBCGAvsTCGA_Her2_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_Her2_TCGA <- CBCGAvsTCGA_Her2_TCGA[order(CBCGAvsTCGA_Her2_TCGA$name, CBCGAvsTCGA_Her2_TCGA$Variant_Classification), ]
CBCGAvsTCGA_Her2_CBCGA <- CBCGAvsTCGA_Her2_CBCGA[order(CBCGAvsTCGA_Her2_CBCGA$name, CBCGAvsTCGA_Her2_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_Her2_TCGA$value <- -CBCGAvsTCGA_Her2_TCGA$value


CBCGAClin_Basal <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "Basal")
CBCGAMaf_df_Basal <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_Basal$Tumor_Sample_Barcode)
CBCGAMaf_Basal <- read.maf(CBCGAMaf_df_Basal)
TCGAClin_white_Basal <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "Basal")
TCGAMaf_df_Basal <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_Basal$Tumor_Sample_Barcode)
TCGAMaf_Basal <- read.maf(TCGAMaf_df_Basal, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_Basal <- coBarplot(TCGAMaf_Basal, CBCGAMaf_Basal,
                               genes = rev(genes),
                               m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                               colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_Basal_TCGA <- CBCGAvsTCGA_Basal$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_Basal_CBCGA <- CBCGAvsTCGA_Basal$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_Basal_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_Basal_TCGA$Variant_Classification,
                                                        levels = rev(c(
                                                          "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                          "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                        ))
)
CBCGAvsTCGA_Basal_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_Basal_CBCGA$Variant_Classification,
                                                         levels = rev(c(
                                                           "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                           "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                         ))
)
CBCGAvsTCGA_Basal_TCGA$name <- factor(CBCGAvsTCGA_Basal_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_Basal_CBCGA$name <- factor(CBCGAvsTCGA_Basal_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_Basal_TCGA <- CBCGAvsTCGA_Basal_TCGA[order(CBCGAvsTCGA_Basal_TCGA$name, CBCGAvsTCGA_Basal_TCGA$Variant_Classification), ]
CBCGAvsTCGA_Basal_CBCGA <- CBCGAvsTCGA_Basal_CBCGA[order(CBCGAvsTCGA_Basal_CBCGA$name, CBCGAvsTCGA_Basal_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_Basal_TCGA$value <- -CBCGAvsTCGA_Basal_TCGA$value


CBCGAClin_Normal <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "Normal")
CBCGAMaf_df_Normal <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_Normal$Tumor_Sample_Barcode)
CBCGAMaf_Normal <- read.maf(CBCGAMaf_df_Normal)
TCGAClin_white_Normal <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "Normal")
TCGAMaf_df_Normal <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_Normal$Tumor_Sample_Barcode)
TCGAMaf_Normal <- read.maf(TCGAMaf_df_Normal, vc_nonSyn = vc_nonSyn)

CBCGAvsTCGA_Normal <- coBarplot(TCGAMaf_Normal, CBCGAMaf_Normal,
                                genes = rev(genes),
                                m1Name = "TCGA Caucasian", m2Name = "CBCGA",
                                colors = mycol, yLims = c(50, 50), showPct = F
)
CBCGAvsTCGA_Normal_TCGA <- CBCGAvsTCGA_Normal$m1 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "TCGA")
CBCGAvsTCGA_Normal_CBCGA <- CBCGAvsTCGA_Normal$m2 %>%
  data.frame() %>%
  mutate(Variant_Classification = rownames(.)) %>%
  pivot_longer(!Variant_Classification) %>%
  mutate(Cohort = "CBCGA")
CBCGAvsTCGA_Normal_TCGA$Variant_Classification <- factor(CBCGAvsTCGA_Normal_TCGA$Variant_Classification,
                                                         levels = rev(c(
                                                           "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                           "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                         ))
)
CBCGAvsTCGA_Normal_CBCGA$Variant_Classification <- factor(CBCGAvsTCGA_Normal_CBCGA$Variant_Classification,
                                                          levels = rev(c(
                                                            "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                                            "In_Frame_Ins", "In_Frame_Del", "Splice_Site", "Translation_Start_Site", "Multi_Hit"
                                                          ))
)
CBCGAvsTCGA_Normal_TCGA$name <- factor(CBCGAvsTCGA_Normal_TCGA$name, levels = rev(genes))
CBCGAvsTCGA_Normal_CBCGA$name <- factor(CBCGAvsTCGA_Normal_CBCGA$name, levels = rev(genes))
CBCGAvsTCGA_Normal_TCGA <- CBCGAvsTCGA_Normal_TCGA[order(CBCGAvsTCGA_Normal_TCGA$name, CBCGAvsTCGA_Normal_TCGA$Variant_Classification), ]
CBCGAvsTCGA_Normal_CBCGA <- CBCGAvsTCGA_Normal_CBCGA[order(CBCGAvsTCGA_Normal_CBCGA$name, CBCGAvsTCGA_Normal_CBCGA$Variant_Classification), ]
CBCGAvsTCGA_Normal_TCGA$value <- -CBCGAvsTCGA_Normal_TCGA$value


CBCGAvsTCGA_AKT1_LumA <- rbind(CBCGAvsTCGA_LumA_CBCGA, CBCGAvsTCGA_LumA_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "LumA")
CBCGAvsTCGA_AKT1_LumB <- rbind(CBCGAvsTCGA_LumB_CBCGA, CBCGAvsTCGA_LumB_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "LumB")
CBCGAvsTCGA_AKT1_Her2 <- rbind(CBCGAvsTCGA_Her2_CBCGA, CBCGAvsTCGA_Her2_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "Her2")
CBCGAvsTCGA_AKT1_Basal <- rbind(CBCGAvsTCGA_Basal_CBCGA, CBCGAvsTCGA_Basal_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "Basal")
CBCGAvsTCGA_AKT1_Normal <- rbind(CBCGAvsTCGA_Normal_CBCGA, CBCGAvsTCGA_Normal_TCGA) %>%
  filter(name == "AKT1") %>%
  mutate(Subtype = "Normal")
CBCGAvsTCGA_AKT1 <- rbind(
  CBCGAvsTCGA_AKT1_LumA, CBCGAvsTCGA_AKT1_LumB, CBCGAvsTCGA_AKT1_Her2,
  CBCGAvsTCGA_AKT1_Basal, CBCGAvsTCGA_AKT1_Normal
)
CBCGAvsTCGA_AKT1$Subtype <- factor(CBCGAvsTCGA_AKT1$Subtype, levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))


p2c <- ggplot() +
  geom_bar(aes(Subtype, value, fill = Variant_Classification),
           data = subset(CBCGAvsTCGA_AKT1, Cohort == "CBCGA"), position = "stack", stat = "identity"
  ) +
  geom_bar(aes(Subtype, value, fill = Variant_Classification),
           data = subset(CBCGAvsTCGA_AKT1, Cohort == "TCGA"), position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "", angle = 45),
    y = guide_axis(title = "Prevalence"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = abs) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
p2c
graph2pdf(p2c,file = '/results/Fig2C.pdf', width = 4.5, height = 3.5)





# Figure 2D ---------------------------------------------------------------
MET_Clin_LumA <- subset(MET_Clin, Pam50_SUBTYPE == "LumA")
MET_Clin_nonLumA <- subset(MET_Clin, Pam50_SUBTYPE %in% c("LumB", "Her2", "Basal", "Normal"))
METMaf_df <- subset(METMaf_df, Tumor_Sample_Barcode %in% MET_Clin$PATIENT_ID)
METMaf_df_LumA <- subset(METMaf_df, Tumor_Sample_Barcode %in% MET_Clin_LumA$PATIENT_ID)
METMaf_df_nonLumA <- subset(METMaf_df, Tumor_Sample_Barcode %in% MET_Clin_nonLumA$PATIENT_ID)

MET_all <- length(unique(METMaf_df$Tumor_Sample_Barcode))
MET_all_LumA <- length(unique(METMaf_df_LumA$Tumor_Sample_Barcode))
MET_all_nonLUmA <- length(unique(METMaf_df_nonLumA$Tumor_Sample_Barcode))

METMaf <- read.maf(METMaf_df, vc_nonSyn = vc_nonSyn)
METMaf_LumA <- read.maf(METMaf_df_LumA, vc_nonSyn = vc_nonSyn)
METMaf_nonLumA <- read.maf(METMaf_df_nonLumA, vc_nonSyn = vc_nonSyn)

MET_AKT1 <- genesToBarcodes(METMaf, genes = "AKT1", justNames = T)[[1]] %>% length() # 78/1866
MET_AKT1_LumA <- genesToBarcodes(METMaf_LumA, genes = "AKT1", justNames = T)[[1]] %>% length() # 38/687
MET_AKT1_nonLumA <- genesToBarcodes(METMaf_nonLumA, genes = "AKT1", justNames = T)[[1]] %>% length() # 38/1173


AKT1_df <- data.frame(
  PAM50 = c(rep(c("All", "Luminal A", "non-Luminal A"), 3), "All", "All"),
  Cohort = c(rep("CBCGA", 3), rep("TCGA Caucasian", 3), rep("METABRIC", 3), "GDPH", "NCCH"),
  Freq = c(40 / 624, 22 / 182, 18 / 425, 12 / 474, 10 / 229, 2 / 233, MET_AKT1 / MET_all, MET_AKT1_LumA / MET_all_LumA, MET_AKT1_nonLumA / MET_all_nonLUmA, 0.0706, 0.0740)
)
AKT1_df$Freq2 <- paste0(round(AKT1_df$Freq, 3) * 100, "%")
AKT1_df$Cohort <- factor(AKT1_df$Cohort, levels = c("CBCGA", "GDPH", "NCCH", "TCGA Caucasian", "METABRIC"))

p  = ggplot(AKT1_df, aes(PAM50, Freq, fill = PAM50)) +
  geom_bar(position = "stack", stat = "identity", colour = "black") +
  geom_text(aes(y = Freq + 0.005, label = Freq2)) +
  scale_fill_manual(values = c("#F2AF00", "#009CA6", "#CECECE")) +
  scale_y_continuous("Alteration frequency (%)",
                     labels = scales::percent, limits = c(0, 0.15),
                     expand = expansion(mult = c(0, 0.03), add = 0)
  ) +
  guides(x = guide_axis(title = ""), fill = guide_legend(title = "")) +
  facet_grid(. ~ Cohort, scales = "free_x", space = "free_x", switch = "x") +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.text.x = element_blank(),
    legend.position = c(0.9, 0.9),
    legend.background = element_blank(),
    strip.background = element_blank()
  )
graph2pdf(p,file = '/results/Fig2D.pdf', width = 4.5, height = 3.5)


# Figure 2E ---------------------------------------------------------------
CBCGAClin_LumA <- subset(CBCGAClin_WES, intrinsic_subtype_pam50 == "LumA")
CBCGAMaf_df_LumA <- subset(CBCGA_WES_Somatic, Tumor_Sample_Barcode %in% CBCGAClin_LumA$Tumor_Sample_Barcode)
CBCGAMaf_LumA <- read.maf(CBCGAMaf_df_LumA, vc_nonSyn = vc_nonSyn)
TCGAClin_white_LumA <- subset(TCGAClin_WES_white, PAM50Call_RNAseq == "LumA")
TCGAMaf_df_LumA <- subset(TCGAMaf_df, Tumor_Sample_Barcode %in% TCGAClin_white_LumA$Tumor_Sample_Barcode)
TCGAMaf_LumA <- read.maf(TCGAMaf_df_LumA)

CBCGAClin_LumA_pre <- subset(CBCGAClin_LumA, menopause == "No")
CBCGAMaf_df_LumA_pre <- subset(CBCGAMaf_df_LumA, Tumor_Sample_Barcode %in% CBCGAClin_LumA_pre$Tumor_Sample_Barcode)
CBCGAMaf_LumA_pre <- read.maf(CBCGAMaf_df_LumA_pre)
TCGAClin_white_LumA_Pre <- subset(TCGAClin_white_LumA, menopause_status == "Pre")
TCGAMaf_df_LumA_pre <- subset(TCGAMaf_df_LumA, Tumor_Sample_Barcode %in% TCGAClin_white_LumA_Pre$Tumor_Sample_Barcode)
TCGAMaf_LumA_pre <- read.maf(TCGAMaf_df_LumA_pre)

CBCGA_pre_all <- unique(CBCGAMaf_df_LumA_pre$Tumor_Sample_Barcode) %>% length()
CBCGA_pre_mut <- genesToBarcodes(CBCGAMaf_LumA_pre, "AKT1", justNames = T)[[1]] %>% length()
CBCGA_pre_nonmut <- CBCGA_pre_all - CBCGA_pre_mut

TCGA_pre_all <- unique(TCGAMaf_df_LumA_pre$Tumor_Sample_Barcode) %>% length()
TCGA_pre_mut <- genesToBarcodes(TCGAMaf_LumA_pre, "AKT1", justNames = T)[[1]] %>% length()
TCGA_pre_nonmut <- TCGA_pre_all - TCGA_pre_mut


CBCGAClin_LumA_post <- subset(CBCGAClin_LumA, menopause == "Yes")
CBCGAMaf_df_LumA_post <- subset(CBCGAMaf_df_LumA, Tumor_Sample_Barcode %in% CBCGAClin_LumA_post$Tumor_Sample_Barcode)
CBCGAMaf_LumA_post <- read.maf(CBCGAMaf_df_LumA_post)
TCGAClin_white_LumA_post <- subset(TCGAClin_white_LumA, menopause_status == "Post")
TCGAMaf_df_LumA_post <- subset(TCGAMaf_df_LumA, Tumor_Sample_Barcode %in% TCGAClin_white_LumA_post$Tumor_Sample_Barcode)
TCGAMaf_LumA_post <- read.maf(TCGAMaf_df_LumA_post)

CBCGA_post_all <- unique(CBCGAMaf_df_LumA_post$Tumor_Sample_Barcode) %>% length()
CBCGA_post_mut <- genesToBarcodes(CBCGAMaf_LumA_post, "AKT1", justNames = T)[[1]] %>% length()
CBCGA_post_nonmut <- CBCGA_post_all - CBCGA_post_mut

TCGA_post_all <- unique(TCGAMaf_df_LumA_post$Tumor_Sample_Barcode) %>% length()
TCGA_post_mut <- genesToBarcodes(TCGAMaf_LumA_post, "AKT1", justNames = T)[[1]] %>% length()
TCGA_post_nonmut <- TCGA_post_all - TCGA_post_mut



AKT1_meno_df <- data.frame(
  Cohort = c(rep("CBCGA", 4), rep("TCGA Caucasian", 4)),
  Meno = rep(c("Pre-menopausal", "Pre-menopausal", "Post-menopausal", "Post-menopausal"), 2),
  Mutation = rep(c("Mut", "Wt"), 4),
  Count = c(
    CBCGA_pre_mut, CBCGA_pre_nonmut, CBCGA_post_mut, CBCGA_post_nonmut,
    TCGA_pre_mut, TCGA_pre_nonmut, TCGA_post_mut, TCGA_post_nonmut
  ),
  Percent = c(
    CBCGA_pre_mut / CBCGA_pre_all, NA, CBCGA_post_mut / CBCGA_post_all, NA,
    TCGA_pre_mut / TCGA_pre_all, NA, TCGA_post_mut / TCGA_post_all, NA
  )
)
AKT1_meno_df$Percent <- round(AKT1_meno_df$Percent, 3)
AKT1_meno_df$Percent2 <- ifelse(is.na(AKT1_meno_df$Percent), NA,
                                paste0(AKT1_meno_df$Percent * 100, "%")
)

AKT1_meno_df$Meno <- factor(AKT1_meno_df$Meno, levels = c("Pre-menopausal", "Post-menopausal"))
AKT1_meno_df$Mutation2 <- ifelse(AKT1_meno_df$Mutation == "Wt", "Wt",
                                 paste(AKT1_meno_df$Cohort, "Mut")
)
AKT1_meno_df$Mutation2 <- factor(AKT1_meno_df$Mutation2, levels = c("Wt", "CBCGA Mut", "TCGA Caucasian Mut"))


p = ggplot(AKT1_meno_df, aes(Cohort, Count)) +
  geom_bar(position = "fill", stat = "identity", aes(fill = Mutation2)) +
  geom_text(aes(y = Percent + 0.03, label = Percent2)) +
  # geom_signif(
  #   data = data.frame(Meno = factor(c("Pre-menopausal", "Post-menopausal"),
  #                                   levels = c("Pre-menopausal", "Post-menopausal")
  #   )),
  #   aes(
  #     y_position = c(1.01, 1.01), xmin = c(1, 1), xmax = c(2, 2),
  #     annotations = c("P = 0.031", "P = 0.144")
  #   ), tip_length = 0.00005, manual = T
  # ) +
  scale_fill_manual(values = c("#F2F2F2", "#0077C9", "#91B720")) +
  scale_y_continuous("AKT1 mutation frequency (%)",
                     labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))
  ) +
  guides(
    x = guide_axis(title = ""),
    fill = guide_legend(title = "")
  ) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    strip.background = element_blank()
  ) +
  facet_wrap(~Meno, ncol = 2, strip.position = "bottom")
graph2pdf(p,file = "/results/Fig2E.pdf", width = 6, height = 4.3)
fisher.test(matrix(c(CBCGA_pre_mut, CBCGA_pre_nonmut, TCGA_pre_mut, TCGA_pre_nonmut), byrow = F, ncol = 2))
fisher.test(matrix(c(CBCGA_post_mut, CBCGA_post_nonmut, TCGA_post_mut, TCGA_post_nonmut), byrow = F, ncol = 2))


# Figure 2F ---------------------------------------------------------------
CBCGAMaf_df_LumA2 <- CBCGAMaf_LumA@data
CBCGAMaf_df_LumA2 <- CBCGAMaf_df_LumA2[!is.na(CBCGAMaf_df_LumA2$aaChange), ]
TCGAMaf_df_LumA2 <- TCGAMaf_LumA@data
TCGAMaf_df_LumA2 <- TCGAMaf_df_LumA2[!is.na(TCGAMaf_df_LumA2$HGVSp_Short), ]

extractpos <- function(maf_aachange) {
  prot.spl <- strsplit(x = as.character(maf_aachange), split = ".", fixed = TRUE)
  prot.conv <- sapply(sapply(prot.spl, function(x) x[length(x)]), "[", 1)
  pos <- gsub(pattern = "Ter.*", replacement = "", x = prot.conv)
  pos <- gsub(pattern = "[[:alpha:]]", replacement = "", x = pos)
  pos <- gsub(pattern = "\\*$", replacement = "", x = pos)
  pos <- gsub(pattern = "^\\*", replacement = "", x = pos)
  pos <- gsub(pattern = "\\*.*", replacement = "", x = pos)
  pos <- as.numeric(sapply(X = strsplit(x = pos, split = "_", fixed = TRUE), FUN = function(x) x[1]))
  aa <- paste0(unlist(regmatches(maf_aachange, gregexpr("p[.].[0-9]+", maf_aachange))), "X")
  mutpos <- data.frame(position = pos, mutation = maf_aachange, aa = aa, stringsAsFactors = F)
  return(mutpos[order(mutpos$position), ])
}
posCBCGA <- extractpos(CBCGAMaf_df_LumA2[CBCGAMaf_df_LumA2$Hugo_Symbol == "AKT1", ]$aaChange)
posTCGA <- extractpos(TCGAMaf_df_LumA2[TCGAMaf_df_LumA2$Hugo_Symbol == "AKT1", ]$HGVSp_Short)

nrposCBCGA <- posCBCGA[!duplicated(posCBCGA), ]
rownames(nrposCBCGA) <- nrposCBCGA$mutation
nrposCBCGA$rate <- 1
nrposCBCGA[names(table(nrposCBCGA$mutation)), ]$rate <- table(posCBCGA$mutation)
nrposCBCGA$side <- "top"
head(nrposCBCGA)
nrposTCGA <- posTCGA[!duplicated(posTCGA), ]
rownames(nrposTCGA) <- nrposTCGA$mutation
nrposTCGA$rate <- 1
nrposTCGA[names(table(nrposTCGA$mutation)), ]$rate <- table(posTCGA$mutation)
head(nrposTCGA)
nrposTCGA$side <- "bottom"
nrpos <- rbind(nrposCBCGA, nrposTCGA)

features <- GRanges(
  "chr14",
  IRanges(
    start = c(1, 5, 150, 409),
    end = c(4, 108, 408, 480),
    names = c("", "PH", "Protein kinase", "AGC_kinase")
  )
)
features$height <- c(0, 0.05, 0.07, 0.07)
features$fill <- c(NA, "#FF8833", "#51C6E6", "#DFA32D")

sample.gr <- GRanges("chr14", IRanges(nrpos$position, width = 1, names = nrpos$mutation))
sample.gr$label.parameter.rot <- 90
sample.gr$label <- as.character(nrpos$rate)
sample.gr$label.col <- "white"
sample.gr$SNPsideID <- nrpos$side
sample.gr$color <- c(rep("#DC5035", nrow(nrposCBCGA)), rep("#0085C4", nrow(nrposTCGA)))
sample.gr$border <- "black"
sample.gr$alpha <- rep(1, nrow(nrpos))
sample.gr$score <- log2(nrpos$rate)

pdf("/results/Fig2F.pdf")
par(mar = c(0, 1, 0, 1))
lolliplot(sample.gr, features,
          xaxis = c(1, 100, 200, 300, 400, 480),
          yaxis = F, ylab = F, type = "circle", label_on_feature = T
)
dev.off()

# Figure 2G ---------------------------------------------------------------
# CBCGAMaf_df_LumA3 <- CBCGAMaf_LumA@data
# CBCGAMaf_df_LumA3$Chromosome <- gsub("chr", "", CBCGAMaf_df_LumA3$Chromosome)
# CBCGAMaf_df_LumA3$mutation_id <- paste(CBCGAMaf_df_LumA3$Tumor_Sample_Barcode, CBCGAMaf_df_LumA3$Chromosome,
#                                        CBCGAMaf_df_LumA3$Start_Position, CBCGAMaf_df_LumA3$Reference_Allele,
#                                        sep = ":"
# )
# 
# DriverGene <- c(
#   "TP53", "AKT1", "PIK3CA", "MAP3K1", "PTEN", "MAP2K4", "CBFB", "GATA3", "KMT2C", "CDH1", "RUNX1", "FOXA1", "CDKN1B", "GPS2", "ARID1A",
#   "ZFP36L1", "NF1", "PIK3R1", "SF3B1", "CTCF", "TBX3", "SMAD4"
# )
# CBCGA_pyclone_mutation_LumA <- subset(CBCGA_pyclone_mutation, Sample %in% CBCGAClin_LumA$Tumor_Sample_Barcode)
# CBCGA_pyclone_mutation_LumA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene %in% DriverGene)
# CBCGA_pyclone_mutation_LumA <- CBCGA_pyclone_mutation_LumA[intersect(rownames(CBCGA_pyclone_mutation_LumA), unique(CBCGAMaf_df_LumA3$mutation_id)), ]
# 
# 
# genePair <- combn(DriverGene, 2) %>%
#   t() %>%
#   data.frame()
# genePair$CBCGA <- NA
# if (T) {
#   for (i in 1:nrow(genePair)) {
#     tmp_ID1_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair[i, 1])$Sample %>% unique()
#     tmp_ID2_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair[i, 2])$Sample %>% unique()
#     tmp_ID_CBCGA <- intersect(tmp_ID1_CBCGA, tmp_ID2_CBCGA) %>% length()
#     genePair[i, "CBCGA"] <- tmp_ID_CBCGA
#   }
#   
#   
#   genePair_CBCGA <- subset(genePair, CBCGA >= 2)
#   for (i in 1:nrow(genePair_CBCGA)) {
#     tmp_ID1_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair_CBCGA[i, 1])$Sample %>% unique()
#     tmp_ID2_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair_CBCGA[i, 2])$Sample %>% unique()
#     tmp_ID_CBCGA <- intersect(tmp_ID1_CBCGA, tmp_ID2_CBCGA)
#     
#     tmp_win <- NULL
#     for (j in tmp_ID_CBCGA) {
#       tmp_VAF1_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair_CBCGA[i, 1] & Sample == j)$PostClusterCCF %>% mean()
#       tmp_VAF2_CBCGA <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene == genePair_CBCGA[i, 2] & Sample == j)$PostClusterCCF %>% mean()
#       tmp_win2 <- ifelse(tmp_VAF1_CBCGA - tmp_VAF2_CBCGA >= 0.05, "CBCGA_win1",
#                          ifelse(tmp_VAF1_CBCGA - tmp_VAF2_CBCGA <= -0.05, "CBCGA_win2", "CBCGA_amb")
#       )
#       tmp_win <- c(tmp_win, tmp_win2)
#     }
#     tmp_df_CBCGA <- table(tmp_win) %>% data.frame()
#     if ("CBCGA_win1" %in% tmp_df_CBCGA$tmp_win) {
#       genePair_CBCGA[i, "CBCGA_win1"] <- tmp_df_CBCGA[tmp_df_CBCGA$tmp_win == "CBCGA_win1", "Freq"]
#     } else {
#       genePair_CBCGA[i, "CBCGA_win1"] <- 0
#     }
#     if ("CBCGA_win2" %in% tmp_df_CBCGA$tmp_win) {
#       genePair_CBCGA[i, "CBCGA_win2"] <- tmp_df_CBCGA[tmp_df_CBCGA$tmp_win == "CBCGA_win2", "Freq"]
#     } else {
#       genePair_CBCGA[i, "CBCGA_win2"] <- 0
#     }
#     if ("CBCGA_amb" %in% tmp_df_CBCGA$tmp_win) {
#       genePair_CBCGA[i, "CBCGA_amb"] <- tmp_df_CBCGA[tmp_df_CBCGA$tmp_win == "CBCGA_amb", "Freq"]
#     } else {
#       genePair_CBCGA[i, "CBCGA_amb"] <- 0
#     }
#   }
#   
#   
#   # Bradley-Terry2
#   colnames(genePair_CBCGA)[c(1, 2, 4, 5)] <- c("Gene1", "Gene2", "Gene1_wins", "Gene2_wins")
#   genePair_CBCGA <- genePair_CBCGA[, c(1, 2, 4, 5)]
#   genePair_CBCGA$Gene1 <- factor(genePair_CBCGA$Gene1, levels = unique(c(genePair_CBCGA$Gene1, genePair_CBCGA$Gene2)))
#   genePair_CBCGA$Gene2 <- factor(genePair_CBCGA$Gene2, levels = unique(c(as.character(genePair_CBCGA$Gene1), genePair_CBCGA$Gene2)))
#   CBCGA_btData <- BTm(cbind(Gene1_wins, Gene2_wins), Gene1, Gene2, data = genePair_CBCGA)
#   CBCGA_btData_ordering <- as.data.frame(summary(CBCGA_btData)$coefficients)
#   CBCGA_btData_ordering$Symbol <- rownames(CBCGA_btData_ordering)
#   CBCGA_btData_ordering$Symbol <- substring(CBCGA_btData_ordering$Symbol, 3)
#   CBCGA_btData_ordering <- subset(CBCGA_btData_ordering, `Std. Error` <= 2)
#   CBCGA_btData_ordering$Pathway <- c(
#     "PI3K", "PI3K", "Other", "Transcription factor", "Transcription factor",
#     "Chromatin histone modifiers", "Splicing", "Transcription factor", "Cell cycle"
#   )
#   
#   mycol2 <- c("#A94543", "#1F97C6", "#75A992", "#D68A46", "#334D54", "#666393")
#   names(mycol2) <- c(
#     "PI3K", "Cell cycle", "Transcription factor", "Chromatin histone modifiers",
#     "Splicing", "Other"
#   )
#   p2g1 <- ggplot(
#     CBCGA_btData_ordering,
#     aes(reorder(factor(Symbol), Estimate),
#         y = Estimate,
#         ymin = (Estimate - `Std. Error`), ymax = (Estimate + `Std. Error`)
#     )
#   ) +
#     geom_pointrange(size = 0.75, aes(color = Pathway)) +
#     scale_color_manual(values = mycol2) +
#     scale_y_reverse() +
#     guides(
#       x = guide_axis(title = "Point Estimate + 95% CI"),
#       y = guide_axis(title = ""),
#       color = guide_legend(title = "")
#     ) +
#     coord_flip() +
#     theme_classic() +
#     theme(
#       axis.line.y = element_blank(),
#       axis.ticks.y = element_blank(),
#       axis.text = element_text(colour = "black", size = 10),
#       axis.title = element_text(size = 12)
#     )
#   graph2pdf(p2g1,file = "/results/Fig2G_1.pdf")
# }
# 
# CBCGA_btData_ordering <- CBCGA_btData_ordering[order(CBCGA_btData_ordering$Estimate, decreasing = T), ]
# CBCGA_pyclone_mutation_LumA2 <- subset(CBCGA_pyclone_mutation_LumA, Gene.refGene %in% CBCGA_btData_ordering$Symbol)
# CBCGA_pyclone_mutation_LumA2$Gene.refGene <- factor(CBCGA_pyclone_mutation_LumA2$Gene.refGene, levels = rev(CBCGA_btData_ordering$Symbol))
# CBCGA_pyclone_mutation_LumA2 <- merge(CBCGA_pyclone_mutation_LumA2, CBCGA_btData_ordering[, 5:6], by.x = "Gene.refGene", by.y = "Symbol", all.x = T)
# p2g2 <- ggplot(CBCGA_pyclone_mutation_LumA2, aes(PostClusterCCF, Gene.refGene, fill = Pathway)) +
#   geom_density_ridges(scale = 1, rel_min_height = 0.01) +
#   scale_fill_manual(values = mycol2) +
#   guides(
#     y = guide_axis(title = ""),
#     fill = guide_legend(title = "")
#   ) +
#   theme_classic()
# ggarrange(p2g2, p2g1, ncol = 2, common.legend = T, legend = "right")
# graph2pdf(p2g2,file = "/results/Fig2G_2.pdf")

# Figure S2A --------------------------------------------------------------
CBCGAClin_Age <- data.frame(Age = CBCGAClin_WES$age_at_surgery, Cohort = "CBCGA", PAM50 = CBCGAClin_WES$intrinsic_subtype_pam50)
CBCGAClin_Age$PAM50[is.na(CBCGAClin_Age$PAM50)] <- ""
TCGAClin_WES_Asian <- subset(TCGAClin_WES, Race.Category == "ASIAN" & histological_type == "Infiltrating Ductal Carcinoma")
TCGAClin_Age <- data.frame(
  Age = c(TCGAClin_WES_white$age_at_initial_pathologic_diagnosis, TCGAClin_WES_Asian$age_at_initial_pathologic_diagnosis),
  Cohort = c(rep("TCGA Caucasian", nrow(TCGAClin_WES_white)), rep("TCGA Asian", nrow(TCGAClin_WES_Asian))),
  PAM50 = c(TCGAClin_WES_white$PAM50Call_RNAseq, TCGAClin_WES_Asian$PAM50Call_RNAseq)
)
MET_Clin_Age <- data.frame(Age = MET_Clin$AGE_AT_DIAGNOSIS, Cohort = "METABRIC", PAM50 = MET_Clin$Pam50_SUBTYPE)
Age_compr <- rbind(CBCGAClin_Age, TCGAClin_Age, MET_Clin_Age)
Age_compr <- subset(Age_compr, PAM50 %in% c("LumA", "LumB", "Her2", "Basal", "Normal"))
Age_compr$PAM50 <- factor(Age_compr$PAM50, levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))
Age_compr$Cohort <- factor(Age_compr$Cohort, levels = c("CBCGA", "TCGA Asian", "TCGA Caucasian", "METABRIC"))
# escape::ridgeEnrichment(Age_compr,
#   gene.set = "Age", group = "Cohort",
#   add.rug = TRUE, facet = "PAM50",
#   colors = c("#0087BE", "#00B140", "#FFC600")
# ) +
#   guides(
#     x = guide_axis(title = "Age"),
#     y = guide_axis(title = "")
#   ) +
#   hrbrthemes::theme_ipsum(
#     base_family = "Arial", axis_title_size = 10, axis_text_size = 10,
#     strip_text_size = 10, axis_title_just = "m"
#   ) +
#   theme(axis.text = element_text(colour = "black"))
# graph2pdf(file = 'Figure/Figure S2A-2', width = 13.5, height = 4.5)
aggregate(Age ~ Cohort + PAM50, Age_compr, median)



# Figure S2B --------------------------------------------------------------
CBCGAvsTCGA_age2 <- data.frame(
  age = c(CBCGAClin_WES$age_at_surgery, CBCGAClin_WES$age_at_surgery, TCGAClin_WES_white$age_at_initial_pathologic_diagnosis, TCGAClin_WES_white$age_at_initial_pathologic_diagnosis),
  cohort = c(rep("CBCGA", nrow(CBCGAClin_WES)), rep("CBCGA", nrow(CBCGAClin_WES)), rep("TCGA", nrow(TCGAClin_WES_white)), rep("TCGA", nrow(TCGAClin_WES_white))),
  ihc = c(rep("All", nrow(CBCGAClin_WES)), CBCGAClin_WES$clinical_subtype, rep("All", nrow(TCGAClin_WES_white)), TCGAClin_WES_white$Clinical_Subtype)
)
CBCGAvsTCGA_age2 <- subset(CBCGAvsTCGA_age2, ihc != "Unknown")
CBCGAvsTCGA_age2$ihc <- paste(CBCGAvsTCGA_age2$ihc, CBCGAvsTCGA_age2$cohort, sep = "_")
CBCGAvsTCGA_age2$ihc <- factor(CBCGAvsTCGA_age2$ihc, levels = rev(c(
  "All_CBCGA", "All_TCGA", "HR+HER2-_CBCGA", "HR+HER2-_TCGA",
  "HR+HER2+_CBCGA", "HR+HER2+_TCGA", "HR-HER2+_CBCGA", "HR-HER2+_TCGA",
  "TNBC_CBCGA", "TNBC_TCGA"
)))

p = ggplot(CBCGAvsTCGA_age2, aes(x = as.numeric(age), y = ihc, fill = cohort)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01) +
  scale_fill_manual(values = c("CBCGA" = "#0077C9", "TCGA" = "#91B720")) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(
    x = guide_axis(title = "Age"),
    fill = guide_legend(title = "")
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(colour = "grey"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "bottom"
  )
graph2pdf(p,file = '/results/FigS2B-2.pdf', width = 6, height = 6.9)

aggregate(age ~ ihc, CBCGAvsTCGA_age2, median)
pval <- NULL
for (i in c("All", "HR+HER2-", "HR+HER2+", "HR-HER2+", "TNBC")) {
  tmp_function <- paste("age", "ihc", sep = "~") %>% as.formula()
  tmp_df <- subset(CBCGAvsTCGA_age2, ihc %in% paste(i, c("CBCGA", "TCGA"), sep = "_"))
  tmp_wcox <- wilcox.test(tmp_function, tmp_df)
  tmp_pval <- tmp_wcox$p.value
  pval <- c(pval, tmp_pval)
}
fdr <- p.adjust(as.numeric(pval), method = "fdr")
fdr

# Figure S2C --------------------------------------------------------------
pS2C <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_LumA_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_LumA_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
graph2pdf(pS2C,file = '/results/FigS2C.pdf', width = 6, height = 6.9)

fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_LumA, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_LumA$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_LumA, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_LumA$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/Luminal A.csv', row.names = F)




# Figure S2D ---------------------------------------------------------------
pS2D <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_LumB_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_LumB_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
graph2pdf(pS2D,file = '/results/FigS2D.pdf', width = 6, height = 6.9)


fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_LumB, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_LumB$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_LumB, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_LumB$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/Luminal B.csv', row.names = F)





# Figure S2E ---------------------------------------------------------------
pS2E <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Her2_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Her2_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(limits = c(-80, 80), breaks = seq(-80, 80, 40), labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
graph2pdf(pS2E,file = '/results/FigS2E.pdf', width = 6, height = 6.9)


fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_Her2, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_Her2$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_Her2, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_Her2$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/HER2 enriched.csv', row.names = F)





# Figure S2F ---------------------------------------------------------------
pS2F <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Basal_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Basal_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(limits = c(-80, 80), breaks = seq(-80, 80, 40), labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
graph2pdf(pS2F,file = '/results/FigS2F.pdf', width = 6, height = 6.9)


fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_Basal, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_Basal$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_Basal, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_Basal$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/Basal.csv', row.names = F)





# Figure S2G ---------------------------------------------------------------
pS2G <- ggplot() +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Normal_CBCGA, position = "stack", stat = "identity"
  ) +
  geom_bar(aes(name, value, fill = Variant_Classification),
           data = CBCGAvsTCGA_Normal_TCGA, position = "stack", stat = "identity"
  ) +
  geom_hline(yintercept = 0) +
  guides(
    x = guide_axis(title = "Prevalence"),
    y = guide_axis(title = "Hugo symbol"),
    fill = guide_legend(title = "")
  ) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, 30), labels = abs) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 10, colour = "black")
  )
pS2G
# ggarrange(pS2C, pS2D, pS2E, pS2F, pS2G,
#           ncol = 2, nrow = 3, labels = c("C", "D", "E", "F", "G"),
#           legend = "none"
# )
graph2pdf(pS2G,file = '/results/FigS2G.pdf', width = 6, height = 6.9)


fisher_table <- data.frame()
for (i in genes) {
  tmp_CBCGA_mut <- genesToBarcodes(CBCGAMaf_Normal, genes = i, justNames = T)[[1]] %>% length()
  tmp_CBCGA_wt <- length(unique(CBCGAMaf_df_Normal$Tumor_Sample_Barcode)) - tmp_CBCGA_mut
  tmp_CBCGA_freq <- tmp_CBCGA_mut / (tmp_CBCGA_mut + tmp_CBCGA_wt)
  tmp_TCGA_mut <- genesToBarcodes(TCGAMaf_Normal, genes = i, justNames = T)[[1]] %>% length()
  tmp_TCGA_wt <- length(unique(TCGAMaf_df_Normal$Tumor_Sample_Barcode)) - tmp_TCGA_mut
  tmp_TCGA_freq <- tmp_TCGA_mut / (tmp_TCGA_mut + tmp_TCGA_wt)
  
  tmp_fisher <- fisher.test(matrix(c(tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_TCGA_mut, tmp_TCGA_wt), byrow = F, ncol = 2))
  tmp_pval <- tmp_fisher$p.value
  tmp_df <- c(i, tmp_CBCGA_mut, tmp_CBCGA_wt, tmp_CBCGA_freq, tmp_TCGA_mut, tmp_TCGA_wt, tmp_TCGA_freq, tmp_pval)
  fisher_table <- rbind(fisher_table, tmp_df)
}
colnames(fisher_table) <- c("Symbol", "CBCGA_mut", "CBCGA_wt", "CBCGA_freq", "TCGA_mut", "TCGA_wt", "TCGA_freq", "pval")
fisher_table[, 2:8] <- fisher_table[, 2:8] %>% mutate_if(is.character, as.numeric)
fisher_table$FDR <- p.adjust(fisher_table$pval, method = "fdr")
# write.csv(fisher_table, 'Figure/Normal.csv', row.names = F)
