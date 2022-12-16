#-----------------------------------------------------------------#
# Fig6_S6 
#-----------------------------------------------------------------#
rm(list = ls())
options(stringsAsFactors = F)
gc()
# library(readxl)
# library(stringr)
# library(dplyr)
library(tidyverse)
library(GSVA)
library(GSEABase)
#library(ggpubr)
library(pheatmap)
####Load_expr_data
load("/data/CBCGA.Extended_MergedData_V2.5_220722.Rdata")
load("/data/Fig6_S6_data.Rdata")

sample_TT <- intersect(colnames(CBCGA.Extended_RNA_913_FPKM_symbol),paste(CBCGA_Cohort.Info$PatientCode,"RNA_T",sep = "_"))
expr_data_TT  <- CBCGA.Extended_RNA_913_FPKM_symbol[,sample_TT]

FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
expr_data_TT_tpm <- as.data.frame(apply(expr_data_TT,2,FPKM2TPM))

expr_data_TT_LOGtpm <- log2(expr_data_TT_tpm+1)


####------------------Fig 6A??????????????????????????????????--------------
##devolution??ssGSEA
dev_signature_CIBERSORT_trans <-list()
for (i in unique(dev_signature_CIBERSORT$CellType)){
  dev_signature_CIBERSORT_trans[[i]] <-  c(dev_signature_CIBERSORT[dev_signature_CIBERSORT$CellType == i,1])
}

table(dev_signature_CIBERSORT$Symbol%in%rownames(expr_data_TT))

res_dev_all <- gsva(as.matrix(expr_data_TT_LOGtpm),dev_signature_CIBERSORT_trans,method="ssgsea",parallel.sz=1, min.sz=1, max.sz=Inf, verbose=T,kcdf='Gaussian',abs.ranking=TRUE)
##cluster
set.seed(123)
res_dev_scale <- res_dev_all
for (i in 1:ncol(res_dev_scale)){
  res_dev_scale[,i] <-scale(res_dev_all[,i])
}
km1 <- kmeans(t(res_dev_scale),3,iter.max = 1000,nstart=25)
table(km1$cluster)
Patient_order <- names(sort(km1$cluster))
res_dev_scale_reorder <- res_dev_all[,Patient_order]

anno_col <- data.frame(
  row.names = Patient_order,
  Immune_cluster = sort(km1$cluster)
)

##heatmap
TILs_data <- CBCGA_TILs_F
rownames(TILs_data) <- paste(TILs_data$PatientCode,"RNA_T",sep = "_")
TILs_data <- TILs_data[rownames(anno_col),]
anno_col$TILs <- as.numeric(TILs_data$sTILs)

pic_data <-res_dev_all[,rownames(anno_col)]
bk = unique(c(seq(-2,2, length=100)))
mycolors <- c(colorRampPalette(c("#4E85AC", "white"))(50), colorRampPalette(c("white", "#C0052A"))(50))

pheatmap(as.matrix(pic_data), color = mycolors, breaks=bk ,fontsize_row = 8, cluster_rows = F,cluster_cols=F,show_colnames = F, 
         clustering_distance_cols = "euclidean", treeheight_col = 50, scale = "row",
         annotation_col = anno_col,
         filename = "/results/Fig6A_1.pdf")

anno_col$Immune_cluster[anno_col$Immune_cluster=="1"] <- "Cluster1"
anno_col$Immune_cluster[anno_col$Immune_cluster=="2"] <- "Cluster2"
anno_col$Immune_cluster[anno_col$Immune_cluster=="3"] <- "Cluster3"
anno_col <- arrange(anno_col,anno_col$Immune_cluster)
#  annotation color
anno_color <-list(
  TILs =  c(colorRampPalette(c("#F7D8D3", "#DC4F35"))(15)),
  Immune_cluster=c("Cluster1"="#637EA4","Cluster2"="#D7AB25","Cluster3"="#A51F27")
  
)

immune_cell_order <- c("T.cells.CD8","T.cells.regulatory..Tregs.","T.cells.CD4.naive","T.cells.follicular.helper",
                       "B.cells.naive","B.cells.memory",
                       "T.cells.gamma.delta","Dendritic.cells.activated","Macrophages.M1","NK.cells.activated",
                       "Plasma.cells","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated",
                       "Mast.cells.activated","NK.cells.resting",
                       "Macrophages.M0","Monocytes","Eosinophils","Macrophages.M2","Mast.cells.resting",
                       "Dendritic.cells.resting","Neutrophils","Endothelial cells","Fibroblasts"
)



pic_data <-res_dev_all[immune_cell_order,rownames(anno_col)]
rownames(pic_data) <- c("CD8+T cells","Regulatory T cells","Naive CD4+T cells","Follicular helper T cells",
                        "Naive B cells","Memory B cells",
                        "?æ? Tcells","Activated dendritic cells","M1 Macrophages","Activated NK cells",
                        "Plasma cells","Resting CD4+ memory T cells","Activated CD4+ memory T cells",
                        "Activated Mast cells","Resting NK cells",
                        "M0 Macrophages","Monocytes","Eosinophils","M2 Macrophages","Resting Mast cells",
                        "Resting dendritic cells","Neutrophils","Endothelial cells","Fibroblasts")




pheatmap(as.matrix(pic_data), color = mycolors, breaks=bk ,fontsize_row = 8, cluster_rows = F,cluster_cols=F,show_colnames = F, 
         clustering_distance_cols = "euclidean", treeheight_col = 50, scale = "row",
         annotation_col = anno_col,annotation_colors=anno_color,
         filename = "/results/Fig6A_1.pdf")




####-----------------????????????????-Fig 6B-----------------------
#-----ICB_siganture

ICB_siganture_list <-list()
for (i in unique(ICB_siganture$Biomarker)){
  ICB_siganture_list[[i]] <-  c(ICB_siganture[ICB_siganture$Biomarker == i,2])
}
ICB_signature <- data.frame(row.names = colnames(expr_data_TT_LOGtpm))
for (i in 1: length(ICB_siganture_list)) {
  gene_list <- ICB_siganture_list[[i]]$Gene
  expre <- expr_data_TT_LOGtpm[gene_list,]
  expre_scale <- na.omit(expre)
  for (j in 1:nrow(expre_scale)){
    expre_scale[j,] <-scale(as.numeric(expre_scale[j,]),center =T)
  }
  sig_score<- as.numeric(colMeans(expre_scale))
  sig_score_scal <- as.data.frame(scale(sig_score,scale = T,center = T))
  colnames(sig_score_scal) <- names(ICB_siganture_list[i])
  ICB_signature <- cbind(ICB_signature,sig_score_scal)
}

##draw
anno_col_1 <- CBCGA_Cohort.Info[,c(1,10)]
anno_col_1$RNA_seqID_T <- paste(anno_col_1$PatientCode,"RNA_T",sep = "_")
anno_col_1 <- subset(anno_col_1,anno_col_1$RNA_seqID_T%in%rownames(anno_col))
anno_col_1 <- cbind(anno_col_1,anno_col[anno_col_1$RNA_seqID_T,])
anno_col_1<- anno_col_1[,c(3,2,4,5)]
colnames(anno_col_1) <- c("RNAseq_ID","PAM50","Immune_cluster","TILs")
pic_data_ICB_sig <- cbind(anno_col_1,ICB_signature[anno_col_1$RNAseq_ID,])

cols_Immune_cluster_anno <- c('Cluster1'='#6682AA','Cluster2'='#DFAF1E','Cluster3'='#AE1426')
colnames(pic_data_ICB_sig) <- c("RNAseq_ID","PAM50","Immune_cluster","TILs","TIS_signature","STAT1_signature")
 
# p <- ggboxplot(pic_data_ICB_sig, x = "PAM50", y = "STAT1_signature",
#                color = "Immune_cluster", palette = cols_Immune_cluster_anno
# )
# ggsave(plot = p,filename = "/results/Fig6B_ImmunE_PAm50_STAT1.pdf",width = 8,height = 6)
# a <- compare_means(STAT1_signature ~ Immune_cluster, pic_data_ICB_sig, group.by = "PAM50",method = "kruskal.test")
# write.csv(a,file = "/results/Fig6B_STAT1_PAM50_signature_KW.csv")
# 
# 
# p <- ggboxplot(pic_data_ICB_sig, x = "PAM50", y = "TIS_signature",
#                color = "Immune_cluster", palette = cols_Immune_cluster_anno
# )
# ggsave(plot = p,filename = "/results/Fig6B_ImmunE_PAm50_TIS.pdf",width = 8,height = 6)
# a <- compare_means(TIS_signature ~ Immune_cluster, pic_data_ICB_sig, group.by = "PAM50",method = "kruskal.test")
# write.csv(a,file = "/results/Fig6B_TIS_PAM50_signature_KW.csv")


####-----------------????????????????-Fig 6C-----------------------

P <-anno_col_1 %>% ggplot(aes(x = PAM50, fill = Immune_cluster)) +
  geom_bar(position = position_fill()) + 
  scale_fill_manual(values =cols_Immune_cluster_anno) + 
  theme(panel.border = element_rect(fill = NA)) +
  labs(y = 'Percent') +
  labs(x='')+coord_flip()
P
ggsave(plot = P,filename = "/results/Fig6C_ImmunE_PAm50_Pro.pdf",width = 8,height = 4)

####-----------------????????????????-Fig 6D-----------------------

# library(dplyr)
# #library(ggpubr)
# Immunogenic_mean <- data.frame()
# Immunogenic_WT_test <-  data.frame()
# tmp_data <- Immunogenic_data_gene
# tmp_data$PAM50_immun <- paste(tmp_data$PAM50,tmp_data$Immune_cluster,sep = "_")
# colnames(tmp_data)
# for (i in c(5:20)) {
#   data <- tmp_data
#   colnames(data)[i] <- "var"
#   data$var <- as.numeric(data$var)
#   WT_test <- compare_means(var ~ PAM50_immun, data = data,method = "wilcox.test")
#   anova_res <- oneway.test(var ~ PAM50_immun, data = data)
#   WT_test$Anno <- colnames(tmp_data)[i]
#   Immunogenic_WT_test <- rbind(Immunogenic_WT_test,WT_test)
#   OutTab <- data.frame(
#     row.names =  colnames(tmp_data)[i],
#     Basal_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Basal_Cluster1",i]))),
#     Basal_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Basal_Cluster2",i]))),
#     Basal_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Basal_Cluster3",i]))),
#     LumA_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumA_Cluster1",i]))),
#     LumA_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumA_Cluster2",i]))),
#     LumA_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumA_Cluster3",i]))),
#     LumB_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumB_Cluster1",i]))),
#     LumB_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumB_Cluster2",i]))),
#     LumB_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumB_Cluster3",i]))),
#     Her2_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Her2_Cluster1",i]))),
#     Her2_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Her2_Cluster2",i]))),
#     Her2_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Her2_Cluster3",i]))),
#     Normal_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Normal_Cluster1",i]))),
#     Normal_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Normal_Cluster2",i]))),
#     Normal_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Normal_Cluster3",i]))),
#     ANOVA_test_pvalue = anova_res$p.value
#   )
#   Immunogenic_mean <- rbind(Immunogenic_mean,OutTab)
# }
# write.csv(Immunogenic_mean,file = "/results/Fig6D_Immunogenic_mean.csv")
# Immunogenic_mean_pic <- Immunogenic_mean[c("TRA_diversity","TRB_diversity","TRA_clonality","TRB_clonality","IGH_diversity","IGL_diversity","IGH_clonality",
#                                            "IGL_clonality","TMB","HRD","missense","FS","Inframe"),1:15]
# 
# rownames(Immunogenic_mean_pic) <- c("TRA diversity","TRB diversity","TRA clonality","TRB clonality","IGH diversity","IGL diversity",
#                                     "IGH clonality","IGL clonality","TMB","HRD","Missense","Frameshift","Inframe")
# 
# pheatmap(Immunogenic_mean_pic,show_colnames = T,cluster_cols = F,cluster_rows = F,breaks = bk, scale = "row",
#          treeheight_row = 0,fontsize = 12,color =mycolors,
#          filename = "/results/Fig6D_Immunogenic_Immune_PAM50.pdf",width =6,height = 6)
# 
# graphics.off()

####-----------------????????????????-Fig 6E-----------------------
rownames(Immunogenic_data_gene) <- substr(Immunogenic_data_gene$RNAseq_ID,1,4)
patient_filter <- rownames(Immunogenic_data_gene)
patient_filter <- intersect(patient_filter,rownames(HLA_matrix)) 

Cus_Clinicaldata <- Immunogenic_data_gene[patient_filter,c("PAM50","Immune_cluster")]

Cus_HLA <- HLA_matrix[patient_filter,]


combine_matrix <- cbind(Cus_Clinicaldata,Cus_HLA)

subtype_HLA_matrix <- matrix(nrow=10,ncol=15) 
colnames(subtype_HLA_matrix) <- c(paste("LumA",c("Cluster1","Cluster2","Cluster3"),sep="_"),
                                  paste("LumB",c("Cluster1","Cluster2","Cluster3"),sep="_"),
                                  paste("Her2",c("Cluster1","Cluster2","Cluster3"),sep="_"),
                                  paste("Basal",c("Cluster1","Cluster2","Cluster3"),sep="_"),
                                  paste("Normal",c("Cluster1","Cluster2","Cluster3"),sep="_")
)
rownames(subtype_HLA_matrix) <- c("sum_number",
                                  paste(c("HLA-A","HLA-B","HLA-C"),"homogeneity",sep="_"),
                                  paste(c("HLA-A","HLA-B","HLA-C"),"LOH",sep="_"),
                                  paste(c("HLA-A","HLA-B","HLA-C"),"AI",sep="_")
)
subtype_HLA_matrix[1,] <- c(as.numeric(table(Cus_Clinicaldata$Immune_cluster[Cus_Clinicaldata$PAM50=="LumA"])),
                            as.numeric(table(Cus_Clinicaldata$Immune_cluster[Cus_Clinicaldata$PAM50=="LumB"])),
                            as.numeric(table(Cus_Clinicaldata$Immune_cluster[Cus_Clinicaldata$PAM50=="Her2"])),
                            as.numeric(table(Cus_Clinicaldata$Immune_cluster[Cus_Clinicaldata$PAM50=="Basal"])),
                            as.numeric(table(Cus_Clinicaldata$Immune_cluster[Cus_Clinicaldata$PAM50=="Normal"]))
)
for (i in c("LumA","LumB","Her2","Basal","Normal")){
  for (j in c("Cluster1","Cluster2","Cluster3")){
    tmp <- combine_matrix[combine_matrix$PAM50==i & combine_matrix$Immune_cluster==j,]
    a <- c(length(rownames(tmp)[tmp$A_homo=="homo"]),length(rownames(tmp)[tmp$B_homo=="homo"]),length(rownames(tmp)[tmp$C_homo=="homo"]),
           length(rownames(tmp)[tmp$A_LOH=="LOH" & is.na(tmp$A_LOH)==F]),length(rownames(tmp)[tmp$B_LOH=="LOH"& is.na(tmp$B_LOH)==F]),length(rownames(tmp)[tmp$C_LOH=="LOH"& is.na(tmp$C_LOH)==F]),
           length(rownames(tmp)[tmp$A_AI=="AI"& is.na(tmp$A_AI)==F]),length(rownames(tmp)[tmp$B_AI=="AI"& is.na(tmp$B_AI)==F]),length(rownames(tmp)[tmp$C_AI=="AI"& is.na(tmp$C_AI)==F])
    )
    subtype_HLA_matrix[2:10,paste(i,j,sep="_")] <- a/subtype_HLA_matrix[1,paste(i,j,sep="_")]
  }
}

write.csv(subtype_HLA_matrix,file="/results/Fig6E_subtype_HLA_matrix.csv")


####-----------------????????????????-Fig 6F-----------------------
OMTransFun                   <- function(x){ifelse(x == 0, 0, 1)}
library(maftools)
CBCGA.Extended_WES_maftools  <- read.maf(CBCGA.Extended_WES_Somatic, useAll = TRUE)

CBCGA.Extended_WES_Muttable            <- mutCountMatrix(CBCGA.Extended_WES_maftools, includeSyn = F, countOnly = NULL,
                                                         removeNonMutated = FALSE)
CBCGA.Extended_WES_Muttable            <- apply(CBCGA.Extended_WES_Muttable, 2, OMTransFun)
CBCGA.Extended_WES_Muttable[1:5,1:5]

patient_filter <- rownames(Immunogenic_data_gene)
patient_filter <- intersect(patient_filter,colnames(CBCGA.Extended_WES_Muttable))

Cus_clinicaldata <- Immunogenic_data_gene[patient_filter,]
Cus_Muttable <- CBCGA.Extended_WES_Muttable[,patient_filter]

Freq <- apply(Cus_Muttable, 1, sum)/ncol(Cus_Muttable)
gene_filter <- names(Freq)[Freq>0.03]

Cus_Muttable <- Cus_Muttable[gene_filter,]


Cluster_mut_Freq <- matrix(nrow=nrow(Cus_Muttable),ncol=6)
rownames(Cluster_mut_Freq) <- rownames(Cus_Muttable)
colnames(Cluster_mut_Freq) <- c("Cluster1","Cluster2","Cluster3","mean","P","FDR")

Cus_Muttable_Cluster1 <- Cus_Muttable[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster1",])]
Cus_Muttable_Cluster2 <- Cus_Muttable[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster2",])]
Cus_Muttable_Cluster3 <- Cus_Muttable[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster3",])]

Cluster1Freqs      <- apply(Cus_Muttable[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster1",])], 1, sum)/ncol(Cus_Muttable[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster1",])])
Cluster2Freqs      <- apply(Cus_Muttable[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster2",])], 1, sum)/ncol(Cus_Muttable[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster2",])])
Cluster3Freqs      <- apply(Cus_Muttable[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster3",])], 1, sum)/ncol(Cus_Muttable[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster3",])])
ClusterFreqs      <- apply(Cus_Muttable, 1, sum)/ncol(Cus_Muttable)


Cluster_mut_Freq[,1] <- c(Cluster1Freqs)
Cluster_mut_Freq[,2] <- c(Cluster2Freqs)
Cluster_mut_Freq[,3] <- c(Cluster3Freqs)
Cluster_mut_Freq[,4] <- c(ClusterFreqs)

Cluster_mutation_table <- cbind(Cus_clinicaldata[,c("Immune_cluster","PAM50")],t(Cus_Muttable))
colnames(Cluster_mutation_table)[1:2] <- c("Immune_cluster","PAM50")
Cluster_mutation_table <- as.data.frame(Cluster_mutation_table)

# library(glmperm)
# for (i in gene_filter){
#   y <- as.numeric(Cluster_mutation_table[,i])
#   x1 <- as.numeric(as.factor(Cluster_mutation_table[,"PAM50"]))
#   x2 <- as.numeric(as.factor(Cluster_mutation_table[,"Immune_cluster"]))
#   f <- prr.test(y~x2+x1,var="x2",family=gaussian)
#   Cluster_mut_Freq[i,5] <- f$p.value.obs
# }
# Cluster_mut_Freq[,6] <- p.adjust(Cluster_mut_Freq[,5],method="fdr")
# 
# write.table(Cluster_mut_Freq,file="/results/Fig6F_Cluster_mut_Freq_adjustPAm50.txt",sep="\t")

####-----------------????????????????-Fig 6G-----------------------
patient_filter <- rownames(Immunogenic_data_gene)
patient_filter <- intersect(patient_filter,colnames(CBCGA_GISTICpeaks.amp.thre))

Cus_clinicaldata <- Immunogenic_data_gene[patient_filter,]
Cus_CNA_amp <- CBCGA_GISTICpeaks.amp.thre[,patient_filter]
Cus_CNA_del <- CBCGA_GISTICpeaks.del.thre[,patient_filter]

Cus_CNA_amp[as.matrix(which(Cus_CNA_amp==2,arr.ind = TRUE))] <- 1

Cus_CNA_del[as.matrix(which(Cus_CNA_del==2,arr.ind = TRUE))] <- 1

max(Cus_CNA_amp)
max(Cus_CNA_del)

Cus_CNA <- rbind(Cus_CNA_amp,Cus_CNA_del)

Cluster_CNA_Freq <- matrix(nrow=nrow(Cus_CNA),ncol=6)
rownames(Cluster_CNA_Freq) <- rownames(Cus_CNA)
colnames(Cluster_CNA_Freq) <- c("Cluster1","Cluster2","Cluster3","mean","P","FDR")

Cus_CNA_Cluster1 <- Cus_CNA[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster1",])]
Cus_CNA_Cluster2 <- Cus_CNA[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster2",])]
Cus_CNA_Cluster3 <- Cus_CNA[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster3",])]

Cluster1Freqs      <- apply(Cus_CNA[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster1",])], 1, sum)/ncol(Cus_CNA[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster1",])])
Cluster2Freqs      <- apply(Cus_CNA[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster2",])], 1, sum)/ncol(Cus_CNA[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster2",])])
Cluster3Freqs      <- apply(Cus_CNA[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster3",])], 1, sum)/ncol(Cus_CNA[,rownames(Cus_clinicaldata[Cus_clinicaldata$Immune_cluster=="Cluster3",])])
ClusterFreqs      <- apply(Cus_CNA, 1, sum)/ncol(Cus_CNA)


Cluster_CNA_Freq[,1] <- c(Cluster1Freqs)
Cluster_CNA_Freq[,2] <- c(Cluster2Freqs)
Cluster_CNA_Freq[,3] <- c(Cluster3Freqs)
Cluster_CNA_Freq[,4] <- c(ClusterFreqs)

Cluster_CNA_table <- cbind(Cus_clinicaldata[,c("Immune_cluster","PAM50")],t(Cus_CNA))
colnames(Cluster_CNA_table)[1:2] <- c("immune_Cluster","PAM50")
Cluster_CNA_table <- as.data.frame(Cluster_CNA_table)



# library(glmperm)
# for (i in rownames(Cus_CNA)){
#   y <- as.numeric(Cluster_CNA_table[,i])
#   x1 <- as.numeric(as.factor(Cluster_CNA_table[,"PAM50"]))
#   x2 <- as.numeric(as.factor(Cluster_CNA_table[,"immune_Cluster"]))
#   f <- prr.test(y~x2+x1,var="x2",family=gaussian)
#   Cluster_CNA_Freq[i,5] <- f$p.value.obs
# }
# Cluster_CNA_Freq[,6] <- p.adjust(Cluster_CNA_Freq[,5],method="fdr")
# 
# write.table(Cluster_CNA_Freq,file="/results/Fig6G_Cluster_CNA_Freq_adjustPAm50.txt",sep="\t")


####-----------------????????????????-Fig S6A-----------------------

library(tidyverse)
CIBERSORTx_data <- as.data.frame(CIBERSORTx_data)
rownames(CIBERSORTx_data) <- CIBERSORTx_data$Mixture
pic_CIBERSORTx <- CIBERSORTx_data[,-1]
pic_CIBERSORTx <- pic_CIBERSORTx[rownames(anno_col),immune_cell_order[1:22]]
CIBERSORTx_cell_order_ano <- c("CD8+T cells","Regulatory T cells","Naive CD4+T cells","Follicular helper T cells",
                               "Naive B cells","Memory B cells",
                               "?æ? Tcells","Activated dendritic cells","M1 Macrophages","Activated NK cells",
                               "Plasma cells","Resting CD4+ memory T cells","Activated CD4+ memory T cells",
                               "Activated Mast cells","Resting NK cells",
                               "M0 Macrophages","Monocytes","Eosinophils","M2 Macrophages","Resting Mast cells",
                               "Resting dendritic cells","Neutrophils")
colnames(pic_CIBERSORTx) <- CIBERSORTx_cell_order_ano



library(paletteer) 
getPalette = colorRampPalette(paletteer_d("ggthemes::Hue_Circle"))
LM22_color = getPalette(22)


dd1 <- pic_CIBERSORTx %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(cols=2:23,
               names_to= "celltype",
               values_to = "Proportion")
dd1$celltype <- factor(dd1$celltype,levels = CIBERSORTx_cell_order_ano)
dd1$sample <- factor(dd1$sample,levels = rownames(anno_col))
p<- ggplot(dd1,aes(sample,Proportion,fill = celltype)) + 
  geom_bar(position = "stack",stat = "identity")+
  theme_bw()+
  scale_fill_manual(values =LM22_color)+
  theme(axis.text.x = element_text(angle = 90,size = 1))

ggsave(plot = p,filename = "/results/FigS6A_CIBERSORTx_propotion.pdf",width = 16,height = 8)

####-----------------????????????????-Fig S6B-----------------------

# gmtFile <- system.file("extdata", "SI_geneset.gmt",package="estimate")
# signature <- getGmt(gmtFile, geneIdType=SymbolIdentifier())
# res_estimate <- gsva(as.matrix(expr_data_TT_LOGtpm),signature,method='ssgsea',kcdf='Gaussian',abs.ranking=T)
# res_estimate_1 <- as.data.frame(t(res_estimate))
# Estimate_merge <- cbind(anno_col_1,res_estimate_1[anno_col_1$RNAseq_ID,])
# 
# colnames(Estimate_merge) 
# p <- ggboxplot(Estimate_merge, x = "PAM50", y = "StromalSignature",
#                color = "Immune_cluster", palette = cols_Immune_cluster_anno
# )
# p
# ggsave(plot = p,filename = "/results/FigS6B_ImmunE_PAm50_StromalSignature.pdf",width = 8,height = 6)
# a <- compare_means(StromalSignature ~ Immune_cluster, Estimate_merge, group.by = "PAM50",method = "kruskal.test")
# write.csv(a,file = "/results/FigS6B_StromalSignature_PAM50_signature_KW.csv")
# 
# 
# p <- ggboxplot(Estimate_merge, x = "PAM50", y = "ImmuneSignature",
#                color = "Immune_cluster", palette = cols_Immune_cluster_anno
# )
# p
# ggsave(plot = p,filename = "/results/FigS6B_ImmunE_PAm50_ImmuneSignature.pdf",width = 8,height = 6)
# a <- compare_means(ImmuneSignature ~ Immune_cluster, Estimate_merge, group.by = "PAM50",method = "kruskal.test")
# write.csv(a,file = "/results/FigS6B_ImmuneSignature_PAM50_signature_KW.csv")

####-----------------????????????????-Fig S6C-----------------------


####-----------------????????????????-Fig S6D-----------------------


####-----------------????????????????-Fig S6E-----------------------

library(tidyverse)
scCIBERSORTx_data <- as.data.frame(scCIBERSORTx_data)
rownames(scCIBERSORTx_data) <- scCIBERSORTx_data$Mixture
pic_scCIBERSORTx <- scCIBERSORTx_data[,-1]
scCIBERSORTx_cell_order_ano <- c("CD8TZFP36","CD8TGZMK","CD8TIFNG","CD8TLAG3","TIFIT1","Treg","CD4TCCR7","Tfh",
                                 "Bn","Bmem", "Plasmablasts",
                                 "cDC2","pDC","cDC1","DC","Macrophage2",#M1 like
                                 "NKT","NK","CD4TIL7R",
                                 "Macrophage1","Macrophage3","LAM1","LAM2",   
                                 "MonoS100A9","MonoFCGR3A","MonoIL1B",
                                 "Myoepithelial","ProLum","MatureLum",    
                                 "HerCancer","LumBCancer","BasalCancer","LumACancer",
                                 "Cycling",
                                 "EndoACKR1","EndoRGS5","EndoCXCL12","EndoLymph",
                                 "PVL3", "PVL2","PVL1",
                                 "iCAF1","iCAF2","tCAF3", "myCAF4","myCAF5")

getPalette = colorRampPalette(paletteer_d("ggthemes::Hue_Circle"))
sc46_color = getPalette(46)
getPalette = colorRampPalette(paletteer_d("colorBlindness::PairedColor12Steps"))
sc46_color = getPalette(46)

dd1 <- pic_scCIBERSORTx %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(cols=2:47,
               names_to= "celltype",
               values_to = "Proportion")
dd1$celltype <- factor(dd1$celltype,levels = scCIBERSORTx_cell_order_ano)
dd1$sample <- factor(dd1$sample,levels = rownames(anno_col))
p<- ggplot(dd1,aes(sample,Proportion,fill = celltype)) + 
  geom_bar(position = "stack",stat = "identity")+
  theme_bw()+
  scale_fill_manual(values =sc46_color)+
  theme(axis.text.x = element_text(angle = 90,size = 1))
p
ggsave(plot = p,filename = "/results/FigS6E_scCIBERSORTx_propotion_2.pdf",width = 16,height = 8)




####-----------------????????????????-Fig S6F-----------------------
# MHC_data <- c("HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M",
#               "HLA-DPA1","HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DRB1","HLA-DRB5","HLA-DRB6",
#               "HLA-E","HLA-F","HLA-G","HLA-H")
# MHC_DATA_RNA <- expr_data_TT_LOGtpm[MHC_data,anno_col_1$RNAseq_ID]
# MHC_DATA_RNA <- cbind(t(MHC_DATA_RNA),anno_col_1)
# MHC_DATA_RNA$PAM50_immun <- paste(MHC_DATA_RNA$PAM50,MHC_DATA_RNA$Immune_cluster,sep = "_")
# MHC_DATA_RNA_mean <- data.frame()
# 
# resu_KW <- data.frame()
# resu_WT <- data.frame()
# for (i in c(1:length(MHC_data))) {
#   tmp <- MHC_DATA_RNA
#   colnames(tmp)[i]<- "gene"
#   cluster_WT_test <- compare_means(gene~Immune_cluster,tmp)
#   cluster_WT_test$gene <- colnames(MHC_DATA_RNA)[i]
#   cluster_KW_test <- compare_means(gene~Immune_cluster,tmp,method = "kruskal.test")
#   cluster_KW_test$gene <- colnames(MHC_DATA_RNA)[i]
#   resu_KW <- rbind(resu_KW,cluster_KW_test)
#   resu_WT <- rbind(resu_WT,cluster_WT_test)}
# write.csv(resu_KW,file = "/results/FigS6F_MHC_DATA_RNA_resu_KW.csv")
# write.csv(resu_WT,file = "/results/FigS6F_MHC_DATA_RNA_resu_WT.csv")
# 
# 
# for (i in c(1:length(MHC_data))) {
#   data <- MHC_DATA_RNA
#   OutTab <- data.frame(
#     row.names =  colnames(MHC_DATA_RNA)[i],
#     Basal_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Basal_Cluster1",i]))),
#     Basal_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Basal_Cluster2",i]))),
#     Basal_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Basal_Cluster3",i]))),
#     LumA_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumA_Cluster1",i]))),
#     LumA_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumA_Cluster2",i]))),
#     LumA_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumA_Cluster3",i]))),
#     LumB_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumB_Cluster1",i]))),
#     LumB_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumB_Cluster2",i]))),
#     LumB_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumB_Cluster3",i]))),
#     Her2_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Her2_Cluster1",i]))),
#     Her2_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Her2_Cluster2",i]))),
#     Her2_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Her2_Cluster3",i]))),
#     Normal_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Normal_Cluster1",i]))),
#     Normal_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Normal_Cluster2",i]))),
#     Normal_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Normal_Cluster3",i])))  )
#   MHC_DATA_RNA_mean <- rbind(MHC_DATA_RNA_mean,OutTab)
# }
# write.csv(MHC_DATA_RNA_mean,file = "/results/FigS6F_MHC_DATA_RNA_mean.csv")
# 
# pheatmap(MHC_DATA_RNA_mean,show_colnames = T,cluster_cols = F,cluster_rows = F,breaks = bk, scale = "row",
#          treeheight_row = 0,fontsize = 12,color =mycolors
#          ,filename = "/results/FigS6F_MHC_Immune_PAM50.pdf",width =6,height = 6)
# 
# graphics.off()


####-----------------????????????????-Fig S6G-----------------------

Innate_immune <- c("CGAS","TMEM173","IRF3",
                   "MYD88","TICAM1","TLR1","TLR10","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR9",
                   "DDX58","IFIH1","MAVS",
                   "CLEC7A","CLEC6A","CLEC4E","CD209","CLEC10A",
                   "NLRP3","NLRP6","NLRP12","AIM2","PYCARD")
Innate_immune_DATA_RNA <- expr_data_TT_LOGtpm[Innate_immune,anno_col_1$RNAseq_ID]
Innate_immune_DATA_RNA <- cbind(t(Innate_immune_DATA_RNA),anno_col_1)
Innate_immune_DATA_RNA$PAM50_immun <- paste(Innate_immune_DATA_RNA$PAM50,Innate_immune_DATA_RNA$Immune_cluster,sep = "_")
Innate_immune_DATA_RNA_mean <- data.frame()
resu_KW <- data.frame()
resu_WT <- data.frame()
# for (i in c(1:length(Innate_immune))) {
#   tmp <- Innate_immune_DATA_RNA
#   colnames(tmp)[i]<- "gene"
#   cluster_WT_test <- compare_means(gene~Immune_cluster,tmp)
#   cluster_WT_test$gene <- colnames(Innate_immune_DATA_RNA)[i]
#   cluster_KW_test <- compare_means(gene~Immune_cluster,tmp,method = "kruskal.test")
#   cluster_KW_test$gene <- colnames(Innate_immune_DATA_RNA)[i]
#   resu_KW <- rbind(resu_KW,cluster_KW_test)
#   resu_WT <- rbind(resu_WT,cluster_WT_test)}
# write.csv(resu_KW,file = "/results/FigS6G_Innate_immune_RNA_resu_KW.csv")
# write.csv(resu_WT,file = "/results/FigS6G_Innate_immune_RNA_resu_WT.csv")
# 
# for (i in c(1:length(Innate_immune))) {
#   data <- Innate_immune_DATA_RNA
#   OutTab <- data.frame(
#     row.names =  colnames(Innate_immune_DATA_RNA)[i],
#     Basal_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Basal_Cluster1",i]))),
#     Basal_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Basal_Cluster2",i]))),
#     Basal_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Basal_Cluster3",i]))),
#     LumA_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumA_Cluster1",i]))),
#     LumA_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumA_Cluster2",i]))),
#     LumA_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumA_Cluster3",i]))),
#     LumB_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumB_Cluster1",i]))),
#     LumB_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumB_Cluster2",i]))),
#     LumB_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="LumB_Cluster3",i]))),
#     Her2_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Her2_Cluster1",i]))),
#     Her2_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Her2_Cluster2",i]))),
#     Her2_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Her2_Cluster3",i]))),
#     Normal_1 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Normal_Cluster1",i]))),
#     Normal_2 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Normal_Cluster2",i]))),
#     Normal_3 = mean(na.omit(as.matrix(data[data$PAM50_immun=="Normal_Cluster3",i])))  )
#   Innate_immune_DATA_RNA_mean <- rbind(Innate_immune_DATA_RNA_mean,OutTab)
# }
# write.csv(Innate_immune_DATA_RNA_mean,file = "/results/FigS6G_Innate_immune_DATA_RNA_mean.csv")
# 
# 
# 
# pheatmap(Innate_immune_DATA_RNA_mean,show_colnames = T,cluster_cols = F,cluster_rows = F,breaks = bk, scale = "row",
#          treeheight_row = 0,fontsize = 12,color =mycolors
#          ,filename = "/results/FigS6G_Innate_Immune_PAM50.pdf",width =6,height = 6)
# 
# graphics.off()



####-----------------????????????????-Fig S6H-----------------------
virus_gene<- c("ADAR","BST2","GBP2","IFI35","IFIT1","IFIT2","IFIT3","IFITM3",'IFNB1','IRF1','IRF7','IRF9','ISG15','ISG20','MX1','MX2',
               'NLRC5','OAS2','OAS3','PSMB8','RSAD2','SAMHD1','SP100','STAT1','STAT2',"TREX1",'USP18','XAF1','ZBP1')
virus_mi_DATA_RNA <- t(expr_data_TT_LOGtpm[virus_gene,anno_col_1$RNAseq_ID])
result_zscore <- data.frame()
for (k in 1:length(virus_gene)){
  data <- as.numeric(virus_mi_DATA_RNA[,colnames(virus_mi_DATA_RNA)[k]])
  z_scores <- (data-mean(data))/sd(data)
  result_zscore <- rbind(result_zscore,z_scores)
}
rownames(result_zscore)<- colnames(virus_mi_DATA_RNA)[1:29]
colnames(result_zscore)<- rownames(virus_mi_DATA_RNA)
virus_sig <- colMeans(result_zscore)
names(virus_sig) <- colnames(result_zscore)
table(names(virus_sig)==rownames(anno_col_1))
anno_col_1$virus_mimicy <- virus_sig
# p <- ggboxplot(anno_col_1, x = "PAM50", y = "virus_mimicy",
#                color = "Immune_cluster", palette = cols_Immune_cluster_anno
# )
# p
# ggsave(plot = p,filename = "/results/FigS6H_ImmunE_PAm50_virus_mimicy.pdf",width = 8,height = 6)
# a <- compare_means(virus_mimicy ~ Immune_cluster, anno_col_1, group.by = "PAM50",method = "kruskal.test")
# write.csv(a,file = "/results/FigS6H_virus_mimicy_PAM50_signature_KW.csv")
