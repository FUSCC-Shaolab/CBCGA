#-----------------------------------------------------------------#
# Fig5 
#-----------------------------------------------------------------#
rm(list = ls()) ; graphics.off()
options(stringsAsFactors = F)
set.seed(123)


library(Rtsne)
#library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyverse)
####Load data
load("/data/CBCGA.Extended_MergedData_V2.5_220722.Rdata")
load("/data/Fig5_S5_data.Rdata")

length(unique(metabolic_gene_list$b)) #4134
mypro <- intersect(unique(metabolic_gene_list$b),rownames(CBCGA.Extended_PRO_normalized)) #2252
#####################################################################################################################
#####################################################################################################################
####------------------Fig 5B--------------####
#####################################################################################################################
#####################################################################################################################
T_N_label <- as.data.frame(matrix(ncol=2,nrow=ncol(CBCGA.Extended_pol)))
colnames(T_N_label) <- c("ID","T_N")
T_N_label[,1] <- colnames(CBCGA.Extended_pol)
T_N_label[,2] <- substr(T_N_label[,1],10,10)
T_N_label$color <- ifelse (T_N_label$T_N=="T","yellow","green")

T_N_pol <- scale(t(CBCGA.Extended_pol))
#for (i in 1:5){
i = 5   
set.seed(12345)
tsne1 <- Rtsne(T_N_pol, dimS=3, perplexity=i*10, theta=0.0, verbose=TRUE, max_iter=500)
pdf("/results/Fig5B_1.pdf")
plot(tsne1$Y, pch=16,col=T_N_label$color ,main="tsne",cex=1.5)
dev.off()
#}

T_N_lip <- scale(t(CBCGA.Extended_lip))
#for (i in 1:5){
i = 3
set.seed(12345)
tsne1 <- Rtsne(T_N_lip, dimS=3, perplexity=i*10, theta=0.0, verbose=TRUE, max_iter=500)
pdf("/results/Fig5B_2.pdf")
plot(tsne1$Y, pch=16,col=T_N_label$color ,main="tsne",cex=1.5)
dev.off()
#}

#####################################################################################################################
#####################################################################################################################
####------------------Fig 5C-D--------------####
#####################################################################################################################
#####################################################################################################################
## 1. Data filtering and cleaning
#####################################################################################################################

CBCGA_ID <- rownames(CBCGA.Extended_Cohort.Info[CBCGA.Extended_Cohort.Info$CBCGA_Cohort == "Yes",]) # 773

myclinicaldata <- Add_cohort.info[,c(1,5)]
colnames(myclinicaldata) <- c("PatientCode","PAM50_classifier")
myclinicaldata <- rbind(myclinicaldata,CBCGA.Extended_Cohort.Info[CBCGA_ID,c("PatientCode","PAM50_classifier")])
rownames(myclinicaldata) <- myclinicaldata$PatientCode

Basal_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "Basal",1] # 173
Her2_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "Her2",1] # 181
LumA_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "LumA",1] # 244
LumB_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "LumB",1] # 242
Normal_ID <- myclinicaldata[myclinicaldata$PAM50_classifier == "Normal",1] # 75

Basal_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(Basal_ID,"_PRO_T",sep=""))]) # 59
Her2_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(Her2_ID,"_PRO_T",sep=""))]) # 59
LumA_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(LumA_ID,"_PRO_T",sep=""))]) # 56
LumB_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(LumB_ID,"_PRO_T",sep=""))]) # 77
Normal_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(Normal_ID,"_PRO_T",sep=""))]) # 20
PT_PRO <- as.data.frame(CBCGA.Extended_PRO_normalized[,intersect(colnames(CBCGA.Extended_PRO_normalized),paste0(myclinicaldata$PatientCode,"_PRO_N",sep=""))]) # 86

Basal_PRO <- Basal_PRO[mypro,] # 59
colnames(Basal_PRO) <- substr(colnames(Basal_PRO),1,4)
Her2_PRO <- Her2_PRO[mypro,] # 59
colnames(Her2_PRO) <- substr(colnames(Her2_PRO),1,4)
LumA_PRO <- LumA_PRO[mypro,] # 56
colnames(LumA_PRO) <- substr(colnames(LumA_PRO),1,4)
LumB_PRO <- LumB_PRO[mypro,] # 77
colnames(LumB_PRO) <- substr(colnames(LumB_PRO),1,4)
Normal_PRO <- Normal_PRO[mypro,] # 20
colnames(Normal_PRO) <- substr(colnames(Normal_PRO),1,4)

TT_PRO <- cbind(Basal_PRO,Her2_PRO,LumA_PRO,LumB_PRO,Normal_PRO) # 271
PT_PRO <- PT_PRO[mypro,] # 86
colnames(PT_PRO) <- substr(colnames(PT_PRO),1,4)


Basal_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(Basal_ID,"_pol_T",sep=""))] # 92
colnames(Basal_pol) <- substr(colnames(Basal_pol),1,4)
Her2_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(Her2_ID,"_pol_T",sep=""))] # 110
colnames(Her2_pol) <- substr(colnames(Her2_pol),1,4)
LumA_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(LumA_ID,"_pol_T",sep=""))] # 120
colnames(LumA_pol) <- substr(colnames(LumA_pol),1,4)
LumB_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(LumB_ID,"_pol_T",sep=""))] # 144
colnames(LumB_pol) <- substr(colnames(LumB_pol),1,4)
Normal_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste0(Normal_ID,"_pol_T",sep=""))] # 35
colnames(Normal_pol) <- substr(colnames(Normal_pol),1,4)

TT_pol <- cbind(Basal_pol,Her2_pol,LumA_pol,LumB_pol,Normal_pol) # 501
PT_pol <- CBCGA.Extended_pol[,intersect(colnames(CBCGA.Extended_pol),paste(myclinicaldata$PatientCode,"_pol_N",sep=""))] # 76
colnames(PT_pol) <- substr(colnames(PT_pol),1,4)


Basal_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(Basal_ID,"_pol_T",sep=""))] # 92
colnames(Basal_lip) <- substr(colnames(Basal_lip),1,4)
Her2_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(Her2_ID,"_pol_T",sep=""))] # 110
colnames(Her2_lip) <- substr(colnames(Her2_lip),1,4)
LumA_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(LumA_ID,"_pol_T",sep=""))] # 120
colnames(LumA_lip) <- substr(colnames(LumA_lip),1,4)
LumB_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(LumB_ID,"_pol_T",sep=""))] # 144
colnames(LumB_lip) <- substr(colnames(LumB_lip),1,4)
Normal_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste0(Normal_ID,"_pol_T",sep=""))] # 35
colnames(Normal_lip) <- substr(colnames(Normal_lip),1,4)

TT_lip <- cbind(Basal_lip,Her2_lip,LumA_lip,LumB_lip,Normal_lip) # 501
PT_lip <- CBCGA.Extended_lip[,intersect(colnames(CBCGA.Extended_lip),paste(myclinicaldata$PatientCode,"_pol_N",sep=""))] # 76
colnames(PT_lip) <- substr(colnames(PT_lip),1,4)


#####################################################################################################################
## 2. Metabolic protein and polar metabolite subtype-specific analysis
#####################################################################################################################

##### metabolic protein #####
comparison_matrix <- matrix(ncol=14,nrow=length(mypro))
colnames(comparison_matrix) <- c("LumA","LumB","Her2","Basal","Normal","max_min",
                                 "SD","P","FDR","Tumor","PT","T_PT","P_TN","FDR_TN")
rownames(comparison_matrix) <- mypro

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"LumA"] <- mean(as.numeric(LumA_PRO[i,]),na.rm = T)
  comparison_matrix[i,"LumB"] <- mean(as.numeric(LumB_PRO[i,]),na.rm = T)
  comparison_matrix[i,"Her2"] <- mean(as.numeric(Her2_PRO[i,]),na.rm = T)
  comparison_matrix[i,"Basal"] <- mean(as.numeric(Basal_PRO[i,]),na.rm = T)
  comparison_matrix[i,"Normal"] <- mean(as.numeric(Normal_PRO[i,]),na.rm = T)
  comparison_matrix[i,"max_min"] <- max(as.numeric(comparison_matrix[i,1:5]))-min(as.numeric(comparison_matrix[i,1:5]))
  comparison_matrix[i,"SD"] <- sd(TT_PRO[i,],na.rm = T)
  a <- kruskal.test(list(as.numeric(LumA_PRO[i,]),as.numeric(LumB_PRO[i,]),as.numeric(Her2_PRO[i,]),
                         as.numeric(Basal_PRO[i,]),as.numeric(Normal_PRO[i,])))
  comparison_matrix[i,"P"] <- a$p.value
  comparison_matrix[i,"Tumor"] <- mean(as.numeric(TT_PRO[i,]),na.rm = T)
  comparison_matrix[i,"PT"] <- mean(as.numeric(PT_PRO[i,]),na.rm = T)
  comparison_matrix[i,"T_PT"] <- as.numeric(comparison_matrix[i,"Tumor"])-as.numeric(comparison_matrix[i,"PT"])
  b <- wilcox.test(as.numeric(TT_PRO[i,]),as.numeric(PT_PRO[i,]))
  comparison_matrix[i,"P_TN"] <- b$p.value
}
comparison_matrix[,"FDR"] <- p.adjust(comparison_matrix[,"P"],method="fdr")
comparison_matrix[,"FDR_TN"] <- p.adjust(comparison_matrix[,"P_TN"],method="fdr")
comparison_matrix <- as.data.frame(comparison_matrix,stringsAsFactors = F)
CBCGA_PRO_PAM50 <- comparison_matrix
# write.csv(CBCGA_PRO_PAM50,file="CBCGA_PRO_PAM50.csv")


##### Polar metabolite #####
comparison_matrix <- matrix(ncol=14,nrow=nrow(CBCGA.Extended_pol))
colnames(comparison_matrix) <- c("LumA","LumB","Her2","Basal","Normal","max_min",
                                 "SD","P","FDR","Tumor","PT","T_PT","P_TN","FDR_TN")
rownames(comparison_matrix) <- rownames(CBCGA.Extended_pol)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"LumA"] <- mean(as.numeric(LumA_pol[i,]),na.rm = T)
  comparison_matrix[i,"LumB"] <- mean(as.numeric(LumB_pol[i,]),na.rm = T)
  comparison_matrix[i,"Her2"] <- mean(as.numeric(Her2_pol[i,]),na.rm = T)
  comparison_matrix[i,"Basal"] <- mean(as.numeric(Basal_pol[i,]),na.rm = T)
  comparison_matrix[i,"Normal"] <- mean(as.numeric(Normal_pol[i,]),na.rm = T)
  comparison_matrix[i,"max_min"] <- max(as.numeric(comparison_matrix[i,1:5]))-min(as.numeric(comparison_matrix[i,1:5]))
  comparison_matrix[i,"SD"] <- sd(TT_pol[i,],na.rm = T)
  a <- kruskal.test(list(as.numeric(LumA_pol[i,]),as.numeric(LumB_pol[i,]),as.numeric(Her2_pol[i,]),
                         as.numeric(Basal_pol[i,]),as.numeric(Normal_pol[i,])))
  comparison_matrix[i,"P"] <- a$p.value
  comparison_matrix[i,"Tumor"] <- mean(as.numeric(TT_pol[i,]),na.rm = T)
  comparison_matrix[i,"PT"] <- mean(as.numeric(PT_pol[i,]),na.rm = T)
  comparison_matrix[i,"T_PT"] <- as.numeric(comparison_matrix[i,"Tumor"])-as.numeric(comparison_matrix[i,"PT"])
  b <- wilcox.test(as.numeric(TT_pol[i,]),as.numeric(PT_pol[i,]))
  comparison_matrix[i,"P_TN"] <- b$p.value
}
comparison_matrix[,"FDR"] <- p.adjust(comparison_matrix[,"P"],method="fdr")
comparison_matrix[,"FDR_TN"] <- p.adjust(comparison_matrix[,"P_TN"],method="fdr")
comparison_matrix <- as.data.frame(comparison_matrix,stringsAsFactors = F)
CBCGA_Pol_PAM50 <- comparison_matrix
# write.csv(CBCGA_Pol_PAM50,file="CBCGA_Pol_PAM50.csv")


##### Lipid #####
lipid_cat <- unique(CBCGA_lip_anno$Subclass)
Cus_lipid_cat_TT <- matrix(ncol=ncol(TT_lip),nrow=length(lipid_cat))
rownames(Cus_lipid_cat_TT) <- lipid_cat
colnames(Cus_lipid_cat_TT) <- colnames(TT_lip)

for (i in rownames(Cus_lipid_cat_TT)){
  peak <- rownames(CBCGA_lip_anno)[CBCGA_lip_anno$Subclass==i]
  Cus_lipid_cat_TT[i,] <- c(apply(TT_lip[peak,],2,mean))
}

LumA_matrix <- Cus_lipid_cat_TT[,intersect(substr(colnames(CBCGA.Extended_lip),1,4),LumA_ID)] # 120
LumB_matrix <- Cus_lipid_cat_TT[,intersect(substr(colnames(CBCGA.Extended_lip),1,4),LumB_ID)] # 144
Her2_matrix <- Cus_lipid_cat_TT[,intersect(substr(colnames(CBCGA.Extended_lip),1,4),Her2_ID)] # 110
Basal_matrix <- Cus_lipid_cat_TT[,intersect(substr(colnames(CBCGA.Extended_lip),1,4),Basal_ID)] # 92
Normal_matrix <- Cus_lipid_cat_TT[,intersect(substr(colnames(CBCGA.Extended_lip),1,4),Normal_ID)] # 35

comparison_matrix <- matrix(ncol=9,nrow=nrow(Cus_lipid_cat_TT))
colnames(comparison_matrix) <- c("LumA","LumB","Her2","Basal","Normal","max_min",
                                 "SD","P","FDR")
rownames(comparison_matrix) <- rownames(Cus_lipid_cat_TT)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"LumA"] <- mean(as.numeric(LumA_matrix[i,]),na.rm = T)
  comparison_matrix[i,"LumB"] <- mean(as.numeric(LumB_matrix[i,]),na.rm = T)
  comparison_matrix[i,"Her2"] <- mean(as.numeric(Her2_matrix[i,]),na.rm = T)
  comparison_matrix[i,"Basal"] <- mean(as.numeric(Basal_matrix[i,]),na.rm = T)
  comparison_matrix[i,"Normal"] <- mean(as.numeric(Normal_matrix[i,]),na.rm = T)
  comparison_matrix[i,"max_min"] <- max(as.numeric(comparison_matrix[i,1:5]))-min(as.numeric(comparison_matrix[i,1:5]))
  comparison_matrix[i,"SD"] <- sd(Cus_lipid_cat_TT[i,],na.rm = T)
  a <- kruskal.test(list(as.numeric(LumA_matrix[i,]),as.numeric(LumB_matrix[i,]),
                         as.numeric(Her2_matrix[i,]),as.numeric(Basal_matrix[i,]),
                         as.numeric(Normal_matrix[i,])))
  comparison_matrix[i,"P"] <- a$p.value
}
comparison_matrix[,"FDR"] <- p.adjust(comparison_matrix[,"P"],method="fdr")
comparison_matrix <- as.data.frame(comparison_matrix,stringsAsFactors = F)
CBCGA_lip_cat_PAM50 <- comparison_matrix
# write.csv(CBCGA_lip_cat_PAM50,file="CBCGA_lip_cat_PAM50.csv")


#####################################################################################################################
## 3. Metabolic protein and metabolite correlation network construction
#####################################################################################################################

# It takes several time to finish part 3/4 and their figures are not drawn by R....... temporarily skip them 

# ## Metabolic protein correlation network construction
# pro_Cor.Res <- matrix(nrow=nrow(TT_PRO),ncol=nrow(TT_PRO))
# pro_P.val <- matrix(nrow=nrow(TT_PRO),ncol=nrow(TT_PRO))
# colnames(pro_Cor.Res) <- rownames(TT_PRO)
# rownames(pro_Cor.Res) <- rownames(TT_PRO)
# colnames(pro_P.val) <- rownames(TT_PRO)
# rownames(pro_P.val) <- rownames(TT_PRO)
# 
# TT_PRO <- t(TT_PRO)
# 
# #pb <- txtProgressBar(style=3)
# for (i in 1:ncol(TT_PRO)){
#   for (j in 1:ncol(TT_PRO)){
#     TEMP_Res <- cor.test(TT_PRO[,i],TT_PRO[,j], method = "spearman")
#     pro_Cor.Res[i,j] <- TEMP_Res$estimate
#     pro_P.val[i,j] <- TEMP_Res$p.value
#   }
#   #setTxtProgressBar(pb, i/ncol(TT_PRO))
# }
# 
# pro_FDR <- matrix(p.adjust(pro_P.val,method="fdr"),ncol=2252)
# colnames(pro_FDR) <- colnames(pro_P.val)
# rownames(pro_FDR) <- rownames(pro_P.val)
# 
# #save(pro_Cor.Res,pro_P.val,pro_FDR,file="pro_Cor_spearman.Rdata")
# 
# 
# ## Polar metabolite correlation network construction
# pol_Cor.Res <- matrix(nrow=nrow(TT_pol),ncol=nrow(TT_pol))
# pol_P.val <- matrix(nrow=nrow(TT_pol),ncol=nrow(TT_pol))
# colnames(pol_Cor.Res) <- rownames(TT_pol)
# rownames(pol_Cor.Res) <- rownames(TT_pol)
# colnames(pol_P.val) <- rownames(TT_pol)
# rownames(pol_P.val) <- rownames(TT_pol)
# 
# TT_pol <- t(TT_pol)
# 
# #pb <- txtProgressBar(style=3)
# for (i in 1:ncol(TT_pol)){
#   for (j in 1:ncol(TT_pol)){
#     TEMP_Res <- cor.test(TT_pol[,i],TT_pol[,j], method = "spearman")
#     pol_Cor.Res[i,j] <- TEMP_Res$estimate
#     pol_P.val[i,j] <- TEMP_Res$p.value
#   }
#   #setTxtProgressBar(pb, i/ncol(TT_pol))
# }
# 
# pol_FDR <- matrix(p.adjust(pol_P.val,method="fdr"),ncol=669)
# colnames(pol_FDR) <- colnames(pol_P.val)
# rownames(pol_FDR) <- rownames(pol_P.val)
# 
# #save(pol_Cor.Res,pol_P.val,pol_FDR,file="pol_Cor_spearman.Rdata")
# 
# 
# ## Metabolic protein and Polar metabolite correlation network construction
# CBCGA_PRO_POL_ID <- intersect(rownames(TT_PRO),rownames(TT_pol)) # 236
# 
# pro_pol_matrix <- cbind(TT_PRO[CBCGA_PRO_POL_ID,],TT_pol[CBCGA_PRO_POL_ID,]) #236*2921
# dim(pro_pol_matrix)
# 
# pro_pol_Cor.Res <- matrix(nrow=ncol(pro_pol_matrix),ncol=ncol(pro_pol_matrix))
# pro_pol_P.val <- matrix(nrow=ncol(pro_pol_matrix),ncol=ncol(pro_pol_matrix))
# colnames(pro_pol_Cor.Res) <- colnames(pro_pol_matrix)
# rownames(pro_pol_Cor.Res) <- colnames(pro_pol_matrix)
# colnames(pro_pol_P.val) <- colnames(pro_pol_matrix)
# rownames(pro_pol_P.val) <- colnames(pro_pol_matrix)
# 
# #pb <- txtProgressBar(style=3)
# for (i in 1:ncol(pro_pol_matrix)){
#   for (j in 1:ncol(pro_pol_matrix)){
#     TEMP_Res <- cor.test(pro_pol_matrix[,i],pro_pol_matrix[,j], method = "spearman")
#     pro_pol_Cor.Res[i,j] <- TEMP_Res$estimate
#     pro_pol_P.val[i,j] <- TEMP_Res$p.value
#   }
#   #setTxtProgressBar(pb, i/ncol(pro_pol_matrix))
# }
# 
# pro_pol_FDR <- matrix(p.adjust(pro_pol_P.val,method="fdr"),ncol=2921)
# colnames(pro_pol_FDR) <- colnames(pro_pol_P.val)
# rownames(pro_pol_FDR) <- rownames(pro_pol_P.val)
# 
# #save(pro_pol_Cor.Res,pro_pol_P.val,pro_pol_FDR,file="pro_pol_Cor_spearman.Rdata")
# 
# 
# ## Refer to Gephi


#####################################################################################################################
## 4. Several important metabolic pathways presentation
#####################################################################################################################

# key_gene <- c("SLC3A2","IDO1","KYNU","GLS","ACY1","CD38","NNMT","AOX1","CERS4","CERS6","ASAH1")
# key_pol <- c("M205T282_POS","M209T283_POS","M146T478_NEG","M188T410_NEG","M133T521_POS","M664T457_POS","M123T60_POS","M300T73_POS","M302T43_POS")
# 
# TableS4F <- matrix(ncol = 19,nrow = 24)
# colnames(TableS4F) <- c("Protein","Polar_metabolite","Metabolite_name","Cor","P.value",
#                         "LumA_polar","LumB_polar","Her2_polar","Basal_polar","Normal_polar","P.value","FDR",
#                         "LumA_Pro","LumB_Pro","Her2_Pro","Basal_Pro","Normal_Pro","P.value","FDR")
# TableS4F <- as.data.frame(TableS4F)
# TableS4F[,1] <- c("SLC3A2","SLC3A2","IDO1","IDO1","KYNU","KYNU","GLS","GLS","GLS","ACY1","ACY1","ACY1","CD38","CD38","NNMT","NNMT","AOX1","AOX1","CERS4","CERS4","CERS6","CERS6","ASAH1","ASAH1")
# TableS4F[,2] <- c("M205T282_POS","M209T283_POS","M205T282_POS","M209T283_POS","M205T282_POS","M209T283_POS","M146T478_NEG","M188T410_NEG","M133T521_POS","M146T478_NEG","M188T410_NEG","M133T521_POS",
#                   "M664T457_POS","M123T60_POS","M664T457_POS","M123T60_POS","M664T457_POS","M123T60_POS","M300T73_POS","M302T43_POS","M300T73_POS","M302T43_POS","M300T73_POS","M302T43_POS")
# TableS4F[,3] <- CBCGA_pol_anno[TableS4F[,2],"Putative_metabolite_name"]
# 
# for ( i in 1:nrow(TableS4F)){
#   TableS4F[i,4] <- pro_pol_Cor.Res[TableS4F[i,1],TableS4F[i,2]]
#   TableS4F[i,5] <- pro_pol_FDR[TableS4F[i,1],TableS4F[i,2]]
# }
# TableS4F[,6:10] <- CBCGA_Pol_PAM50[TableS4F[,2],1:5]
# TableS4F[,11:12] <- CBCGA_Pol_PAM50[TableS4F[,2],8:9]
# TableS4F[,13:17] <- CBCGA_PRO_PAM50[TableS4F[,1],1:5]
# TableS4F[,18:19] <- CBCGA_PRO_PAM50[TableS4F[,1],8:9]
# write.csv(TableS4F,"/results/TableS4F.csv")

#####################################################################################################################
#####################################################################################################################
####------------------Fig 5E-G--------------####
#####################################################################################################################
#####################################################################################################################

###Fig 5E
#Add_cohort.info <- read.csv(file.choose())
#Add_cohort.info <- read.csv("/results/TableS4.csv")
#colnames(Add_cohort.info) <- Add_cohort.info[1,]
#Add_cohort.info <- Add_cohort.info[-1,]
CBCGA_ID <- rownames(CBCGA.Extended_Cohort.Info[CBCGA.Extended_Cohort.Info$CBCGA_Cohort == "Yes",]) # 773

myclinicaldata <- Add_cohort.info[,c(1,5)]
colnames(myclinicaldata) <- c("PatientCode","PAM50_classifier")
myclinicaldata <- rbind(myclinicaldata,CBCGA.Extended_Cohort.Info[CBCGA_ID,c("PatientCode","PAM50_classifier")])
rownames(myclinicaldata) <- myclinicaldata$PatientCode

CBCGA_Cohort_META <- myclinicaldata
row.names(CBCGA_Cohort_META) <- paste0(myclinicaldata$PatientCode,"_pol_T","")
lipid_data_supplementary <- CBCGA.Extended_lip

TT_ID_supple <- colnames(lipid_data_supplementary[, !str_detect(colnames(lipid_data_supplementary), fixed("_N"))])
lipid_data_supplementary_TT <- lipid_data_supplementary[,TT_ID_supple]


Basal_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="Basal",])
LumA_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="LumA",])
LumB_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="LumB",])
Her2_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="Her2",])
Normal_TT_ID <- row.names(CBCGA_Cohort_META[CBCGA_Cohort_META$PAM50_classifier=="Normal",])

Basal_TT_ID <- intersect(TT_ID_supple,Basal_TT_ID)
LumA_TT_ID <- intersect(TT_ID_supple,LumA_TT_ID)
LumB_TT_ID <- intersect(TT_ID_supple,LumB_TT_ID)
Her2_TT_ID  <- intersect(TT_ID_supple,Her2_TT_ID )
Normal_TT_ID  <- intersect(TT_ID_supple,Normal_TT_ID )


lipid_TT_data <- lipid_data_supplementary_TT
TT_ID_with_meta <- colnames(lipid_TT_data)
Basal_TT_ID_meta <- intersect(Basal_TT_ID,TT_ID_with_meta)#92
LumA_TT_ID_meta <- intersect(LumA_TT_ID,TT_ID_with_meta)#120
LumB_TT_ID_meta <- intersect(LumB_TT_ID,TT_ID_with_meta)#144
Her2_TT_ID_meta <- intersect(Her2_TT_ID,TT_ID_with_meta)#110
Normal_TT_ID_meta <- intersect(Normal_TT_ID,TT_ID_with_meta)#35
nonbasal_TT_ID_meta <- union(LumA_TT_ID_meta,LumB_TT_ID_meta)
nonbasal_TT_ID_meta <- union(nonbasal_TT_ID_meta,Her2_TT_ID_meta)
nonbasal_TT_ID_meta <- union(nonbasal_TT_ID_meta,Normal_TT_ID_meta)#409

####map/PAM50
#matrix basal vs non-basal
comparison_matrix5 <- matrix(ncol=5,nrow= nrow(lipid_TT_data))
colnames(comparison_matrix5) <- c("Basal_mean","LumA_mean","LumB_mean","Her2_mean","Normal_mean")
row.names(comparison_matrix5) <- row.names(lipid_TT_data)

for (i in rownames(comparison_matrix5)){
  comparison_matrix5[i,1] <- mean(as.numeric(lipid_TT_data [i,Basal_TT_ID_meta]))
  comparison_matrix5[i,2] <- mean(as.numeric(lipid_TT_data [i,LumA_TT_ID_meta]))
  comparison_matrix5[i,3] <- mean(as.numeric(lipid_TT_data [i,LumB_TT_ID_meta]))
  comparison_matrix5[i,4] <- mean(as.numeric(lipid_TT_data [i,Her2_TT_ID_meta]))
  comparison_matrix5[i,5] <- mean(as.numeric(lipid_TT_data [i,Normal_TT_ID_meta]))
}
comparison_matrix5 <- as.data.frame(comparison_matrix5)

id_oxpi <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPI"),])
id_oxpc <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPC"),])
id_oxpe <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPE"),])
id_oxpg <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPG"),])
id_oxps <- row.names(CBCGA_lip_anno[which(CBCGA_lip_anno$Subclass=="OxPS"),])
id_ox_lip <- union(id_oxpi,id_oxpc)
id_ox_lip <- union(id_ox_lip,id_oxpe)
id_ox_lip <- union(id_ox_lip,id_oxpg)
id_ox_lip <- union(id_ox_lip,id_oxps)

lipid_mapping_new_order <- CBCGA_lip_anno[id_ox_lip,]
lipid_id_new_order <- row.names(lipid_mapping_new_order)
lipid_comparison_new_order <- comparison_matrix5[lipid_id_new_order,]
###############
Plot_Data <- lipid_comparison_new_order
Plot_Data <- lipid_comparison_new_order[, c("LumA_mean","LumB_mean","Her2_mean","Basal_mean","Normal_mean")]
cluster_data <- Plot_Data
for (i in rownames(cluster_data)){
  cluster_data[i,] <- scale(as.numeric(cluster_data[i,]), center = T, scale = T)
}

########plot
library(pheatmap)
Plot_Data       <- cluster_data
Plot_Data       <- data.frame(Plot_Data,stringsAsFactors = F) 
datexpr2  <- as.data.frame(lapply(Plot_Data,as.numeric))
row.names(datexpr2) <- row.names(Plot_Data)
Plot_Data <- datexpr2

library(RColorBrewer)
pos      <- max(Plot_Data)+0.1 ; neg <- min(Plot_Data)-0.1
poscut   <- 100             ; negcut <- 100
mybreaks1<- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^ 2
mybreaks2<- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^ 2
mybreaks <- c(mybreaks1, mybreaks2)
mycolors <- c(colorRampPalette(c("#053061", "white"))(negcut), colorRampPalette(c("white", "#67001F"))(poscut))
#plot 

library(ggplot2)
library(pheatmap)

mycolors <- c(colorRampPalette(c("#0077C9", "white"))(negcut), colorRampPalette(c("white", "#E4002C"))(poscut))

pheatmap(Plot_Data, color = mycolors, breaks = mybreaks ,
         fontsize_row = 8,fontsize_col = 4, cluster_rows = F,cluster_cols=F,
         clustering_distance_cols = "euclidean", clustering_method = "ward.D", treeheight_col = 50,
         filename ="/results/Fig5E_right_ox_lipids.pdf"
         #         show_colnames = T, show_rownames = T, labels_col <- anno$metabolite_mapping_name
)





####46 class 

####input lipid_anno_new; reset the order of metabolites;saved in "ferroptosis_V3_20220720.Rdata"######

lipid_anno_new <- CBCGA_lip_anno
lipid_TT_test<- lipid_TT_data[row.names(CBCGA_lip_anno),]

lipid_TT_test <- cbind(lipid_TT_test,lipid_anno_new$Subclass)
colnames(lipid_TT_test)[604] <- 'Abbreviation'
test_value_1 <- tapply(lipid_TT_test[,1], lipid_TT_test$Abbreviation, mean)
test_value_2 <- tapply(lipid_TT_test[,2], lipid_TT_test$Abbreviation, mean)
#######matrix
inputdata_for_matrix <- lipid_TT_test
comparison_matrix <- matrix(ncol=604,nrow= 46)
colnames(comparison_matrix) <- colnames(inputdata_for_matrix)
#############
test_matrix <- aggregate(lipid_TT_test[,c("KPIY_pol_T","KFBP_pol_T")],by=list(lipid_TT_test$Abbreviation),mean)
###################
comparison_matrix[,1] <- test_value_1

test_matrix_2 <- as.data.frame(test_value_1)

for (i in 1:603){
  comparison_matrix[,i] <- tapply(lipid_TT_test[,i], lipid_TT_test$Abbreviation, mean)
}
comparison_matrix <- as.data.frame(comparison_matrix)
rownames(comparison_matrix) <- rownames(test_matrix_2)
#write.csv(comparison_matrix,"comparison_matrix_46class_mean_per_sample_20220806.csv")


transfer_46class_mean_per_sample <- t(comparison_matrix)
transfer_46class_mean_per_sample <- transfer_46class_mean_per_sample[-604,]


CBCGA_Cohort_META_ID <- row.names(CBCGA_Cohort_META)
CBCGA_Cohort_META_TT_ID <- intersect(TT_ID_supple,CBCGA_Cohort_META_ID)
CBCGA_Cohort_META_TT <- CBCGA_Cohort_META[CBCGA_Cohort_META_TT_ID,]
CBCGA_Cohort_META_TT_PAM50_ID <- intersect(CBCGA_Cohort_META_TT_ID,row.names(transfer_46class_mean_per_sample))
transfer_46class_mean_per_sample_PAM50 <- transfer_46class_mean_per_sample[CBCGA_Cohort_META_TT_PAM50_ID,]
CBCGA_Cohort_META_TT <- CBCGA_Cohort_META_TT[row.names(transfer_46class_mean_per_sample_PAM50),]

transfer_46class_mean_per_sample <- cbind(transfer_46class_mean_per_sample_PAM50,CBCGA_Cohort_META_TT$PAM50_classifier)
transfer_46class_mean_per_sample <- as.data.frame(transfer_46class_mean_per_sample) 
colnames(transfer_46class_mean_per_sample)[47] <- 'PAM50'
testdata <- transfer_46class_mean_per_sample
library(ggplot2)
testdata <- na.omit(testdata)
#construct matrix3

inputdata_for_matrix <- t(testdata[,1:46])

comparison_matrix <- matrix(ncol=5,nrow=nrow(inputdata_for_matrix))
rownames(comparison_matrix) <- rownames(inputdata_for_matrix)
LumA_TT_ID_meta

for (i in rownames(comparison_matrix)){
  a <- as.numeric(inputdata_for_matrix[i,LumA_TT_ID_meta])
  b <- as.numeric(inputdata_for_matrix[i,LumB_TT_ID_meta])
  c <- as.numeric(inputdata_for_matrix[i,Her2_TT_ID_meta])
  d <- as.numeric(inputdata_for_matrix[i,Basal_TT_ID_meta ])
  e <- as.numeric(inputdata_for_matrix[i,Normal_TT_ID_meta])
  comparison_matrix[i,1] <- mean(a)
  comparison_matrix[i,2] <- mean(b)
  comparison_matrix[i,3] <- mean(c)
  comparison_matrix[i,4] <- mean(d)
  comparison_matrix[i,5] <- mean(e)
}
comparison_matrix <- as.data.frame(comparison_matrix)
colnames(comparison_matrix) <- c("LumA_TT","LumB_TT","HER2E_TT","Basal_TT","Normal_TT")
comparison_matrix9 <- comparison_matrix
#write.csv(comparison_matrix,"comparison_matrix9_pam50_46class.csv")

Plot_Data <- comparison_matrix9
cluster_data <- Plot_Data


for (i in rownames(cluster_data)){
  cluster_data[i,] <- scale(as.numeric(cluster_data[i,]), center = T, scale = T)
}

########plot
library(pheatmap)
Plot_Data       <- cluster_data
Plot_Data       <- data.frame(Plot_Data,stringsAsFactors = F) 
datexpr2  <- as.data.frame(lapply(Plot_Data,as.numeric))
row.names(datexpr2) <- row.names(Plot_Data)
Plot_Data <- datexpr2
library(RColorBrewer)
pos      <- max(Plot_Data)+0.1 ; neg <- min(Plot_Data)-0.1
poscut   <- 100             ; negcut <- 100
mybreaks1<- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^ 2
mybreaks2<- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^ 2
mybreaks <- c(mybreaks1, mybreaks2)
library(ggplot2)
library(pheatmap)
mycolors <- c(colorRampPalette(c("#0077C9", "white"))(negcut), colorRampPalette(c("white", "#E4002C"))(poscut))
p3 <- pheatmap(Plot_Data, color = mycolors, breaks = mybreaks ,
               fontsize_row = 8,fontsize_col = 4, cluster_rows = T,cluster_cols=F,
               clustering_distance_cols = "euclidean", clustering_method = "ward.D", treeheight_col = 50,
               filename = "/results/Fig5E_left.pdf"
               #         show_colnames = T, show_rownames = T, labels_col <- anno$metabolite_mapping_name
)

#####Fig 5F volcano
####comparison_matrix10 volcano_46class
testdata <- transfer_46class_mean_per_sample
testdata <- na.omit(testdata)
inputdata_for_matrix <- t(testdata[,1:46])
comparison_matrix10 <- matrix(ncol=4,nrow=nrow(inputdata_for_matrix))
rownames(comparison_matrix10) <- rownames(inputdata_for_matrix)

for (i in rownames(comparison_matrix10)){
  a <- as.numeric(inputdata_for_matrix[i,Basal_TT_ID_meta])
  b <- as.numeric(inputdata_for_matrix[i,nonbasal_TT_ID_meta])
  comparison_matrix10[i,1] <- mean(a)
  comparison_matrix10[i,2] <- mean(b)
  c<-t.test(a,b,paired = FALSE)
  comparison_matrix10[i,3]<-c$p.value
}
comparison_matrix10[,4]<-p.adjust(c(comparison_matrix10[,3]),method = "fdr")
comparison_matrix10 <- as.data.frame(comparison_matrix10)
colnames(comparison_matrix10) <- c("basal_mean","nonbasal_mean","p","fdr")

comparison_matrix10$logFC <- comparison_matrix10[,1]-comparison_matrix10[,2]
## valcano plot

Expr             <- as.data.frame(comparison_matrix10)
Expr$label       <- rownames(Expr)
selected         <- matrix(nrow = 5,ncol = 1)
selected[,1] <- c("OxPC","OxPS","OxPI","OxPE","OxPG")
colnames(selected)[1] <- "gene"
row.names(selected) <- selected[,1]
Expr$gene        <- row.names(Expr)
selectgenes      <- merge(selected, Expr, by = "gene")
Expr$logFC <- Expr$basal_mean-Expr$nonbasal_mean
row.names(selectgenes) <- selectgenes$gene
selectgenes$logFC <- selectgenes$basal_mean-selectgenes$nonbasal_mean
library(ggplot2)
library(ggrepel)
#library(ggthemes)
#library(gridExtra)
logFCcut         <- 0.5
adjPCut          <- 0.05
xmin             <- (range(Expr$logFC)[1]- (range(Expr$logFC)[1]+ 1.5))
xmax             <- (range(Expr$logFC)[1]+ (1.5-range(Expr$logFC)[1]))
ymin             <- 0
ymax             <- -log10(Expr$p)[3] * 4
Expr$color       <- ifelse((Expr$p < adjPCut & Expr$logFC > logFCcut), "red", 
                           ifelse((Expr$p < adjPCut &Expr$logFC < -logFCcut), "blue","grey"))
size             <- ifelse((Expr$p < adjPCut & abs(Expr$logFC) > logFCcut), 2, 1)


Volcano          <- ggplot(data=Expr, aes( logFC, -log10(fdr), label = label)) +
  geom_point( size = 8, colour = Expr$color) +
  labs(x=bquote(~Log[2]~"(fold change)-basalvsnon"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(breaks = c(-1, -logFCcut, 0, logFCcut, 1),        
                     labels = c(-1, -logFCcut, 0, logFCcut, 1),
                     limits = c(-2, 2)) +
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", linetype="longdash", lwd = 0.6) +
  geom_hline(yintercept = -log10(adjPCut), color="grey40", linetype="longdash", lwd = 0.6) +
  theme_bw(base_size = 12) +
  theme(panel.grid=element_blank()) +
  geom_point(data = selectgenes, size= 8, shape = 1, color = "black") +
  guides(color=guide_legend(title = NULL))

ggsave(Volcano,file="/results/Fig5F_Volcano.pdf")

Expr_2 <- Expr
Expr_2[47,] <- NA
row.names(Expr_2)[47]<- "Ox"
a <- CBCGA.Extended_lip[id_ox_lip,Basal_TT_ID_meta]
a1 <- unlist(a)
a2 <- as.numeric(a1)

b <- CBCGA.Extended_lip[id_ox_lip,nonbasal_TT_ID_meta]
b1 <- unlist(b)
b2 <- as.numeric(b1)

Expr_2[47,1] <- mean(a2)
Expr_2[47,2] <- mean(b2)
c<-t.test(a2,b2,paired = FALSE)
Expr_2[47,3]<-c$p.value
Expr_2[47,5]<- Expr_2[47,1]-Expr_2[47,2]
#write.csv(Expr_2,"volcano_data_2022.csv")

###protein level Fig5G
####ferroptosis metabolism genes####
#####protein
protein_selected <- c("ACSL4","TFRC","LPCAT3","ALOX5")
protein_data_analysis <- CBCGA.Extended_PRO_normalized[protein_selected,]
dim(protein_data_analysis)
library(stringr)

CBCGA_Cohort_PRO <- myclinicaldata
row.names(CBCGA_Cohort_PRO) <- paste0(myclinicaldata$PatientCode,"_PRO_T","")

CBCGA_PRODATA <- CBCGA.Extended_PRO_normalized

library(stringr)
TT_ID_PRO <- colnames(CBCGA_PRODATA[, !str_detect(colnames(CBCGA_PRODATA), fixed("_N"))])
CBCGA_PRODATA_TT <- CBCGA_PRODATA[,TT_ID_PRO]



library(stringr)

Basal_TT_ID_pro <- row.names(CBCGA_Cohort_PRO[CBCGA_Cohort_PRO$PAM50_classifier=="Basal",])
LumA_TT_ID_pro <- row.names(CBCGA_Cohort_PRO[CBCGA_Cohort_PRO$PAM50_classifier=="LumA",])
LumB_TT_ID_pro <- row.names(CBCGA_Cohort_PRO[CBCGA_Cohort_PRO$PAM50_classifier=="LumB",])
Her2_TT_ID_pro <- row.names(CBCGA_Cohort_PRO[CBCGA_Cohort_PRO$PAM50_classifier=="Her2",])
Normal_TT_ID_pro <- row.names(CBCGA_Cohort_PRO[CBCGA_Cohort_PRO$PAM50_classifier=="Normal",])

Basal_TT_ID_pro <- intersect(TT_ID_PRO,Basal_TT_ID_pro)#59
LumA_TT_ID_pro <- intersect(TT_ID_PRO,LumA_TT_ID_pro)#56
LumB_TT_ID_pro <- intersect(TT_ID_PRO,LumB_TT_ID_pro)#77
Her2_TT_ID_pro  <- intersect(TT_ID_PRO,Her2_TT_ID_pro )#59
Normal_TT_ID_pro  <- intersect(TT_ID_PRO,Normal_TT_ID_pro)#20

LUM_ID_PRO <- union(LumA_TT_ID_pro,LumB_TT_ID_pro)
her2_normal_PRo <- union(Her2_TT_ID_pro,Normal_TT_ID_pro)
nonbasal_TT_ID_pro <- union(LUM_ID_PRO,her2_normal_PRo)
total_TT_PRO <- union(Basal_TT_ID_pro,nonbasal_TT_ID_pro)

protein_data_analysis_TT <- protein_data_analysis[,total_TT_PRO]
test4 <- protein_data_analysis[,Basal_TT_ID_pro]
##COMPARISON MATRIX 12 ---Fig 5G ferroptosis protein comparison
#matrix basal vs non-basal ACSL4_TFRC_ALOX5
tmp_data2 <- as.data.frame(protein_data_analysis_TT)
dim(tmp_data2)
tmp_data3 <- t(tmp_data2)
tmp_data3 <- na.omit(tmp_data3)
id_selected  <- row.names(tmp_data3)
tmp_data2 <- na.omit(tmp_data2)

comparison_matrix12 <- matrix(ncol=5,nrow= nrow(tmp_data2))
colnames(comparison_matrix12) <- c("Basal_mean","LumA_mean","LumB_mean","Her2_mean","Normal_mean")
row.names(comparison_matrix12) <- row.names(tmp_data2)

for (i in rownames(comparison_matrix12)){
  comparison_matrix12[i,1] <- mean(as.numeric(tmp_data2 [i,Basal_TT_ID_pro]))
  comparison_matrix12[i,2] <- mean(as.numeric(tmp_data2 [i,LumA_TT_ID_pro]))
  comparison_matrix12[i,3] <- mean(as.numeric(tmp_data2 [i,LumB_TT_ID_pro]))
  comparison_matrix12[i,4] <- mean(as.numeric(tmp_data2 [i,Her2_TT_ID_pro]))
  comparison_matrix12[i,5] <- mean(as.numeric(tmp_data2 [i,Normal_TT_ID_pro]))
}
comparison_matrix12 <- as.data.frame(comparison_matrix12)

#write.csv(comparison_matrix12,"comparison_matrix12_ACSL4_TFRC_ALOX5_PAM50.csv")

#LPCAT3(exclude NA)
tmp_data3 <- as.data.frame(t(tmp_data3))
comparison_matrix13 <- matrix(ncol=5,nrow= nrow(tmp_data3))
colnames(comparison_matrix13) <- c("Basal_mean","LumA_mean","LumB_mean","Her2_mean","Normal_mean")
row.names(comparison_matrix13) <- row.names(tmp_data3)
Basal_TT_ID_pro <- intersect(id_selected,Basal_TT_ID_pro)#48
LumA_TT_ID_pro <- intersect(id_selected,LumA_TT_ID_pro)#50

LumB_TT_ID_pro <- intersect(id_selected,LumB_TT_ID_pro)#69
Her2_TT_ID_pro <- intersect(id_selected,Her2_TT_ID_pro)#53
Normal_TT_ID_pro <- intersect(id_selected,Normal_TT_ID_pro)#17

for (i in rownames(comparison_matrix13)){
  comparison_matrix13[i,1] <- mean(as.numeric(tmp_data3 [i,Basal_TT_ID_pro]))
  comparison_matrix13[i,2] <- mean(as.numeric(tmp_data3 [i,LumA_TT_ID_pro]))
  comparison_matrix13[i,3] <- mean(as.numeric(tmp_data3 [i,LumB_TT_ID_pro]))
  comparison_matrix13[i,4] <- mean(as.numeric(tmp_data3 [i,Her2_TT_ID_pro]))
  comparison_matrix13[i,5] <- mean(as.numeric(tmp_data3 [i,Normal_TT_ID_pro]))
}
comparison_matrix13 <- as.data.frame(comparison_matrix13)
#write.csv(comparison_matrix13,"comparison_matrix13_LPCAT3_PAM50.csv")


