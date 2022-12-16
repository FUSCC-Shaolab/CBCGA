rm(list = ls())
load("/data/FigS4A-C_FigS1C_data.Rdata")
###################################
## Fig S4A: NAs in each Batch
###################################
na_batch_index<-unique(CBCGA_PRO_batch_anno$BatchId)
x<-c()
z<-c()
for (i in na_batch_index){
  na_batch_matrix<-CBCGA_PRO_ori_matrix[,rownames(CBCGA_PRO_batch_anno[CBCGA_PRO_batch_anno$BatchId %in% i,])[1]]
  y<-as.data.frame(table(is.na(na_batch_matrix)))[1,2]
  x<-c(x,y)
  dim_matrix<-length(na_batch_matrix)
  z<-c(z,y/dim_matrix)
}
names(x)<-unique(CBCGA_PRO_batch_anno$BatchId)
names(z)<-unique(CBCGA_PRO_batch_anno$BatchId)

pdf("/results/FigS4A.pdf",heigh=5,width = 15)
barplot(x)
title("NA in Batchs")
dev.off()

#######################################
##median intensity ratio filtering
#######################################
QC_matrix <- as.matrix(CBCGA_PRO_ori_matrix)

#calculating median intensity ratio in each sample
QC_matrix <- as.matrix(t(QC_matrix))
QC_data <- data.frame("PatientCode_long"=rownames(QC_matrix),"median"="")
rownames(QC_data) <- QC_data$PatientCode_long

QC_data$median <- apply(QC_matrix,1,median,na.rm = T)

##filtering out outliers for tumor samples

QC_data_tumor <- QC_data[substr(QC_data$PatientCode_long,10,10) %in% c("T"),]
temp_check <-  (QC_data_tumor$median < mean(QC_data_tumor$median)+2*IQR(QC_data_tumor$median)) &  (QC_data_tumor$median > mean(QC_data_tumor$median)-2*IQR(QC_data_tumor$median))

QC_data_tumor_ourliers_filtered <- QC_data_tumor[temp_check,]

QC_data_PT <- QC_data[substr(QC_data$PatientCode_long,10,10) %in% c("N"),]
temp_check <-  (QC_data_PT$median < mean(QC_data_PT$median)+2*IQR(QC_data_PT$median)) &  (QC_data_PT$median > mean(QC_data_PT$median)-2*IQR(QC_data_PT$median))
table(temp_check)

QC_data_PT_ourliers_filtered <- QC_data_PT[temp_check,]


QC_matrix_median_filtered <- QC_matrix[c(QC_data_tumor_ourliers_filtered$PatientCode_long,QC_data_PT_ourliers_filtered$PatientCode_long),]

QC_matrix_median_filtered <- t(QC_matrix_median_filtered)
CBCGA_PRO_batch_anno_median_filtered <- CBCGA_PRO_batch_anno[colnames(QC_matrix_median_filtered),]


#######################################
##normalization
#######################################
QC_matrix_median_filtered_log2 <- log2(QC_matrix_median_filtered)
QC_matrix_median_filtered_log2_median_centered <- scale(QC_matrix_median_filtered_log2) ##scale by columns

QC_matrix_median_filtered_log2_median_centered_rebatch <- limma::removeBatchEffect(QC_matrix_median_filtered_log2_median_centered,CBCGA_PRO_batch_anno_median_filtered[colnames(QC_matrix_median_filtered_log2_median_centered),]$BatchId)

select_matrix <- QC_matrix_median_filtered_log2_median_centered_rebatch

x<-c()
for (i in 1:nrow(select_matrix)){
  na_batch_matrix<-select_matrix[i,]
  y=0
  for( j in 1:ncol(select_matrix)){
    if (is.na(na_batch_matrix[j])){
      y<-y+1
    }
    z<-y/ncol(select_matrix)
  }
  x<-c(x,z)
}
barplot(x)
na_protein_per<-x
title("NA% in Proteins")

na_check<-na_protein_per<0.3

na_30<-rownames(select_matrix)[na_check]

QC_matrix_median_filtered_log2_median_centered_rebatch_na_30<-QC_matrix_median_filtered_log2_median_centered_rebatch[na_30,]

#######################################
## Fig S4B: batch effect evaluation
#######################################
exp_temp<-QC_matrix_median_filtered_log2_median_centered_rebatch_na_30
anno<-CBCGA_PRO_batch_anno_median_filtered[colnames(exp_temp),]

{library(ggfortify)
  
  dat<-cbind(paste("Batch",stringr::str_pad(anno$BatchId,2,side="left","0")),apply(t(exp_temp), 2, as.numeric))
  colnames(dat)[1]<-c("label")
  dat<-t(na.omit(t(dat)))
  
  out_pca <- prcomp(apply((dat[,-1]), 2, as.numeric))
  dat[,1]<-as.character(dat[,1])
  
  p<-autoplot(out_pca,data=dat,colour='label',size=3,label=F,label.size=3,palette="aaas")
  
}
pdf("/results/FigS4B.pdf",height = 4.3,width = 6)
print(p)
dev.off()


#############################
##PCA filtering
#############################
##set pca_confidence_level
pca_confidence_level<-0.9

test_matrix <- QC_matrix_median_filtered_log2_median_centered_rebatch_na_30

test_anno <- CBCGA_PRO_batch_anno_median_filtered


library(ggfortify)
dat<-cbind(test_anno$DiseaseType,t(test_matrix))
colnames(dat)[1]<-c("label")
dat[,1]<-as.character(dat[,1])
dat<-t(na.omit(t(dat)))

out_pca <- prcomp(apply(dat[,-1],2,as.numeric))

# library(ggplot2)
# library(ggord)
# 
# PCA1<-out_pca[["x"]]
# rownames(PCA1)<-colnames(test_matrix)
# rownames(PCA1)<-rownames(test_anno)
# 
# cancer_samples_test_matrix<-test_anno[test_anno$DiseaseType %in% c("C"),"PatientCode_long"]
# 
# 
# PT_samples_test_matrix<-test_anno[test_anno$DiseaseType %in% c("P"),"PatientCode_long"]
# 
# ##PT filtering
# PCA1_range_PT<-car::dataEllipse(PCA1[PT_samples_test_matrix,"PC1"],
#                                 PCA1[PT_samples_test_matrix,"PC2"],
#                                 levels=pca_confidence_level)
# 
# within_polygon_check_PT<-sp::point.in.polygon(point.x=PCA1[PT_samples_test_matrix,"PC1"],
#                                               point.y=PCA1[PT_samples_test_matrix,"PC2"], 
#                                               pol.x=PCA1_range_PT[,1], 
#                                               pol.y=PCA1_range_PT[,2])
# 
# within_polygon_check_PT<-(within_polygon_check_PT %in% c(1,2,3))
# 
# PT_samples_test_matrix<-PT_samples_test_matrix[within_polygon_check_PT]
# 
# 
# ##T filtering
# PCA1_range_T<-car::dataEllipse(PCA1[cancer_samples_test_matrix,"PC1"],
#                                PCA1[cancer_samples_test_matrix,"PC2"],
#                                levels=pca_confidence_level,)
# 
# within_polygon_check_T<-sp::point.in.polygon(point.x=PCA1[cancer_samples_test_matrix,"PC1"],
#                                              point.y=PCA1[cancer_samples_test_matrix,"PC2"], 
#                                              pol.x=PCA1_range_T[,1], 
#                                              pol.y=PCA1_range_T[,2])
# 
# within_polygon_check_T<-(within_polygon_check_T %in% c(1,2,3))
# 
# 
# cancer_samples_test_matrix<-cancer_samples_test_matrix[(within_polygon_check_T)]
# 
# sample_by_PCA <- na.omit(c(cancer_samples_test_matrix,PT_samples_test_matrix))
# 
# test_matrix_PCA_0.9 <- test_matrix[,sample_by_PCA]
# 
# #################
# ### Fig S4C: technical replicates correlation
# #################
# PCA_0.9_sample_names <- colnames(test_matrix_PCA_0.9)
# 
# rep_tec_samples <- PCA_0.9_sample_names[stringr::str_detect(PCA_0.9_sample_names,"REP_TEC")]
# rep_tec_samples_all <- rep_tec_samples
# 
# rep_tec_samples <- unique(substr(rep_tec_samples,1,10))
# 
# library(stringr)
# 
# exp <- test_matrix_PCA_0.9
# library(ggplot2)
# 
# pdf("/results/FigS4C.pdf")
# for (i in rep_tec_samples){
#   temp_sample <- PCA_0.9_sample_names[stringr::str_detect(PCA_0.9_sample_names,substr(i,1,10))]
#   temp_sample <- rep_tec_samples_all[stringr::str_detect(rep_tec_samples_all,substr(i,1,10))]
#   temp<-exp[,temp_sample]
#   colnames(temp) <- gsub("F","Batch",CBCGA_PRO_batch_anno[colnames(temp),"label"])
#   kk <- GGally::ggpairs(temp,title=paste("Technical replicate: ",substr(i,1,4),sep = ""),upper = list(continuous = "cor"))
#   print(kk)
# }
# dev.off()
# 
# #################
# ###merge replicates 
# #################
# test_matrix_PCA_0.9_norep <- test_matrix_PCA_0.9
# 
# for (i in rep_tec_samples){
#   temp_sample <- PCA_0.9_sample_names[stringr::str_detect(PCA_0.9_sample_names,substr(i,1,10))]
#   temp<-exp[,temp_sample]
#   temp<-data.frame(temp)
#   temp <- rowMeans(temp,na.rm = T)
#   
#   if (! substr(i,1,10) %in% colnames(test_matrix_PCA_0.9_norep)){
#     test_matrix_PCA_0.9_norep <- data.frame(test_matrix_PCA_0.9_norep,temp)
#     colnames(test_matrix_PCA_0.9_norep)[ncol(test_matrix_PCA_0.9_norep)] <- substr(i,1,10)
#   }else{test_matrix_PCA_0.9_norep[,substr(i,1,10)] <- as.numeric(temp)}
#   
#   print(i)
# }
# 
# 
# rep_bio_samples <- PCA_0.9_sample_names[stringr::str_detect(PCA_0.9_sample_names,"REP_BIO")]
# for (i in unique(substr(rep_bio_samples,1,10))){
#   temp_sample <- PCA_0.9_sample_names[stringr::str_detect(PCA_0.9_sample_names,substr(i,1,10))]
#   temp<-exp[,temp_sample]
#   temp<-data.frame(temp)
#   temp <- rowMeans(temp,na.rm = T)
#   
#   if (! substr(i,1,10) %in% colnames(test_matrix_PCA_0.9_norep)){
#     test_matrix_PCA_0.9_norep <- data.frame(test_matrix_PCA_0.9_norep,temp)
#     colnames(test_matrix_PCA_0.9_norep)[ncol(test_matrix_PCA_0.9_norep)] <- substr(i,1,10)
#   }else{test_matrix_PCA_0.9_norep[,substr(i,1,10)] <- as.numeric(temp)}
#   
#   print(i)
# }
# 
# 
# test_matrix_PCA_0.9_norep <- test_matrix_PCA_0.9_norep[,nchar(colnames(test_matrix_PCA_0.9_norep)) ==10]
# 
# #################
# ###gene annotation and merge genes 
# #################
# CBCGA_PRO_gene_annotation <- na.omit(CBCGA_PRO_gene_annotation)
# 
# test_matrix_PCA_0.9_norep_symbol <- test_matrix_PCA_0.9_norep[rownames(test_matrix_PCA_0.9_norep) %in% rownames(CBCGA_PRO_gene_annotation),]
# 
# CBCGA_PRO_gene_annotation_sub <- CBCGA_PRO_gene_annotation[rownames(test_matrix_PCA_0.9_norep_symbol),]
# 
# CBCGA_PRO_gene_annotation_sub_dup <- CBCGA_PRO_gene_annotation_sub$Gene[duplicated(CBCGA_PRO_gene_annotation_sub$Gene)]
# 
# test_matrix_PCA_0.9_norep_symbol2 <- test_matrix_PCA_0.9_norep_symbol
# 
# exp2 <- test_matrix_PCA_0.9_norep_symbol
# 
# for (i in CBCGA_PRO_gene_annotation_sub_dup){
#   temp_gene <- CBCGA_PRO_gene_annotation_sub[CBCGA_PRO_gene_annotation_sub$Gene %in% i,"Label"]
#   
#   temp<-exp2[temp_gene,]
#   temp<-data.frame(temp)
#   temp <- colMeans(temp,na.rm = T)
#   
#   test_matrix_PCA_0.9_norep_symbol2[temp_gene[1],] <- as.numeric(temp)
#   test_matrix_PCA_0.9_norep_symbol2 <- test_matrix_PCA_0.9_norep_symbol2[! rownames(test_matrix_PCA_0.9_norep_symbol2) %in% temp_gene[-1],]
#   
#   print(i)
# }
# 
# rownames(test_matrix_PCA_0.9_norep_symbol2) <- CBCGA_PRO_gene_annotation_sub[rownames(test_matrix_PCA_0.9_norep_symbol2) ,"Gene"]
# 
# #write.csv(file = )

#############################
## Fig S1C:  heatmap for CBCGA
#############################

{library(genefu)
  library(org.Hs.eg.db)
  data("AIMSmodel")
  
  AIMS_GENE<-AIMSmodel$all.pairs
  AIMS_GENE<-strsplit(AIMS_GENE,"<")
  AIMS_GENE<-data.frame(AIMS_GENE)
  AIMS_GENE<-apply(AIMS_GENE, 2, as.numeric)
  AIMS_GENE<-as.vector(AIMS_GENE)
  
  AIMS_GENE<-clusterProfiler::bitr(AIMS_GENE,fromType = "ENTREZID",toType = "SYMBOL",org.Hs.eg.db)
  AIMS_GENE<-as.vector(AIMS_GENE[,2])}


exprSet<-test_matrix_PCA_0.9_norep_symbol2
exprSet <- exprSet[,substr(colnames(exprSet),10,10) %in% "T"]
exprSet <- exprSet[,CBCGA_PRO_batch_anno[substr(colnames(exprSet),1,10),"CBCGA_cohort"] %in% "Yes"]
choose_gene<-AIMS_GENE[AIMS_GENE %in% rownames(exprSet)]


Color_PAM50 = c('LumA' = '#0077c9', 'LumB' = '#74d2e8','HER2-enriched' = '#7552cd', 'Basal-like' = '#e4002c',
                'Normal-like' = '#cecece')

Anno_Color <- list( Intrinsic_subtype_PAM50 = Color_PAM50,
                    Intrinsic_subtype_AIMS = Color_PAM50)


methods <- c("ward.D2")
choose_matrix = exprSet[choose_gene , ]
choose_matrix <- na.omit(choose_matrix)
choose_matrix = t( scale( t( choose_matrix ) ) )

annotation_col = data.frame(Intrinsic_subtype_PAM50 = CBCGA_PRO_batch_anno[substr(colnames(choose_matrix),1,10),"Intrinsic_subtype_PAM50"],
                            Intrinsic_subtype_AIMS = CBCGA_PRO_batch_anno[substr(colnames(choose_matrix),1,10),"Intrinsic_subtype_AIMS"])
annotation_col$Intrinsic_subtype_PAM50[annotation_col$Intrinsic_subtype_PAM50 == "Her2"] = 'HER2-enriched'
annotation_col$Intrinsic_subtype_PAM50[annotation_col$Intrinsic_subtype_PAM50 == "Basal"] = 'Basal-like'
annotation_col$Intrinsic_subtype_PAM50[annotation_col$Intrinsic_subtype_PAM50 == "Normal"] = 'Normal-like'
annotation_col$Intrinsic_subtype_AIMS[annotation_col$Intrinsic_subtype_AIMS == "Her2"] = 'HER2-enriched'
annotation_col$Intrinsic_subtype_AIMS[annotation_col$Intrinsic_subtype_AIMS == "Basal"] = 'Basal-like'
annotation_col$Intrinsic_subtype_AIMS[annotation_col$Intrinsic_subtype_AIMS == "Normal"] = 'Normal-like'
rownames(annotation_col) <- colnames(choose_matrix)

pos <- 3 + 0.1
neg <- -3.1 + 0.1
poscut <- 50
negcut <- 50
mybreaks1 <- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^2
mybreaks2 <- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^2
mybreaks  <- c(mybreaks1, mybreaks2)
col =  c(colorRampPalette(c("#3A53A4", "#FFFFFF"))(51), colorRampPalette(c( "#FFFFFF", "#E21E26"))(50))


library(pheatmap)
pheatmap(choose_matrix, fontsize = 6,
         color =col,
         annotation_col = annotation_col,
         cutree_cols = 4,
         show_rownames = T, show_colnames = F, breaks = mybreaks,
         annotation_colors = Anno_Color,
         clustering_method = methods,
         annotation_legend = T, cluster_cols = T, main=methods,
         filename = "/results/FigS1C.pdf")




