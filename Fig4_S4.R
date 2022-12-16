############################################################################################
#                                                                                          #
#                               CBCGA Fig4 & Fig S4                                  #
#                                                                                          #
############################################################################################


rm(list = ls())

load("/data/CBCGA.Extended_MergedData_V2.5_220722.Rdata")
load("/data/Fig4_S4_data.RData")

# part1：样本筛选及矩阵处理 ---------------------------------------------------------------

CBCGA_clinical    <-  CBCGA.Extended_Cohort.Info

patient <-  rownames(CBCGA_clinical[which(CBCGA_clinical$CBCGA_Cohort == "Yes"),])        # 全亚型顺反式

# patient <-  rownames(CBCGA_clinical[which(CBCGA_clinical$CBCGA_Cohort == "Yes" &        # 各PAM50亚型顺反式
#                                           CBCGA_clinical$PAM50_classifier == "LumA"),])

# CNV矩阵

CBCGA_cna          <- CBCGA_GISTICgene.alldata[,which(colnames(CBCGA_GISTICgene.alldata)%in%patient)] 

# RNA矩阵

RNA_matrix <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[,which(substring(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode),10,10) == "T")]

colnames(RNA_matrix) <- substring(colnames(RNA_matrix),1,4)

CBCGA_mRNA <- RNA_matrix[,which(colnames(RNA_matrix)%in%patient)] 

CBCGA_mRNA <- log2(CBCGA_mRNA+1)

# 蛋白矩阵

pro_matrix <- CBCGA.Extended_PRO_normalized[,which(substring(colnames(CBCGA.Extended_PRO_normalized),10,10) == "T")]

colnames(pro_matrix) <- substring(colnames(pro_matrix),1,4)

CBCGA_pro <- pro_matrix[,which(colnames(pro_matrix)%in%patient)] 

# 选择有三组学数据的样本及基因

common_sample <- Reduce(intersect,list(colnames(CBCGA_cna),colnames(CBCGA_mRNA),colnames(CBCGA_pro)))
common_gene <- Reduce(intersect,list(rownames(CBCGA_cna),rownames(CBCGA_mRNA),rownames(CBCGA_pro)))

CBCGA_cna <- CBCGA_cna[common_gene,common_sample]
CBCGA_mRNA <- CBCGA_mRNA[common_gene,common_sample]
CBCGA_pro <- CBCGA_pro[common_gene,common_sample]

table(colnames(CBCGA_cna)==colnames(CBCGA_mRNA))
table(colnames(CBCGA_cna)==colnames(CBCGA_pro))
table(colnames(CBCGA_mRNA)==colnames(CBCGA_pro))

# part2：按参考染色体顺序排列基因 ------------------------------------------------------------

gene_order <- gene_order_duplicated$gene_id
gene_order <- gene_order[gene_order%in%common_gene]

CBCGA_cna <- as.matrix(CBCGA_cna[gene_order,])
CBCGA_mRNA <- as.matrix(CBCGA_mRNA[gene_order,])
CBCGA_pro <- as.matrix(CBCGA_pro[gene_order,])

table(rownames(CBCGA_cna)==rownames(CBCGA_mRNA))
table(rownames(CBCGA_cna)==rownames(CBCGA_pro))
table(rownames(CBCGA_mRNA)==rownames(CBCGA_pro))  #7398基因

# part3：CNV-RNA顺反式计算 --------------------------------------------------------------------

cna_mRNA_Cor.Res <- matrix(ncol = length(gene_order),nrow = length(gene_order))
cna_mRNA_P.Vals <- matrix(ncol = length(gene_order),nrow = length(gene_order))

rownames(cna_mRNA_Cor.Res) <- gene_order
rownames(cna_mRNA_P.Vals) <- gene_order
colnames(cna_mRNA_Cor.Res) <- gene_order
colnames(cna_mRNA_P.Vals) <- gene_order

for(i in 1:length(gene_order)) {
  
  TEMP_CNA <- CBCGA_cna[i,]
  
  for(j in 1:length(gene_order)) {
    
    TEMP_Exp <- CBCGA_mRNA[j, ]
    
    cna_mRNA_Cor.Res[i,j] <- cor.test(TEMP_CNA, TEMP_Exp, method = "spearman")$estimate
    cna_mRNA_P.Vals[i,j]  <- cor.test(TEMP_CNA, TEMP_Exp, method = "spearman")$p.value
    
    if(j/500 == round(j/500)) {print(100*round(j/length(gene_order),2))}
  }
  print(paste("i = ", i, sep = ""))
}

# 转置后保存(此时行为CNA,列为EXP)

cna_mRNA_Cor.Res  <-  t(cna_mRNA_Cor.Res)
cna_mRNA_P.Vals   <-  t(cna_mRNA_P.Vals)

# 顺反式结果分别校正

p_matrix<-cna_mRNA_P.Vals
r_matirx<-cna_mRNA_Cor.Res

# 顺式校正
cis_p <- diag(p_matrix)
cis_fdr <- p.adjust(cis_p,method = "fdr")

# 反式按列校正

cna_mRNA_P.Vals_FDR <- p_matrix

for (i in 1:ncol(cna_mRNA_P.Vals_FDR)) {
  
  cna_mRNA_P.Vals_FDR[,i] <- p.adjust(cna_mRNA_P.Vals_FDR[,i],method = "fdr")
  
  print(i)
}
diag(cna_mRNA_P.Vals_FDR)  <- cis_fdr

# 染色体位置注释

anno_col <- data.frame(chro=as.factor(gene_order_duplicated[gene_order,"seqnames"]))
rownames(anno_col) <- gene_order

color_24=rep(rainbow(7)[-4],5)[1:24]
names(color_24)<-unique(anno_col$chro)
ann_colors = list(chro = color_24)

# 热图展示(Fig 4A)
library(pheatmap)

data <- r_matirx

for (i in 1:ncol(data)){
  
  for (j in 1:nrow(data)){
    
    if (cna_mRNA_P.Vals_FDR[j,i]>=0.05 & !is.na(cna_mRNA_P.Vals_FDR[j,i])){ data[j,i]<-NA }
  }
  print(i)
}

data[data==0] <- NA
data[data>0] <- 1
data[data<0] <- -1

data <- data[rev(1:nrow(data)),]

png('/results/Fig4A_left_All_CNV_mRNA_heatmap_v2.5.png',width = 1500, height =  1500,res = 800)

pheatmap(data, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("#3653a5", "white", "#e21a21"))(15),na_col = "White",
         annotation_col = anno_col,annotation_colors = ann_colors,annotation_legend = F,
         show_rownames = F,show_colnames = F,legend = F)

dev.off()

# part3：CNV-pro顺反式计算 --------------------------------------------------------------------

cna_pro_Cor.Res <- matrix(ncol = length(gene_order),nrow = length(gene_order))
cna_pro_P.Vals <- matrix(ncol = length(gene_order),nrow = length(gene_order))

rownames(cna_pro_Cor.Res) <- gene_order
rownames(cna_pro_P.Vals) <- gene_order
colnames(cna_pro_Cor.Res) <- gene_order
colnames(cna_pro_P.Vals) <- gene_order

for(i in 1:length(gene_order)) {
  
  TEMP_CNA <- CBCGA_cna[i,]
  
  for(j in 1:length(gene_order)) {
    
    TEMP_Exp <- CBCGA_pro[j,]
    
    cna_pro_Cor.Res[i,j] <- cor.test(TEMP_CNA, TEMP_Exp, method = "spearman")$estimate
    cna_pro_P.Vals[i,j]  <- cor.test(TEMP_CNA, TEMP_Exp, method = "spearman")$p.value
    
    if(j/500 == round(j/500)) {print(100*round(j/length(gene_order),2))}
  }
  print(paste("i = ", i, sep = ""))
}

# 转置后保存(此时行为CNA,列为EXP)

cna_pro_Cor.Res  <-  t(cna_pro_Cor.Res)
cna_pro_P.Vals   <-  t(cna_pro_P.Vals)

# 顺反式结果分别校正

p_matrix <- cna_pro_P.Vals
r_matirx <- cna_pro_Cor.Res

# 顺式校正
cis_p <- diag(p_matrix)
cis_fdr <- p.adjust(cis_p,method = "fdr")

# 反式按列校正

cna_pro_P.Vals_FDR <- p_matrix

for (i in 1:ncol(cna_pro_P.Vals_FDR)) {
  
  cna_pro_P.Vals_FDR[,i] <- p.adjust(cna_pro_P.Vals_FDR[,i],method = "fdr")
  
  print(i)
}
diag(cna_pro_P.Vals_FDR)  <- cis_fdr

# 染色体位置注释

anno_col <- data.frame(chro=as.factor(gene_order_duplicated[gene_order,"seqnames"]))
rownames(anno_col) <- gene_order

color_24=rep(rainbow(7)[-4],5)[1:24]
names(color_24)<-unique(anno_col$chro)
ann_colors = list(chro = color_24)

# 热图展示(Fig 4A)
library(pheatmap)

data <- cna_pro_Cor.Res

for (i in 1:ncol(data)){
  
  for (j in 1:nrow(data)){
    
    if (cna_pro_P.Vals_FDR[j,i]>=0.05 & !is.na(cna_pro_P.Vals_FDR[j,i])){ data[j,i]<-NA }
  }
  print(i)
}

data[data==0]   <- NA
data[data>0]  <- 1
data[data<0]  <- -1

data <- data[rev(1:nrow(data)),]

png('/results/Fig4A_right_CNV_pro_heatmap_v2.5.png',width = 1500, height =  1500,res = 800)

pheatmap(data, cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(c("#3653a5", "white", "#e21a21"))(15),na_col = "White",
         annotation_col = anno_col,annotation_colors = ann_colors,annotation_legend = F,
         show_rownames = F,show_colnames = F,legend = F)

dev.off()

# part4：全乳癌顺式直方图(Fig 4B upper) -----------------------------------------------

rm(list = ls())

#加载顺式相关性及p值矩阵

# 参数修改
FDR_filter_for_cis_p <- 0.05
file_names <- "CBCGA_All_"

# 顺式维度：CNV-mRNA及CNV-pro

omics <- "CNA_mRNA"
p_matrix <- cna_mRNA_P.Vals
r_matirx <- cna_mRNA_Cor.Res

# omics <- "CNA_pro"
# p_matrix <- cna_pro_P.Vals
# r_matirx <- cna_pro_Cor.Res

# 以下无需调整

cis_p <- diag(p_matrix)
cis_fdr <- p.adjust(cis_p,method = "fdr") #校正
cis_fdr <- -log10(cis_fdr)

cis_r <- diag(r_matirx)

gene_list <- colnames(p_matrix)
index <- 1:ncol(p_matrix)

chr_info <- gene_order_duplicated[gene_list,"seqnames"]

plot_dataframe<-data.frame(gene_id = gene_list, chr_c = chr_info, chr_n = chr_info, chr_info_num_total = chr_info,
                           log_cis_fdr = cis_fdr,cis_r = cis_r,gene_order = index)

# 只展示显著正相关

plot_dataframe[which(plot_dataframe$cis_r<0),"log_cis_fdr"] <- NA

# 染色体位置
gene_label <- unique(plot_dataframe$chr_n)
gene_label <-sort(as.numeric(gene_label),decreasing = F)


gene_label_order <- rep(0,length(gene_label))

for (i in gene_label){
  gene_label_order[i] <- round(median(plot_dataframe[plot_dataframe$chr_info_num_total %in% i,"gene_order"]))
}

gene_label[gene_label %in% 23] <- "X"

# 染色体末尾
chr_end_order <- plot_dataframe[!duplicated(plot_dataframe$chr_info_num_total),"gene_order"]
chr_end_order <- c(chr_end_order,max(index))

# 直方图
library(ggplot2)

pdf(paste("/results/Fig4B_upper",file_names,omics,"_Histogram_cis_FDR_v2.5_poscor",".pdf",sep=""),width = 9,height = 1.5)

ggplot(plot_dataframe, aes(x = gene_order, y = log_cis_fdr),fill=gene_order) +
  geom_bar(stat = "identity",  size = 0.1, colour = "#F17162") +
  scale_fill_manual(values = c("#F17162","#F17162"))+
  scale_x_continuous(breaks = gene_label_order,labels = gene_label)+
  scale_y_continuous(limits = c(0,30),breaks = c(0,10,20,30),labels = c(0,10,20,30))+
  geom_vline(xintercept = chr_end_order,lty=1,lwd=0.6,alpha=1,color="#9FA0A0")+
  geom_hline(yintercept =c(-log10(FDR_filter_for_cis_p)),lty=1,lwd=0.6,alpha=1,color="Red")+
  labs(x="Chromosome",y="-log10(FDR)")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

dev.off()

# part5：全乳癌及各亚型反式直方图 (Fig 4B lower/Fig 4G/Fig S4F/Table S3B)-------------------------------------------------------------------

rm(list=ls())

# 加载反式相关性及p值矩阵

# 参数修改
FDR_filter_for_trans_R_number <- 0.05
file_names <- "CBCGA_Normal_"

# 反式维度：CNV-mRNA/CNV-pro

omics <- "CNA_mRNA"
p_matrix <- cna_mRNA_P.Vals
r_matirx <- cna_mRNA_Cor.Res

# omics <- "CNA_pro"
# p_matrix <- cna_pro_P.Vals
# r_matirx <- cna_pro_Cor.Res

# 以下无需调整

gene_list <- colnames(p_matrix)
index <- 1:ncol(p_matrix)

# 顺式校正
cis_p <- diag(p_matrix)
cis_fdr <- p.adjust(cis_p,method = "fdr")

# 反式按列校正

FDR <- p_matrix

for (i in 1:ncol(FDR)) {
  
  FDR[,i] <- p.adjust(FDR[,i],method = "fdr")
  
  print(i)
}
diag(FDR) <- cis_fdr


# 统计反式数量及比率
trans_r_number <- rep(NA,ncol(p_matrix))

for (i in index){  
  
  trans_r_number[i] <- table(FDR[,i] < FDR_filter_for_trans_R_number)[2]  
  trans_r_number[i] <- trans_r_number[i]-as.numeric(FDR[i,i] < FDR_filter_for_trans_R_number)  #若对角线被纳入,从总数中减去
  print(i)
  
}

# Table S3B

chr_info <- gene_order_duplicated[gene_list,"seqnames"]
plot_dataframe <- data.frame(gene_id = gene_list,chr_c = chr_info, chr_n = chr_info, chr_info_num_total = chr_info,
                             trans_r_number = trans_r_number,gene_order = index)
write.csv(file=paste("/results/TableS3B_",file_names,omics,"_trans_number_table",".csv",sep=""),plot_dataframe,row.names = F)

# 绘图展示(Fig 4B/4G/S4F)

trans_r_number<-trans_r_number/max(index)
trans_r_number<-trans_r_number*100

plot_dataframe <- data.frame(gene_id = gene_list,chr_c = chr_info, chr_n = chr_info, chr_info_num_total = chr_info,
                             trans_r_number = trans_r_number,gene_order = index)

# 染色体位置

gene_label <- unique(plot_dataframe$chr_n)
gene_label <-sort(as.numeric(gene_label),decreasing = F)

gene_label_order <- rep(0,length(gene_label))

for (i in gene_label){
  gene_label_order[i] <- round(median(plot_dataframe[plot_dataframe$chr_info_num_total %in% i,"gene_order"]))
}

gene_label[gene_label %in% 23] <- "X"

# 染色体末尾
chr_end_order <- plot_dataframe[!duplicated(plot_dataframe$chr_info_num_total),"gene_order"]
chr_end_order <- c(chr_end_order,max(index))

# 直方图
library(ggplot2)
pdf(paste("/results/Fig4B_lower_",file_names,omics,"_Histogram_trans_freq_v2.5_poscor",".pdf",sep=""),width = 9,height = 1.5)
# pdf(paste("Fig4G_",file_names,omics,"_Histogram_trans_freq_v2.5_poscor",".pdf",sep=""),width = 9,height = 1.5)
# pdf(paste("FigS4F_",file_names,omics,"_Histogram_trans_freq_v2.5_poscor",".pdf",sep=""),width = 9,height = 1.5)

ggplot(plot_dataframe, aes(x = gene_order, y = trans_r_number),fill=gene_order) +
  geom_bar(stat = "identity",  size = 0.1, colour = "#317EB8") +
  scale_fill_manual(values = c("#317EB8","#317EB8"))+
  scale_x_continuous(breaks = gene_label_order,labels = gene_label)+
  scale_y_continuous(limits = c(0,40),breaks = c(0,10,20,30,40),labels = c(0,10,20,30,40))+   #Fig4B
  # scale_y_continuous(limits = c(0,10),breaks = c(0,5,10),labels = c(0,5,10))+               #Fig4G
  # scale_y_continuous(limits = c(0,20),breaks = c(0,10,20),labels = c(0,10,20))+             #FigS4F
  geom_vline(xintercept = chr_end_order,lty=1,lwd=0.6,alpha=1,color="#9FA0A0")+
  labs(x="Chromosome",y="Frequency(%)")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

dev.off()

# part6：统计各亚型CNA_RNA/CNA_pro顺式数量(Fig 4C) ------------------------------------------------
rm(list=ls())
# 加载各亚型顺反式结果

# 顺式校正
cis_p  <- diag(cna_mRNA_P.Vals)
cis_fdr  <- p.adjust(cis_p,method = "fdr")
cis_r <- diag(cna_mRNA_Cor.Res)
RNA_result <- cbind(cis_fdr,cis_r)

cis_sig_CNV_RNA <- rownames(RNA_result[RNA_result[,"cis_fdr"]<0.05&RNA_result[,"cis_r"]>0,])

cis_p  <- diag(cna_pro_P.Vals)
cis_fdr  <- p.adjust(cis_p,method = "fdr")
cis_r <- diag(cna_pro_Cor.Res)
pro_result <- cbind(cis_fdr,cis_r)

cis_sig_CNV_pro <- rownames(pro_result[pro_result[,"cis_fdr"]<0.05&pro_result[,"cis_r"]>0,])

# 韦恩图统计结果并计算比率
# library(grid)
# library(futile.logger)
# library(VennDiagram)
# venn.plot <- venn.diagram(
#   x = list('cis_RNA' = unique(cis_sig_CNV_RNA),'cis_pro' = unique(cis_sig_CNV_pro)),filename = NULL
# )
# 
# pdf(file="Fig4C_LumA_cis_gene_venn_v2.5_poscor.pdf")
# grid.draw(venn.plot)
# dev.off()

# part7：各亚型顺式结果计算及基因筛选(Table S3A) ----------------------------------------------

rm(list=ls())
# 加载大礼包及各亚型CNV-RNA/CNA-pro顺反式结果

# 参数修改
file_names <-"CBCGA_Normal"
subtype <- "Normal"

FDR_mRNA_filter <- 0.05     #fdr
FDR_Protein_filter <- 0.05

R_mRNA_filter <- 0.4        #rho
R_Protein_filter <- 0.4 

P_mRNA_filter <- 0.05       #p.value
P_Protein_filter <- 0.05

AMP_filter_alldata <- 0.5   #amp
AMP_filter_thre <- c(2,"2")
amp_percentage_filter <- 0.1

DEL_filter_alldata <- -0.5  #del
DEL_filter_thre <- c(-2,"-2")
del_percentage_filter <- 0.05 #0.01的太少了

# 以下无需调整

# 筛选各亚型样本(同part1)
CBCGA_clinical    <-  CBCGA.Extended_Cohort.Info

patient <-  rownames(CBCGA_clinical[which(CBCGA_clinical$CBCGA_Cohort == "Yes" &
                                            CBCGA_clinical$PAM50_classifier == subtype),])

CBCGA_cna <- CBCGA_GISTICgene.alldata[,which(colnames(CBCGA_GISTICgene.alldata)%in%patient)] 

RNA_matrix <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[,which(substring(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode),10,10) == "T")]
colnames(RNA_matrix) <- substring(colnames(RNA_matrix),1,4)
CBCGA_mRNA <- RNA_matrix[,which(colnames(RNA_matrix)%in%patient)] 

pro_matrix <- CBCGA.Extended_PRO_normalized[,which(substring(colnames(CBCGA.Extended_PRO_normalized),10,10) == "T")]
colnames(pro_matrix) <- substring(colnames(pro_matrix),1,4)
CBCGA_pro <- pro_matrix[,which(colnames(pro_matrix)%in%patient)] 

sample <- Reduce(intersect,list(colnames(CBCGA_cna),colnames(CBCGA_mRNA),colnames(CBCGA_pro)))
CNV <- CBCGA_GISTICgene.thre[colnames(cna_mRNA_P.Vals),sample]

# 构建结果表格

RP_matrix_all <- matrix(nrow = ncol(cna_mRNA_P.Vals),ncol =13)
rownames(RP_matrix_all) <- colnames(cna_mRNA_P.Vals)
colnames(RP_matrix_all) <- c("Chro_pos",
                             "mRNA_R","mRNA_P","mRNA_FDR",
                             "Protein_R","Protein_P","Protein_FDR",
                             "AMP_percent","AMP_check","DEL_percent","DEL_check",
                             "AMP_peak","DEL_peak")

RP_matrix_all <- as.data.frame(RP_matrix_all)

RP_matrix_all$mRNA_R <- diag(cna_mRNA_Cor.Res)
RP_matrix_all$mRNA_P <- diag(cna_mRNA_P.Vals)
RP_matrix_all$mRNA_FDR <- p.adjust(RP_matrix_all$mRNA_P,method = "fdr")

RP_matrix_all$Protein_R <- diag(cna_pro_Cor.Res)
RP_matrix_all$Protein_P <- diag(cna_pro_P.Vals)
RP_matrix_all$Protein_FDR <- p.adjust(RP_matrix_all$Protein_P,method = "fdr")


amp_fun <-function(x){
  y <-length(x)-table(as.numeric(x) %in% c("2",2))[1]
  return(y)
}

RP_matrix_all$AMP_percent <- round(apply(CNV,1,amp_fun)/ncol(CNV),3)

del_fun<-function(x){
  y<-length(x)-table(as.numeric(x) %in% c("-2",-2))[1]
  return(y)
}

RP_matrix_all$DEL_percent <- round(apply(CNV,1,del_fun)/ncol(CNV),3)

RP_matrix_all$Chro_pos <- CBCGA_GISTICgene_Cytoband[rownames(RP_matrix_all)]

RP_matrix_all[rownames(RP_matrix_all) %in% rownames(Genes_in_amp_peak),"AMP_peak"]<-Genes_in_amp_peak[rownames(RP_matrix_all[rownames(RP_matrix_all) %in% rownames(Genes_in_amp_peak),]),2] 
RP_matrix_all[rownames(RP_matrix_all) %in% rownames(Genes_in_del_peak),"DEL_peak"]<-Genes_in_del_peak[rownames(RP_matrix_all[rownames(RP_matrix_all) %in% rownames(Genes_in_del_peak),]),2] 

# 基因筛选
RP_matrix_all[,"AMP_check"] <- RP_matrix_all$AMP_percent>amp_percentage_filter
RP_matrix_all[,"DEL_check"] <- RP_matrix_all$DEL_percent>del_percentage_filter

# Table S3A
write.csv(RP_matrix_all, file = paste("/results/TableS3A_CBCGA_Cis_",file_names,"_all_data_RP_matrix_R_FDR","_AMP_percent_",amp_percentage_filter,"_DEL_percent_",del_percentage_filter,"_poscor",".csv",sep = ""))

check1 <- RP_matrix_all$mRNA_FDR < FDR_mRNA_filter
check2 <- RP_matrix_all$Protein_FDR < FDR_Protein_filter
check3 <- RP_matrix_all$mRNA_R>0       #正相关
check4 <- RP_matrix_all$Protein_R>0

check_FDR <- check1&check2&check3&check4&(RP_matrix_all[,"AMP_check"] | RP_matrix_all[,"DEL_check"])

RP_matrix_FDR <- RP_matrix_all[check_FDR,]

write.csv(RP_matrix_FDR, file = paste("/results/CBCGA_Cis_",file_names,"_all_filtered_RP_matrix_R_FDR","_AMP_percent_",amp_percentage_filter,"_DEL_percent_",del_percentage_filter,"_poscor",".csv",sep = ""))


# part8：各亚型peaks顺式基因热图展示(Fig 4D) ---------------------------------------

rm(list=ls())

# 读取顺式表格

read_index<-c(paste("CBCGA",c("LumA","LumB","Her2","Basal","Normal"),sep="_"))

cis <- list()

for (i in read_index) {
  temp <- read.csv(paste("/results/CBCGA_Cis_",i,"_all_filtered_RP_matrix_R_FDR_AMP_percent_0.1_DEL_percent_0.05.csv",sep=""),row.names=1)
  temp <- temp[!is.na(temp$AMP_peak),]
  cis[[i]] <- temp
}
common_peak <- unique(c(cis[[1]][,12],cis[[2]][,12],cis[[3]][,12],cis[[4]][,12],cis[[5]][,12]))


# peaks里的基因

genes_in_common_peak <- as.data.frame(Genes_in_amp_peak[Genes_in_amp_peak[,2]%in%common_peak,])

table(genes_in_common_peak$V2)

# 11p13 11q13.3 11q14.1 16p13.3   17q12 17q23.1  1q21.1  1q32.1 20q13.2 8p11.23  8q21.3 
# 24       1      16     211       2      19      60      70       9       3      14 

# 排除16p13.3(包含基因过多)
genes_in_common_peak <- genes_in_common_peak[which(genes_in_common_peak$V2!="16p13.3"),1]


# 判断上述基因测到多少
temp <- read.csv(paste("/results/TableS3A_CBCGA_Cis_",read_index[1],"_all_data_RP_matrix_R_FDR_AMP_percent_0.1_DEL_percent_0.05_poscor.csv",sep=""),row.names=1)

genes_in_common_peak <- rownames(temp)[rownames(temp)%in%genes_in_common_peak] #59个

# 提取FDR绘图表格

FDR_plot <- list()
FDR_save <- list()

read_index<-c(paste("CBCGA",c("LumA","LumB","Her2","Basal","Normal"),sep="_"))

for (i in read_index) {
  temp <- read.csv(paste("/results/TableS3A_CBCGA_Cis_",i,"_all_data_RP_matrix_R_FDR_AMP_percent_0.1_DEL_percent_0.05_poscor.csv",sep=""),row.names=1)
  
  FDR <- temp[genes_in_common_peak,c("mRNA_FDR","Protein_FDR")]
  FDR2 <- temp[genes_in_common_peak,c("mRNA_FDR","Protein_FDR","mRNA_R","Protein_R","AMP_percent")]
  colnames(FDR)<-paste(i,colnames(FDR),sep="_")
  colnames(FDR2)<-paste(i,colnames(FDR2),sep="_")
  
  FDR_plot[[i]] <- FDR
  FDR_save[[i]] <- FDR2
}

FDR_plot <-cbind(FDR_plot[[1]],FDR_plot[[2]],FDR_plot[[3]],FDR_plot[[4]],FDR_plot[[5]])
FDR_save <-cbind(FDR_save[[1]],FDR_save[[2]],FDR_save[[3]],FDR_save[[4]],FDR_save[[5]])

# 统计每一个基因显著的亚型个数
FDR_check <- function(x){
  y <- sum(x<0.05)
  y <- y>=1
  return(y)
}

FDR_plot_check <- apply(FDR_plot,1,FDR_check)

# 至少有1个亚型显著的基因
FDR_plot <- FDR_plot[FDR_plot_check,]
FDR_save <- FDR_save[FDR_plot_check,]

annotation_row <- as.data.frame(Genes_in_amp_peak[rownames(FDR_save),2])
colnames(annotation_row) <- "peaks"

cut_col <- c(rep("A",2),
             rep("B",2),
             rep("C",2))

# (inf为4 不scale 上限为5)

exprSet<-FDR_plot

exprSet_loged<-exprSet

exprSet_loged[exprSet_loged == 0.000000e+00] <- 10^(-4)

exprSet_loged<- -log10(exprSet_loged)

library( "pheatmap" )

# FDR热图
pdf('/results/Fig4D_CBCGA_AMP_peak_genes_sig_trans_heatmap_PAM50_v2.5_1_poscor.pdf',width = 3,height = 10)

choose_matrix = exprSet_loged

choose_matrix[choose_matrix > 4] = 4

#染色体anno
annotation_row2= data.frame(Chro_peaks =annotation_row[,1])
rownames(annotation_row2)<-rownames(annotation_row)

annotation_row2[,1]<-paste('Chr',annotation_row2[,1],sep = "")

choose_matrix=as.matrix(choose_matrix)

pheatmap::pheatmap( fontsize = 6, choose_matrix, annotation_row = annotation_row2,
                    cutree_cols = 4,
                    show_rownames = T, show_colnames = T, breaks = NA,
                    annotation_legend = T, cluster_cols = F,cluster_rows = F,
                    color = colorRampPalette(c( "white", "firebrick3"))(50),border=FALSE,
                    display_numbers = matrix(ifelse(choose_matrix > -log10(0.05), "*", ""), nrow(choose_matrix)),number_color="Black",fontsize_number=8,
                    legend_breaks = c(0:4), legend_labels = c("0","1","2","3",">=4"))

dev.off()

# 扩增频率热图

Amp <- FDR_save[,c(5,10,15,20,25)]

choose_matrix = Amp
choose_matrix[choose_matrix > 0.20] = 0.20
choose_matrix[choose_matrix < 0.10] = 0.10

annotation_row2= data.frame(Chro_peaks =annotation_row[,1])
rownames(annotation_row2)<-rownames(annotation_row)

annotation_row2[,1]<-paste('Chr',annotation_row2[,1],sep = "")

choose_matrix=as.matrix(choose_matrix)

pdf('/results/Fig4D_CBCGA_AMP_peak_genes_sig_trans_heatmap_PAM50_v2.5_amp_percentage.pdf',width = 3,height = 10)
pheatmap::pheatmap( fontsize = 6, choose_matrix, annotation_row = annotation_row2,
                    cutree_cols = 4,
                    show_rownames = T, show_colnames = T, breaks = NA,
                    annotation_legend = T, cluster_cols = F,cluster_rows = F,
                    color = colorRampPalette(c("#fff7ed","#fdaf61"))(200),border=FALSE,
                    legend_breaks = c(0.1:0.2), legend_labels = c("<=0.1"))

dev.off()

# part9：WWP1亚型对比(Fig 4E/4F/Fig S4E) ----------------------------------------------------------

# 加载大礼包结果

gene <- "WWP1"
# gene <- "CCND1"

CBCGA_clinical    <-  CBCGA.Extended_Cohort.Info

patient <-  rownames(CBCGA_clinical[which(CBCGA_clinical$CBCGA_Cohort == "Yes"),]) 

CNV <- as.data.frame(t(CBCGA_GISTICgene.thre[gene,which(colnames(CBCGA_GISTICgene.thre)%in%patient)]))
CNV$PatientCode <- rownames(CNV)

RNA <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[,which(substring(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode),10,10) == "T")]
RNA <- log2(RNA+1)
colnames(RNA) <- substring(colnames(RNA),1,4)
RNA <- as.data.frame(t(RNA[gene,which(colnames(RNA)%in%patient)]))
RNA$PatientCode <- rownames(RNA)

pro <- CBCGA.Extended_PRO_normalized[,which(substring(colnames(CBCGA.Extended_PRO_normalized),10,10) == "T")]
colnames(pro) <- substring(colnames(pro),1,4)
pro <- as.data.frame(t(pro[gene,which(colnames(pro)%in%patient)])) 
pro$PatientCode <- rownames(pro)

data <- merge(CNV,RNA,all.x = T,all.y = T,by = "PatientCode")
data <- merge(data,pro,all.x = T,all.y = T,by = "PatientCode")
data <- merge(data,CBCGA.Extended_Cohort.Info[,c(1,7)],by = "PatientCode")

data <- data[(data$WWP1.x==0|data$WWP1.x==1|data$WWP1.x==2)&is.na(data$WWP1.x)==F & is.na(data$PAM50_classifier)==F,]

data$WWP1.x <- as.factor(data$WWP1.x)

# 箱线图展示
library(ggplot2)
pdf(file = "/results/Fig4E_CBCGA_WWP1_GISTIC_RNA_boxplot_v2.5_log2.pdf",width = 6,height = 4)
# pdf(file = "FigS4E_CBCGA_CCND1_GISTIC_RNA_boxplot_v2.5_log2.pdf",width = 6,height = 4)

ggplot(data, aes(x = PAM50_classifier,y = WWP1.y,color = WWP1.x))+
  geom_boxplot(width=0.4,outlier.shape = NA)+
  scale_color_manual(values = c("#34495E","#E67F22","#D60000"))+
  scale_x_discrete(limits= c("LumA","LumB","Her2","Basal","Normal"))+
  # ylim(c(0,500))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  xlab(NULL)+ylab("WWP1 mRNA")

dev.off()

pdf(file = "/results/Fig4F_CBCGA_WWP1_GISTIC_pro_boxplot_v2.5.pdf",width = 6,height = 4)
# pdf(file = "FigS4E_CBCGA_CCND1_GISTIC_pro_boxplot_v2.5.pdf",width = 6,height = 4)

ggplot(data, aes(x = PAM50_classifier,y = WWP1,color = WWP1.x))+
  geom_boxplot(width=0.4,outlier.shape = NA)+
  scale_color_manual(values = c("#34495E","#E67F22","#D60000"))+
  scale_x_discrete(limits= c("LumA","LumB","Her2","Basal","Normal"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  xlab(NULL)+ylab("WWP1 Protein")

dev.off()

# 各亚型内检验
# CNV-RNA
kruskal.test(WWP1.y~WWP1.x,data[data$PAM50_classifier=="Normal",])

# CNV-pro
kruskal.test(WWP1~WWP1.x,data[data$PAM50_classifier=="Normal",])

# part10：LumB 8q反式调控基因富集分析(Fig 4H/Fig S4E) --------------------------------------------------------------------

# 加载LumB RNA及蛋白顺反式结果

# 顺反式分别校正

p_matrix <- cna_mRNA_P.Vals
r_matrix <- cna_mRNA_Cor.Res

# p_matrix <- cna_pro_P.Vals
# r_matrix <- cna_pro_Cor.Res

# 顺式校正
cis_p <- diag(p_matrix)
cis_fdr <- p.adjust(cis_p,method = "fdr")

# 反式按列校正

FDR <- p_matrix

for (i in 1:ncol(FDR)) {
  
  FDR[,i] <- p.adjust(FDR[,i],method = "fdr")
  
  print(i)
}
diag(FDR) <- cis_fdr

FDR <- ifelse(FDR>=0.05,0,1)  #显著反式为1,其余为0
diag(FDR) <- 0 #对角线为顺式

table(FDR)
# RNA
#     0         1 
# 54454198   276206 
# pro
#     0         1 
# 54644666    85738

# 挑选8q反式调控基因

gene_8q <- CBCGA_GISTICgene_Cytoband[substring(CBCGA_GISTICgene_Cytoband,1,2)=="8q"]
gene_8q <- names(gene_8q)


FDR <- FDR[,which(colnames(FDR)%in%gene_8q)]
trans_gene_8q <- rownames(FDR[which(rowSums(FDR)>0),])

trans_gene_8q <- unique(trans_gene_8q)
trans_gene_8q <- trans_gene_8q[!trans_gene_8q%in%colnames(FDR)] #排除顺式基因

# 富集分析
library(org.Hs.eg.db)
library(clusterProfiler)
# library(dplyr)
# library(ggplot2)
library(tidyverse)

gene <- mapIds(org.Hs.eg.db,trans_gene_8q,column = 'ENTREZID', keytype = 'SYMBOL',multiVals = 'filter')

GO <- enrichGO(gene = na.omit(gene),OrgDb=org.Hs.eg.db,ont="ALL",pAdjustMethod = "fdr",qvalueCutoff  = 0.05)   
GO <- as.data.frame(GO)

GO_plot <- rbind(GO[GO$ONTOLOGY=="BP",][1:5,],GO[GO$ONTOLOGY=="CC",][1:5,],GO[GO$ONTOLOGY=="MF",][1:5,])

# 画图
pdf(file = "/results/Fig4H_CBCGA_LumB_8q_trans_pro_enrichment_v2.5.pdf")
# pdf(file = "FigS4G_CBCGA_LumB_8q_trans_RNA_enrichment_v2.5.pdf")

ggplot(data= GO_plot, aes(x = number,y = Count,fill = ONTOLOGY))+
  geom_bar(stat="identity", width=0.8) + coord_flip()+ 
  scale_fill_manual(values = c("#5C9F52","#DA4F08","#4E84AA")) + theme_test()+ 
  scale_x_discrete(labels=rev(GO_plot$Description))+xlab(NULL)                         
theme(axis.text=element_text(size = 8,face = "bold", color="black"))

dev.off()

