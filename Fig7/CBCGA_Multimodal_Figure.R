rm(list=ls())

library(stringr)
library(dplyr)
library(glmnet)
library(pROC)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(PMCMRplus)
library(pheatmap)
library(ggthemes)
library(export)
library(survcomp)
library(Hmisc)

setwd("CBCGA_Multimodal_Code_20221216")


Modal_Res_Ord <- read.csv("Result/CBCGA_Multimodal_Cindex.csv")
colnames(Modal_Res_Ord)[1] <- "Abbr"
Modal_Res_Ord$Abbr <- factor(Modal_Res_Ord$Abbr, levels = c(Modal_Res_Ord$Abbr))



####### C-index Barplot #######


col <- c(rep(c("#6281a5","#594289","#ca9d9d","#b33e38","#dccf24","#ddb529","#478242","#beced2"),4),"#6281a5","#594289","#ca9d9d","#b33e38","#dccf24")
col <- c(rep(c("#beced2"),3),rep(c("#6281a5"),11),rep(c("#ca9d9d"),19),rep(c("#b33e38"),4))

ggplot(Modal_Res_Ord,aes(x=Abbr,y=C_index_test)) +
  geom_bar(stat='identity',fill = col,width = 0.7) +
  theme_base() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+ 
  coord_cartesian(ylim=c(0.45,0.8))



######## Model ########

source("CBCGA_Multimodal.R") # output all multimodality results

source("MultiModal/307 Met Path Rad Clin.R") # We source the specific model we select
# C_index 0.7407834



######## Set Color ########

mypal <- pal_startrek("uniform",alpha=1)(2)


anno_color <- list(

  T_stage = c("pT1"=brewer.pal(3,"Blues")[1], "pT2"=brewer.pal(3,"Blues")[2], "pT3"=brewer.pal(3,"Blues")[3]),
  N_stage = c("pN0"=brewer.pal(4,"Oranges")[1], "pN1"=brewer.pal(4,"Oranges")[2], "pN2"=brewer.pal(4,"Oranges")[3], "pN3"=brewer.pal(4,"Oranges")[4]),
  IHC_subtype = c("HR+HER2-" ="#0085c4", "HR+HER2+" ="#7ab801", "HR-HER2+" ="#f2af01","TNBC"="#dc5035"),
  PAM50_subtype = c("LumA"="#0077c9", "LumB"="#74d2e8", "Her2"="#7552cd", "Basal"="#e4002c", "Normal"="#cecece","Unknown"="#f2f2f2"),
  riskscore = c("#edf8fb", "#2ca25f"),
  risk = c("low" = "#6281a5", "high" = "#b33e38"),
  dataset = c("Train"="#beced2", "Test"="#6281a5"),
  RFS_status = c("0" ="#cecece", "1" ="black"),
  Feat_Cate = c("Meatb" = "#6281a5", "Path" = "#b33e38", "Rad" = "#ddb529")

)



####### Survival Analysis ########


riskscore_train_MPRC <- predict(res.cox, newdata = Trainset_Met_Path_Rad_Clin, type = "risk")
risk_train_MPRC <- cbind(riskscore_train_MPRC,Trainset_Met_Path_Rad_Clin)
colnames(risk_train_MPRC)[1] <- "riskscore"

riskscore_test_MPRC <- predict(res.cox, newdata = Testset_Met_Path_Rad_Clin, type = "risk")
risk_test_MPRC <- cbind(riskscore_test_MPRC,Testset_Met_Path_Rad_Clin)
colnames(risk_test_MPRC)[1] <- "riskscore"

hist(sort(risk_train_MPRC$riskscore))

risk_train_MPRC$riskscore <- log2(risk_train_MPRC$riskscore)
risk_test_MPRC$riskscore <- log2(risk_test_MPRC$riskscore)


train_tmp <- risk_train_MPRC
test_tmp <- risk_test_MPRC


p_test <- c()

for (i in seq(-5,5,0.1)) {
  
  train_tmp$risk <- ifelse(train_tmp$riskscore>=i,"high","low")

  
  surv_diff_train_MPRC_tmp <- survdiff(Surv(RFS_time, RFS_status) ~ risk, data = train_tmp)
  
  p1 = 1 - pchisq(surv_diff_train_MPRC_tmp$chisq, length(surv_diff_train_MPRC_tmp$n) - 1) 
  
  train_highrisk = sum(train_tmp$risk=="high")
  train_lowrisk = sum(train_tmp$risk=="low")
  
  
  
  test_tmp$risk <- ifelse(test_tmp$riskscore>=i,"high","low")
  
  surv_diff_test_MPRC_tmp <- survdiff(Surv(RFS_time, RFS_status) ~ risk, data = test_tmp)
  
  p2 = 1 - pchisq(surv_diff_test_MPRC_tmp$chisq, length(surv_diff_test_MPRC_tmp$n) - 1) 
  
  test_highrisk = sum(test_tmp$risk=="high")
  test_lowrisk = sum(test_tmp$risk=="low")
  
  p_test <- rbind(p_test,c(i,train_highrisk,train_lowrisk,p1,test_highrisk,test_lowrisk,p2))
  
  
}

p_test <- as.data.frame(p_test)
colnames(p_test) <- c("cutoff","train_highrisk","train_lowrisk","p_train","test_highrisk","train_lowrisk","p_test")
median(train_tmp$riskscore) 


###### Risk Score in Multi Trainset

Cutoff = 5

risk_train_MPRC$risk <- ifelse(risk_train_MPRC$riskscore>=Cutoff,"high","low")


####### KM-plot in Multi Trainset

# Model Fit
surv_model_train_MPRC <- survfit(Surv(RFS_time, RFS_status) ~ risk, data = risk_train_MPRC)
# Log-rank test
surv_diff_train_MPRC <- survdiff(Surv(RFS_time, RFS_status) ~ risk, data = risk_train_MPRC)

1 - pchisq(surv_diff_train_MPRC$chisq, length(surv_diff_train_MPRC$n) - 1) 

# KMplot
survPlot_MPRC_Train <- ggsurvplot(surv_model_train_MPRC, 
                                 data = risk_train_MPRC,
                                 conf.int = FALSE,  
                                 palette = mypal, 
                                 pval = TRUE,
                                 pval.coord = c(0, 0.6),
                                 pval.size = 5,
                                 pval.method = TRUE,
                                 legend = c(0.8,0.2),
                                 legend.title = "Risk Score",
                                 linetype = "solid",
                                 size = 2,
                                 ylim = c(0.3,1),
                                 xlab = "Follow up (Months)",
                                 ylab = "RFS Probability",
                                 break.time.by = NULL,
                                 risk.table = FALSE,
                                 censor = TRUE,
                                 censor.shape = 124,
                                 censor.size = 5
)

survPlot_MPRC_Train



###### Risk Score in Multi Testset


risk_test_MPRC$risk <- ifelse(risk_test_MPRC$riskscore>=Cutoff,"high","low")


####### KM-plot in Multi testset

# Model Fit
surv_model_test_MPRC <- survfit(Surv(RFS_time, RFS_status) ~ risk, data = risk_test_MPRC)
# Log-rank test
surv_diff_test_MPRC <- survdiff(Surv(RFS_time, RFS_status) ~ risk, data = risk_test_MPRC)

1 - pchisq(surv_diff_test_MPRC$chisq, length(surv_diff_test_MPRC$n) - 1) 

# KMplot
survPlot_MPRC_Test <- ggsurvplot(surv_model_test_MPRC, 
                                data = risk_test_MPRC,
                                conf.int = FALSE,  
                                palette = mypal, 
                                pval = TRUE,
                                pval.coord = c(0, 0.6),
                                pval.size = 5,
                                pval.method = TRUE,
                                legend = c(0.8,0.2),
                                legend.title = "Risk Score",
                                linetype = "solid",
                                size = 2,
                                ylim = c(0.3,1),
                                xlab = "Follow up (Months)",
                                ylab = "RFS Probability",
                                break.time.by = NULL,
                                risk.table = FALSE,
                                censor = TRUE,
                                censor.shape = 124,
                                censor.size = 5
)

survPlot_MPRC_Test



########### Feature Heatmap ############

risk_train_MPRC_tmp <- risk_train_MPRC
risk_train_MPRC_tmp$dataset <- "Train"

risk_test_MPRC_tmp <- risk_test_MPRC
risk_test_MPRC_tmp$dataset <- "Test"

risk_MPRC <- rbind(risk_train_MPRC_tmp,risk_test_MPRC_tmp)
risk_MPRC <- na.omit(risk_MPRC)
risk_MPRC_scale <- risk_MPRC

risk_MPRC_scale <- risk_MPRC_scale[order(risk_MPRC_scale$riskscore,decreasing = F),]


for (i in c(4:25)) {
  
  risk_MPRC_scale[,i] <- as.numeric(scale(risk_MPRC_scale[,i]))
}


colnames(risk_MPRC_scale)
anno_Clin <- risk_MPRC_scale[,c("riskscore","T_stage","N_stage","RFS_status","risk","dataset")]

anno_Clin$IHC_subtype <- CBCGA_Cohort.Info$`Clinical subtype`[match(rownames(anno_Clin),CBCGA_Cohort.Info$PatientCode)]
anno_Clin$PAM50_subtype <- CBCGA_Cohort.Info$`Intrinsic subtype (PAM50)`[match(rownames(anno_Clin),CBCGA_Cohort.Info$PatientCode)]
anno_Clin$PAM50_subtype[is.na(anno_Clin$PAM50_subtype)] <- "Unknown"
anno_Clin <- anno_Clin[,c("riskscore","N_stage","T_stage","IHC_subtype","PAM50_subtype","RFS_status","risk","dataset")]

risk_MPRC_feat <- as.data.frame(t(risk_MPRC_scale[,c(4:25)]))



cor_test <- c()

for (i in 4:25) {
  a <- cor.test(risk_MPRC_scale[,1],risk_MPRC_scale[,i], method = "spearman")
  cor_test <- rbind(cor_test,c(a$estimate,a$p.value))
}

cor_test <- as.data.frame(cor_test)
colnames(cor_test) <- c("Rho","p")
rownames(cor_test) <- colnames(risk_MPRC_scale)[c(4:25)]

Hm_feat_order <- rownames(cor_test)[order(cor_test$Rho,decreasing = T)]
risk_MPRC_feat_neword <- risk_MPRC_feat[Hm_feat_order,]



anno_Feat <- cor_test
anno_Feat$Feat_Cate <- c(rep("Meatb",17),rep("Path",3),rep("Rad",2))
anno_Feat <- anno_Feat[,-2]


#### Heatmap

range(risk_MPRC_feat)

negcut <- 200
poscut <- 200

mybreaks1 <- seq(-2  , 0, by = 0.01)
mybreaks2 <- seq(0.01, 2, by = 0.01)


mybreaks <- c(mybreaks1, mybreaks2)
mycolors <- c(colorRampPalette(c("#83BCD9", "white"))(negcut), colorRampPalette(c("white", "#CB4A42"))(poscut))



pheatmap(as.matrix(risk_MPRC_feat_neword), color = mycolors, breaks = mybreaks,
         annotation_col = anno_Clin, annotation_row = anno_Feat, annotation_color = anno_color, 
         show_rownames = TRUE, show_colnames = FALSE,
         fontsize = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE)


