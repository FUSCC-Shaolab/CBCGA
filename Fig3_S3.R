rm(list = ls()) ; graphics.off()
#------------------------------------------------------------------------------------#
# package prepare
#------------------------------------------------------------------------------------#
library(coin)
library(ggrepel)
library(plyr)
library(tidyverse)
library(export)
library(pheatmap)
load("/data/CBCGA.Extended_MergedData_V2.5_220722.Rdata")
load("/data/Fig3_S3_data.Rdata")
#--------------------------------------------------------------------
# Fig3A
#--------------------------------------------------------------------
########### data input -------
##### CBCGA
clinic_CBCGA = CBCGA.Extended_Cohort.Info 
clinic_CBCGA = clinic_CBCGA[!is.na(clinic_CBCGA$PAM50_classifier),]
clinic_CBCGA = clinic_CBCGA[clinic_CBCGA$CBCGA_Cohort == "Yes",]
clinic_CBCGA$Histological_type = CBCGA_Cohort.Info[rownames(clinic_CBCGA),"Histological type"]
clinic_CBCGA = clinic_CBCGA[clinic_CBCGA$Histological_type == "IDC",]

##### TCGA
clinic_TCGA = clinic_TCGA[!is.na(clinic_TCGA$Race.Category),]
clinic_TCGA$IHC_subtype = "Unknown"
clinic_TCGA$IHC_subtype[(clinic_TCGA$ER_Final == "Positive" |clinic_TCGA$PR_Final == "Positive") & clinic_TCGA$HER2_Final == "Negative"] = "HR+HER2-"
clinic_TCGA$IHC_subtype[(clinic_TCGA$ER_Final == "Positive" |clinic_TCGA$PR_Final == "Positive") & clinic_TCGA$HER2_Final == "Positive"] = "HR+HER2+"
clinic_TCGA$IHC_subtype[(clinic_TCGA$ER_Final == "Negative" & clinic_TCGA$PR_Final == "Negative") & clinic_TCGA$HER2_Final == "Positive"] = "HR-HER2+"
clinic_TCGA$IHC_subtype[(clinic_TCGA$ER_Final == "Negative" & clinic_TCGA$PR_Final == "Negative") & clinic_TCGA$HER2_Final == "Negative"] = "TNBC"
clinic_TCGA = clinic_TCGA[clinic_TCGA$IHC_subtype != "Unknown",]

clinic_TCGA = clinic_TCGA[clinic_TCGA$histological_type == "Infiltrating Ductal Carcinoma" , ]
clinic_TCGA = clinic_TCGA[substring(rownames(clinic_TCGA),14,15) == "01",]

rownames(clinic_TCGA) = substring(rownames(clinic_TCGA),1,12)
clinic_TCGA = clinic_TCGA[clinic_TCGA$PAM50Call_RNAseq != "",]

clinic_TCGA_Asian = clinic_TCGA_white = clinic_TCGA

clinic_TCGA_white = clinic_TCGA[clinic_TCGA$Race.Category == "WHITE" , ]
clinic_TCGA_Asian = clinic_TCGA[clinic_TCGA$Race.Category == "ASIAN", ]
clinic_TCGA = rbind(clinic_TCGA_white,clinic_TCGA_Asian)

########### draw all patient ---------
mat = data.frame(row.names = c(rownames(clinic_TCGA),rownames(clinic_CBCGA)) ,
                 race = c(clinic_TCGA$Race.Category,rep("Asian_CBCGA",nrow(clinic_CBCGA))),
                 PAM50 = c(clinic_TCGA$PAM50Call_RNAseq, clinic_CBCGA$PAM50_classifier)
)
cluster.merge4plot = as.data.frame(table(mat$race, mat$PAM50))
colnames(cluster.merge4plot)[1:2] = c("race","PAM50")

cluster.merge4plot2 = ddply(cluster.merge4plot,"race", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$PAM50 = as.character(cluster.merge4plot2$PAM50)
cluster.merge4plot2$race = as.character(cluster.merge4plot2$race)

cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Her2"] = "HER2-enriched"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "LumB"] = "Luminal B"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Basal"] = "Basal-like"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Normal"] = "Normal-like"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "LumA"] = "Luminal A"

cluster.merge4plot2$race[cluster.merge4plot2$race == "ASIAN"] = paste0("TCGA Asian (N = ",nrow(clinic_TCGA_Asian),")")
cluster.merge4plot2$race[cluster.merge4plot2$race == "WHITE"] = paste0("TCGA Caucasian (N = ",nrow(clinic_TCGA_white),")")
cluster.merge4plot2$race[cluster.merge4plot2$race == "Asian_CBCGA"] = paste0("CBCGA (N = ",nrow(clinic_CBCGA),")")

cluster.merge4plot2$PAM50 = factor(cluster.merge4plot2$PAM50, levels = c("Normal-like","Luminal B","Luminal A","HER2-enriched","Basal-like"))

color = c("#1D76BC","#76CFE6","#6E59A6","#E11D2E","#CDCFD0")
names(color) = c("Luminal A","Luminal B","HER2-enriched","Basal-like","Normal-like")

tmp = fisher.test(table(mat$race, mat$PAM50),simulate.p.value=TRUE)
p = ggplot(cluster.merge4plot2, aes(x=race, y = percent, fill = PAM50 )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5), 
        axis.ticks = element_line(size = 0.5),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",tmp[["p.value"]]) ) 
graph2pdf(p, file = "/results/Fig3A_part1.pdf")

# write.csv(cluster.merge4plot2,file ="/results/race_PAM50_percent_barplot.csv")

########### draw subtype ---------
mat0 = data.frame(row.names = c(rownames(clinic_TCGA_white),rownames(clinic_CBCGA)) ,
                  race = c(clinic_TCGA_white$Race.Category,rep("Asian_CBCGA",nrow(clinic_CBCGA))),
                  PAM50 = c(clinic_TCGA_white$PAM50Call_RNAseq, clinic_CBCGA$PAM50_classifier),
                  IHC_subtype = c(clinic_TCGA_white$IHC_subtype, clinic_CBCGA$Clinical_Subtype)
)

for( i in names(table(mat0$IHC_subtype))){
  mat = mat0[mat0$IHC_subtype == i,]
  cluster.merge4plot = as.data.frame(table(mat$race, mat$PAM50))
  colnames(cluster.merge4plot)[1:2] = c("race","PAM50")
  
  cluster.merge4plot2 = ddply(cluster.merge4plot,"race", transform,
                              percent = Freq / sum(Freq) *100)
  cluster.merge4plot2$PAM50 = as.character(cluster.merge4plot2$PAM50)
  cluster.merge4plot2$race = as.character(cluster.merge4plot2$race)
  
  cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Her2"] = "HER2-enriched"
  cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "LumB"] = "Luminal B"
  cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Basal"] = "Basal-like"
  cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Normal"] = "Normal-like"
  cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "LumA"] = "Luminal A"
  
  cluster.merge4plot2$race[cluster.merge4plot2$race == "WHITE"] = 
    paste0("TCGA Caucasian (N = ", nrow(mat[mat$race == "WHITE" ,]) ,")")
  cluster.merge4plot2$race[cluster.merge4plot2$race == "Asian_CBCGA"] = 
    paste0("CBCGA (N = ", nrow( mat[mat$race == "Asian_CBCGA", ]) ,")")
  
  cluster.merge4plot2$PAM50 = factor(cluster.merge4plot2$PAM50, levels = c("Normal-like","Luminal B","Luminal A","HER2-enriched","Basal-like"))
  
  color = c("#1D76BC","#76CFE6","#6E59A6","#E11D2E","#CDCFD0")
  names(color) = c("Luminal A","Luminal B","HER2-enriched","Basal-like","Normal-like")
  
  tmp = fisher.test(table(mat$race, mat$PAM50),simulate.p.value=TRUE)
  p = ggplot(cluster.merge4plot2, aes(x=race, y = percent, fill = PAM50 )) + geom_bar(stat = "identity") +
    scale_fill_manual(values = color)+
    theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.text.x=element_text(size = 10),
          axis.text.y=element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    ) +
    ggtitle(paste0("p = ",tmp[["p.value"]]) ) 
  graph2pdf(p, file = paste0("/results/Fig3A_part2_",i,".pdf"))
  
  #write.csv(cluster.merge4plot2,file =paste0("/results/",i,"_race_PAM50_percent_barplot.csv"))
}
### change the P-value to FDR manually 

#--------------------------------------------------------------------
# Fig3B-C, S3B-C part1: the whole population
#--------------------------------------------------------------------
load("/data/Fig3_S3_data.Rdata")
########### data input -------
##### CBCGA
thre = CBCGA_GISTICgene.thre
gene.loc = gene.loc.noX
thre = thre[rownames(gene.loc),]

clinic_CBCGA = CBCGA_Cohort.Info
clinic_CBCGA = clinic_CBCGA[intersect(rownames(clinic_CBCGA), colnames(thre)),]
clinic_CBCGA = clinic_CBCGA[clinic_CBCGA$`Histological type` == "IDC",]

subtype_CBCGA = list(HRposHER2neg = rownames(clinic_CBCGA)[(clinic_CBCGA$`ER status` == "Positive" | clinic_CBCGA$`PR status` == "Positive") & clinic_CBCGA$`HER2 status (combined)` == "Negative"],
                     HRposHER2pos = rownames(clinic_CBCGA)[(clinic_CBCGA$`ER status` == "Positive" | clinic_CBCGA$`PR status` == "Positive") & clinic_CBCGA$`HER2 status (combined)` == "Positive"],
                     HRnegHER2pos = rownames(clinic_CBCGA)[(clinic_CBCGA$`ER status` == "Negative" & clinic_CBCGA$`PR status` == "Negative") & clinic_CBCGA$`HER2 status (combined)` == "Positive"],
                     HRnegHER2neg = rownames(clinic_CBCGA)[ clinic_CBCGA$`ER status` == "Negative" & clinic_CBCGA$`PR status` == "Negative" & clinic_CBCGA$`HER2 status (combined)` == "Negative"]
)

##### TCGA
clinic_TCGA = clinic_TCGA[clinic_TCGA$Race.Category == "WHITE" & !is.na(clinic_TCGA$Race.Category), ]
clinic_TCGA = clinic_TCGA[substring(rownames(clinic_TCGA),14,15) == "01",]
rownames(clinic_TCGA) = substring(rownames(clinic_TCGA),1,12)
clinic_TCGA = clinic_TCGA[clinic_TCGA$histological_type == "Infiltrating Ductal Carcinoma", ]

tmp = intersect(rownames(TCGA_cnv),rownames(thre))
TCGA_cnv = TCGA_cnv[tmp,]
thre = thre[tmp,]

tmp = intersect(rownames(clinic_TCGA),colnames(TCGA_cnv))
TCGA_cnv = TCGA_cnv[,tmp]
clinic_TCGA = clinic_TCGA[tmp,]

subtype_TCGA = list(HRposHER2neg = na.omit(rownames(clinic_TCGA)[(clinic_TCGA$ER_Final == "Positive" | clinic_TCGA$PR_Final == "Positive") & clinic_TCGA$HER2_Final == "Negative"]),
                    HRposHER2pos = na.omit(rownames(clinic_TCGA)[(clinic_TCGA$ER_Final == "Positive" | clinic_TCGA$PR_Final == "Positive") & clinic_TCGA$HER2_Final == "Positive"]),
                    HRnegHER2pos = na.omit(rownames(clinic_TCGA)[(clinic_TCGA$ER_Final == "Negative" & clinic_TCGA$PR_Final == "Negative") & clinic_TCGA$HER2_Final == "Positive"]),
                    HRnegHER2neg = na.omit(rownames(clinic_TCGA)[(clinic_TCGA$ER_Final == "Negative" & clinic_TCGA$PR_Final == "Negative") & clinic_TCGA$HER2_Final == "Negative"])
)


########### DCG & draw ---------
CBCGA = rownames(clinic_CBCGA)
TCGA = rownames(clinic_TCGA)

tmp = intersect(cnvlist,rownames(TCGA_cnv))
TCGA_cnv_amp = TCGA_cnv_loss = TCGA_cnv[tmp,TCGA]
thre_amp = thre_loss = thre[tmp,CBCGA]

for (i in rownames(TCGA_cnv_amp)){
  mat1_loss = mat1_amp = TCGA_cnv_amp[i,]
  mat1_amp[mat1_amp != 2 ] = 0 ; mat1_amp[mat1_amp == 2 ] = 1
  mat1_loss[mat1_loss != -2 ] = 0 ; mat1_loss[mat1_loss == -2 ] = 1
  
  TCGA_cnv_amp[i,] = mat1_amp
  TCGA_cnv_loss[i,] = mat1_loss
  
  mat2_loss = mat2_amp = thre_amp[i,]
  mat2_amp[mat2_amp != 2 ] = 0 ; mat2_amp[mat2_amp == 2 ] = 1
  mat2_loss[mat2_loss != -2 ] = 0 ; mat2_loss[mat2_loss == -2 ] = 1
  
  thre_amp[i,] = mat2_amp
  thre_loss[i,] = mat2_loss
  
  print(i)
}

##### compare ------
CBCGAVStcga.loss.res = CBCGAVStcga.gain.res = as.data.frame( matrix(nrow = nrow(TCGA_cnv_amp), ncol = 6) )
rownames(CBCGAVStcga.loss.res) = rownames(CBCGAVStcga.gain.res)  = row.names(TCGA_cnv_amp)
colnames(CBCGAVStcga.loss.res) = colnames(CBCGAVStcga.gain.res) = c("CBCGA","TCGA","p.val","adj.p","CBCGA_case","TCGA_case")

CBCGAVStcga.gain.res$CBCGA = apply(thre_amp,1,function(x){sum(x)/(length(x))}) 
CBCGAVStcga.gain.res$CBCGA_case = apply(thre_amp,1,function(x){sum(x)}) 
CBCGAVStcga.gain.res$TCGA = apply(TCGA_cnv_amp,1,function(x){sum(x)/(length(x))})
CBCGAVStcga.gain.res$TCGA_case = apply(TCGA_cnv_amp,1,function(x){sum(x)})

CBCGAVStcga.loss.res$CBCGA = apply(thre_loss,1,function(x){sum(x)/(length(x))})
CBCGAVStcga.loss.res$CBCGA_case = apply(thre_loss,1,function(x){sum(x)})
CBCGAVStcga.loss.res$TCGA = apply(TCGA_cnv_loss,1,function(x){sum(x)/(length(x))})
CBCGAVStcga.loss.res$TCGA_case = apply(TCGA_cnv_loss,1,function(x){sum(x)})

## gain DCG -------
for ( i in rownames(CBCGAVStcga.gain.res)){
  tmp1 = t(rbind(rep("CBCGA",ncol(thre_amp)), thre_amp[i,]))
  tmp2= t(rbind(rep("TCGA",ncol(TCGA_cnv_amp)), TCGA_cnv_amp[i,]))
  tmp = as.data.frame(rbind(tmp1,tmp2))
  colnames(tmp) = c("cohort","gene")
  mytab = table(tmp$cohort,tmp$gene)
  if (dim(mytab)[2] != 1 ){
    #if( all(as.data.frame(mytab)[,3] > 5) ){p = chisq.test(mytab)} else {p = fisher.test(mytab)}
    p = fisher.test(mytab)
    CBCGAVStcga.gain.res[i,3] <- p$p.value
  }
}
CBCGAVStcga.gain.res[,4] = p.adjust(CBCGAVStcga.gain.res[,3],method = "fdr")
CBCGAVStcga.gain.res = CBCGAVStcga.gain.res[!is.na(CBCGAVStcga.gain.res$p.val),]
#write.csv(CBCGAVStcga.gain.res,file = paste0("/results/CBCGAVStcga_WHITE_gain_res.csv"))

## loss DCG -------
for ( i in rownames(CBCGAVStcga.loss.res)){
  tmp1 = t(rbind(rep("CBCGA",ncol(thre_loss)), thre_loss[i,]))
  tmp2= t(rbind(rep("TCGA",ncol(TCGA_cnv_loss)), TCGA_cnv_loss[i,]))
  tmp = as.data.frame(rbind(tmp1,tmp2))
  colnames(tmp) = c("cohort","gene")
  mytab = table(tmp$cohort,tmp$gene)
  if (dim(mytab)[2] != 1 ){
    #if( all(as.data.frame(mytab)[,3] > 5) ){p = chisq.test(mytab)} else {p = fisher.test(mytab)}
    p = fisher.test(mytab)
    CBCGAVStcga.loss.res[i,3] <- p$p.value
  }
}
CBCGAVStcga.loss.res[,4] = p.adjust(CBCGAVStcga.loss.res[,3],method = "fdr")
CBCGAVStcga.loss.res = CBCGAVStcga.loss.res[!is.na(CBCGAVStcga.loss.res$p.val),]
#write.csv(CBCGAVStcga.loss.res,file = paste0("/results/CBCGAVStcga_WHITE_loss_res.csv"))

## gain plot -------
#CBCGAVStcga.gain.res = read.csv(paste0("/results/CBCGAVStcga_WHITE_gain_res.csv"),row.names = 1)
p = ggplot(CBCGAVStcga.gain.res,aes(x=CBCGA, y = TCGA)) + 
  geom_point(alpha = 1, size = 2.5, col="#cecece") + 
  geom_point(data = CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05,] , 
             alpha = 1, size = 2.5, col="#C0052A") +
  geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
  xlim(0,1) + ylim(0,1) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
  geom_text_repel(data= CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05 &
                                               (CBCGAVStcga.gain.res$CBCGA > 0.1 & CBCGAVStcga.gain.res$TCGA >0.1),] ,
                  aes(label=rownames(CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05 &
                                                            (CBCGAVStcga.gain.res$CBCGA > 0.1 & CBCGAVStcga.gain.res$TCGA >0.1),])),
                  col="black",alpha = 1,size=3) 
ggsave(p,filename = paste0("/results/Fig3B_CBCGAvsTCGA_WHITE_gain.pdf"),height = 4,width = 4)
#export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_gain.pdf"),height = 4,width = 4)

p = ggplot(CBCGAVStcga.gain.res,aes(x=CBCGA, y = TCGA)) + 
  geom_point(alpha = 1, size = 2.5, col="#cecece") + 
  geom_point(data = CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05,] , 
             alpha = 1, size = 2.5, col="#C0052A") +
  geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
  xlim(0,0.25) + ylim(0,0.25) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
  geom_text_repel(data= CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05 & 
                                               (CBCGAVStcga.gain.res$CBCGA > 0.1 | CBCGAVStcga.gain.res$TCGA > 0.1),] ,
                  aes(label=rownames(CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05 &
                                                            (CBCGAVStcga.gain.res$CBCGA > 0.1 | CBCGAVStcga.gain.res$TCGA > 0.1),])),
                  col="black",alpha = 1,size=3) 
ggsave(p,filename = paste0("/results/Fig3B_CBCGAvsTCGA_WHITE_gain2.pdf"),height = 4,width = 4)
#export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_gain2.pdf"),height = 4,width = 4)

## loss plot-------
#CBCGAVStcga.loss.res = read.csv(paste0("/results/CBCGAVStcga_WHITE_loss_res.csv"),row.names = 1)
p = ggplot(CBCGAVStcga.loss.res,aes(x=CBCGA, y = TCGA)) + 
  geom_point(alpha = 1, size = 2.5, col="#cecece") + 
  geom_point(data = CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05,] , 
             alpha = 1, size = 2.5, col="#4E85AC") +
  geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
  xlim(0,1) + ylim(0,1) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
  geom_text_repel(data= CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 & 
                                               (CBCGAVStcga.loss.res$CBCGA > 0.1 | CBCGAVStcga.loss.res$TCGA > 0.1 ),] ,
                  aes(label=rownames(CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 &
                                                            (CBCGAVStcga.loss.res$CBCGA > 0.1 | CBCGAVStcga.loss.res$TCGA > 0.1 ),])),
                  col="black",alpha = 1,size=3,
                  max.overlaps = 20) 
ggsave(p,filename = paste0("/results/Fig3B_CBCGAvsTCGA_WHITE_loss.pdf"),height = 4,width = 4)
#export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_loss.pdf") ,height = 4,width = 4)

p = ggplot(CBCGAVStcga.loss.res,aes(x=CBCGA, y = TCGA)) + 
  geom_point(alpha = 1, size = 2.5, col="#cecece") + 
  geom_point(data = CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05,] , 
             alpha = 1, size = 2.5, col="#4E85AC") +
  geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
  xlim(0,0.25) + ylim(0,0.25) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
  geom_text_repel(data= CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 & 
                                               (CBCGAVStcga.loss.res$CBCGA > 0.1 | CBCGAVStcga.loss.res$TCGA > 0.1 ),] ,
                  aes(label=rownames(CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 &
                                                            (CBCGAVStcga.loss.res$CBCGA > 0.1 | CBCGAVStcga.loss.res$TCGA > 0.1 ),])),
                  col="black",alpha = 1,size=3,
                  max.overlaps = 20) 
ggsave(p,filename = paste0("/results/Fig3B_CBCGAvsTCGA_WHITE_loss2.pdf"),height = 4,width = 4)
#export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_loss2.pdf") ,height = 4,width = 4)

#--------------------------------------------------------------------
# Fig3B-C, S3B-C part1: IHC subtype
#--------------------------------------------------------------------
for ( k in c("HRposHER2neg","HRposHER2pos","HRnegHER2pos","HRnegHER2neg")){
  
  CBCGA = subtype_CBCGA[[k]]
  TCGA = subtype_TCGA[[k]]
  
  tmp = intersect(cnvlist,rownames(TCGA_cnv))
  TCGA_cnv_amp = TCGA_cnv_loss = TCGA_cnv[tmp,TCGA]
  thre_amp = thre_loss = thre[tmp,CBCGA]
  
  for (i in rownames(TCGA_cnv_amp)){
    mat1_loss = mat1_amp = TCGA_cnv_amp[i,]
    mat1_amp[mat1_amp != 2 ] = 0 ; mat1_amp[mat1_amp == 2 ] = 1
    mat1_loss[mat1_loss != -2 ] = 0 ; mat1_loss[mat1_loss == -2 ] = 1
    
    TCGA_cnv_amp[i,] = mat1_amp
    TCGA_cnv_loss[i,] = mat1_loss
    
    mat2_loss = mat2_amp = thre_amp[i,]
    mat2_amp[mat2_amp != 2 ] = 0 ; mat2_amp[mat2_amp == 2 ] = 1
    mat2_loss[mat2_loss != -2 ] = 0 ; mat2_loss[mat2_loss == -2 ] = 1
    
    thre_amp[i,] = mat2_amp
    thre_loss[i,] = mat2_loss
    
    print(i)
  }
  
  #####compare
  CBCGAVStcga.loss.res = CBCGAVStcga.gain.res = as.data.frame( matrix(nrow = nrow(TCGA_cnv_amp), ncol = 6) )
  rownames(CBCGAVStcga.loss.res) = rownames(CBCGAVStcga.gain.res)  = row.names(TCGA_cnv_amp)
  colnames(CBCGAVStcga.loss.res) = colnames(CBCGAVStcga.gain.res) = c("CBCGA","TCGA","p.val","adj.p","CBCGA_case","TCGA_case")
  
  CBCGAVStcga.gain.res$CBCGA = apply(thre_amp,1,function(x){sum(x)/(length(x))}) 
  CBCGAVStcga.gain.res$CBCGA_case = apply(thre_amp,1,function(x){sum(x)}) 
  CBCGAVStcga.gain.res$TCGA = apply(TCGA_cnv_amp,1,function(x){sum(x)/(length(x))})
  CBCGAVStcga.gain.res$TCGA_case = apply(TCGA_cnv_amp,1,function(x){sum(x)})
  
  CBCGAVStcga.loss.res$CBCGA = apply(thre_loss,1,function(x){sum(x)/(length(x))})
  CBCGAVStcga.loss.res$CBCGA_case = apply(thre_loss,1,function(x){sum(x)})
  CBCGAVStcga.loss.res$TCGA = apply(TCGA_cnv_loss,1,function(x){sum(x)/(length(x))})
  CBCGAVStcga.loss.res$TCGA_case = apply(TCGA_cnv_loss,1,function(x){sum(x)})
  
  ## gain
  for ( i in rownames(CBCGAVStcga.gain.res)){
    tmp1 = t(rbind(rep("CBCGA",ncol(thre_amp)), thre_amp[i,]))
    tmp2= t(rbind(rep("TCGA",ncol(TCGA_cnv_amp)), TCGA_cnv_amp[i,]))
    tmp = as.data.frame(rbind(tmp1,tmp2))
    colnames(tmp) = c("cohort","gene")
    mytab = table(tmp$cohort,tmp$gene)
    if (dim(mytab)[2] != 1 ){
      #if( all(as.data.frame(mytab)[,3] > 5) ){p = chisq.test(mytab)} else {p = fisher.test(mytab)}
      p = fisher.test(mytab)
      CBCGAVStcga.gain.res[i,3] <- p$p.value
    }
  }
  CBCGAVStcga.gain.res[,4] = p.adjust(CBCGAVStcga.gain.res[,3],method = "fdr")
  CBCGAVStcga.gain.res = CBCGAVStcga.gain.res[!is.na(CBCGAVStcga.gain.res$p.val),]
  #write.csv(CBCGAVStcga.gain.res,file = paste0("/results/CBCGAVStcga_WHITE_gain_res_",k,".csv"))
  
  ## loss
  for ( i in rownames(CBCGAVStcga.loss.res)){
    tmp1 = t(rbind(rep("CBCGA",ncol(thre_loss)), thre_loss[i,]))
    tmp2= t(rbind(rep("TCGA",ncol(TCGA_cnv_loss)), TCGA_cnv_loss[i,]))
    tmp = as.data.frame(rbind(tmp1,tmp2))
    colnames(tmp) = c("cohort","gene")
    mytab = table(tmp$cohort,tmp$gene)
    if (dim(mytab)[2] != 1 ){
      #if( all(as.data.frame(mytab)[,3] > 5) ){p = chisq.test(mytab)} else {p = fisher.test(mytab)}
      p = fisher.test(mytab)
      CBCGAVStcga.loss.res[i,3] <- p$p.value
    }
  }
  CBCGAVStcga.loss.res[,4] = p.adjust(CBCGAVStcga.loss.res[,3],method = "fdr")
  CBCGAVStcga.loss.res = CBCGAVStcga.loss.res[!is.na(CBCGAVStcga.loss.res$p.val),]
  #write.csv(CBCGAVStcga.loss.res,file = paste0("/results/CBCGAVStcga_WHITE_loss_res_",k,".csv"))
  
  
  #### draw plot
  #CBCGAVStcga.gain.res = read.csv(paste0("/results/CBCGAVStcga_WHITE_gain_res_",k,".csv"),row.names = 1)
  p = ggplot(CBCGAVStcga.gain.res,aes(x=CBCGA, y = TCGA)) + 
    geom_point(alpha = 1, size = 2.5, col="#cecece") + 
    geom_point(data = CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05,] , 
               alpha = 1, size = 2.5, col="#C0052A") +
    geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
    xlim(0,1) + ylim(0,1) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
    geom_text_repel(data= CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05 & 
                                                 (CBCGAVStcga.gain.res$CBCGA > 0.25 | CBCGAVStcga.gain.res$TCGA > 0.25 ),] ,
                    aes(label=rownames(CBCGAVStcga.gain.res[CBCGAVStcga.gain.res$adj.p < 0.05  & 
                                                              (CBCGAVStcga.gain.res$CBCGA > 0.25 | CBCGAVStcga.gain.res$TCGA > 0.25 ),] )
                    ),
                    col="black",alpha = 1,size=3,
                    max.overlaps = 20) 
  ggsave(p,filename = paste0("/results/Fig3C_S3_CBCGAvsTCGA_WHITE_gain_",k,".pdf"),height = 4,width = 4)
  #export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_gain_",k,".pdf"),height = 4,width = 4)
  
  #CBCGAVStcga.loss.res = read.csv(paste0("/results/CBCGAVStcga_WHITE_loss_res_",k,".csv"),row.names = 1)
  p = ggplot(CBCGAVStcga.loss.res,aes(x=CBCGA, y = TCGA)) + 
    geom_point(alpha = 1, size = 2.5, col="#cecece") + 
    geom_point(data = CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05,] , 
               alpha = 1, size = 2.5, col="#4E85AC") +
    geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.2, color = '#747475')+
    xlim(0,1) + ylim(0,1) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1.2)) +
    geom_text_repel(data= CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 & 
                                                 (CBCGAVStcga.loss.res$CBCGA > 0.25 | CBCGAVStcga.loss.res$TCGA > 0.25 ),] ,
                    aes(label=rownames(CBCGAVStcga.loss.res[CBCGAVStcga.loss.res$adj.p < 0.05 &
                                                              (CBCGAVStcga.loss.res$CBCGA > 0.25 | CBCGAVStcga.loss.res$TCGA > 0.25 ),])),
                    col="black",alpha = 1,size=3,
                    max.overlaps = 20) 
  ggsave(p,filename = paste0("/results/Fig3C_S3_CBCGAvsTCGA_WHITE_loss_",k,".pdf"),height = 4,width = 4)
  #export::graph2ppt(p,file = paste0("/results/CBCGAvsTCGA_WHITE_loss_",k,".pdf"),height = 4,width = 4)
}



#--------------------------------------------------------------------
# Fig3D
#--------------------------------------------------------------------
Tes_Cases <- CBCGA.Extended_Cohort.Info$PatientCode[CBCGA.Extended_Cohort.Info$CBCGA_Cohort == "Yes"]
HRpHER2p_Cases <- CBCGA.Extended_Cohort.Info$PatientCode[CBCGA.Extended_Cohort.Info$Clinical_Subtype == "HR+HER2+" ]

EXP_Tumor <- CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode[, str_detect(colnames(CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode), fixed("_T"))]
colnames(EXP_Tumor) <- substr( colnames(EXP_Tumor), 1, 4)
Tes_EXP   <- EXP_Tumor[, colnames(EXP_Tumor) %in% Tes_Cases]
Tes_EXP <- log2(Tes_EXP + 1)

Tes_PRO <- CBCGA.Extended_PRO_normalized
Tes_PRO <- Tes_PRO[, str_detect(colnames(Tes_PRO), fixed("_T"))]
colnames(Tes_PRO) <- substr( colnames(Tes_PRO), 1, 4)
Tes_PRO <- Tes_PRO[, colnames(Tes_PRO) %in% Tes_Cases]

Tes_CNA <- CBCGA_GISTICgene.alldata

### Generate Matrix
SCNA_Matrix <- Tes_CNA[c("NF1", "ACACA", "LASP1", "CDK12", "STARD3", "ERBB2", "GRB7", "MED24", "TOP2A"), ]
# SCNA_Matrix <- apply(SCNA_Matrix, 1, scale.default)
SCNA_Matrix[SCNA_Matrix > 5] <- 5

mRNA_Matrix <- Tes_EXP[c("NF1", "ACACA", "LASP1", "CDK12", "STARD3", "ERBB2", "GRB7",  "MED24", "TOP2A"), ]
mRNA_Matrix <- t(apply(t(mRNA_Matrix), 2, scale.default))
colnames(mRNA_Matrix) <- colnames(Tes_EXP)

PRO_Matrix  <- Tes_PRO[c("NF1", "ACACA", "LASP1", "CDK12", "STARD3", "ERBB2", "GRB7", "MED24", "TOP2A"), ]
PRO_Matrix  <- t(apply(t(PRO_Matrix), 2, scale.default))
colnames(PRO_Matrix) <- colnames(Tes_PRO)

Tes_Cases <- intersect(Tes_Cases, HRpHER2p_Cases)
Tes_Cases <- intersect( Tes_Cases, colnames(SCNA_Matrix) )
Tes_Cases <- intersect( Tes_Cases, colnames(mRNA_Matrix) )
Tes_Cases <- intersect( Tes_Cases, colnames(PRO_Matrix) ) 

Tes_CNA   <- SCNA_Matrix[, Tes_Cases]
Tes_Cases <- Tes_Cases[order(Tes_CNA["ERBB2", ], decreasing = T)]

Tes_CNA   <- Tes_CNA[, Tes_Cases]
Tes_EXP   <- EXP_Tumor[, Tes_Cases]
Tes_PRO   <- Tes_PRO[, Tes_Cases]

Clin_Matrix <- CBCGA.Extended_Cohort.Info[, c("PAM50_classifier")]

SCNA_Matrix <- SCNA_Matrix[, Tes_Cases]
mRNA_Matrix <- mRNA_Matrix[, Tes_Cases]
PRO_Matrix  <- PRO_Matrix[, Tes_Cases]

Pheatmap.breaks <- function(matrix, type, n.cut)
{
  pos    <- max(matrix) + 0.1 ; neg    <- min(matrix) + 0.1
  poscut <- n.cut
  negcut <- n.cut
  if(type == 1) 
  {mybreaks1 <- seq(neg, 0, (-neg)/negcut); mybreaks2 <- seq(0, pos, pos/poscut)[-1]
  } else {
    mybreaks1 <- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^2
    mybreaks2 <- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^2
  } 
  mybreaks <- c(mybreaks1, mybreaks2)
  return(mybreaks)
}

Anno <- CBCGA_Cohort.Info[colnames(PRO_Matrix),  c("PR positive percentage",  "ER positive percentage", "Intrinsic subtype (PAM50)")]
colnames(Anno) <- c("PR", "ER", "PAM50")
Anno$PR[is.na(Anno$PR)] = 0
Anno$PR_percentage = Anno$PR
Anno$ER_percentage = Anno$ER

Anno$PR_percentage[Anno$PR == 0] = "0"
Anno$PR_percentage[Anno$PR < 10 & Anno$PR > 0 ] = "1-9"
Anno$PR_percentage[Anno$PR >= 10 & Anno$PR < 30] = "10-29"
Anno$PR_percentage[Anno$PR >= 30 & Anno$PR < 70] = "30-69"
Anno$PR_percentage[Anno$PR >= 70 ] = "70-100"

Anno$ER_percentage[Anno$ER < 10] = "1-9"
Anno$ER_percentage[Anno$ER >= 10 & Anno$ER < 30] = "10-29"
Anno$ER_percentage[Anno$ER >= 30 & Anno$ER < 70] = "30-69"
Anno$ER_percentage[Anno$ER >= 70 ] = "70-100"
Anno = Anno[,c("PR_percentage", "ER_percentage", "PAM50")]

Anno_Color <- list(
  PR_percentage =  c( "70-100" = "#000000", "30-69"  = "#666666","10-29" = "#999999" , "1-9" =  "#CCCCCC","0" = "#FFFFFF"), 
  ER_percentage = c( "70-100" = "#000000", "30-69"  = "#666666","10-29" = "#999999" , "1-9" =  "#CCCCCC"), 
  PAM50 = Color_PAM50
)


### Plots
## SCNA
n.cut        <- 50
mybreaks_CNA <- Pheatmap.breaks(SCNA_Matrix, 1, n.cut)
mycolor_CNA  <- c(colorRampPalette(c("#3A53A4", "#FFFFFF"))(n.cut), colorRampPalette(c( "#FFFFFF", "#E21E26"))(n.cut))

pheatmap(SCNA_Matrix, cluster_rows = F, cluster_cols = F, color = mycolor_CNA, breaks = mybreaks_CNA, border_color = NA, 
         annotation_col = Anno, annotation_colors = Anno_Color,
         filename = "/results/Fig3D_CNA.pdf", width = 10, height = 3.4)
## mRNA
n.cut        <- 50
mybreaks_mRNA <- Pheatmap.breaks(mRNA_Matrix, 2, n.cut)
mycolor_mRNA  <- c(colorRampPalette(c("Green", "#000000"))(n.cut), colorRampPalette(c( "#000000", "red"))(n.cut))

pheatmap(mRNA_Matrix, cluster_rows = F, cluster_cols = F, color = mycolor_mRNA, breaks = mybreaks_mRNA, border_color = NA,
         annotation_col = Anno, annotation_colors = Anno_Color,
         filename = "/results/Fig3D_mRNA.pdf", width = 10, height = 3.4)


## Protein
n.cut        <- 50
mybreaks_Pro <- Pheatmap.breaks(PRO_Matrix, 1, n.cut)
mycolor_Pro  <- c(colorRampPalette(c("#3A53A4", "#FFFFFF"))(n.cut), colorRampPalette(c( "#FFFFFF", "#E21E26"))(n.cut))
pheatmap(PRO_Matrix, cluster_rows = F, cluster_cols = F, color = mycolor_Pro, breaks = mybreaks_Pro, border_color = NA,
         annotation_col = Anno, annotation_colors = Anno_Color,
         filename = "/results/Fig3D_Pro.pdf", width = 10, height = 3.4)














