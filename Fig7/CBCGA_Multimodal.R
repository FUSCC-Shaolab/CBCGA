rm(list = ls())

library(dplyr)
library(stringr)
library(survival)
library(survminer)
library(glmnet)
library(survivalROC)
library(pROC)
library(ggplot2)
library(ggtext)
library(ggsci)
library(rms)
library(ggrisk)
library(RColorBrewer)
library(scales)
library(data.table)
library(Hmisc)
library(boot)


for(SampleSize in c(80)) { 
  for (Test_nfolds in c(10)) { 
    for (Test_Seed in c (1)) { 
      for (lambda_Para in c("1se")) {
        

        setwd('CBCGA_Multimodal_Code_20221216/')
        load("Data/CBCGA.Extended_MergedData_V2.6_220802.Rdata")
        load("Data/All_Radiomics_Fusion.Rdata")
        CBCGA_PathAI <- read.csv("Data/AllSampleFeatsDfNoExtreDeNA.csv", row.names = 1, header = T)
        CBCGA_Clin <- read.csv("Data/CBCGA_Lum_TNBC_Clin_Omic_V220823.csv", header = T)
        CBCGA_Pathway_ssGSEA <- data.frame(fread("Data/CBCGA_TT_logTPM_ssGSEA_Hall_TME_FerrYF_rangescale.csv"),row.names = 1)
        
        ## 00 Data Partitioning
        source("SingleModal/00 Data Partitioning V2.R")

        #################################################################################################
        ## Single modal
        #################################################################################################
        ## 01 Single RNA
        source("SingleModal/01 Single RNA.R")
        Singlemodal_RNA   <- c(Train = Train, Test = Test, Num_Index = length(Active_index), Feature = Feat_Name, C_index_train, C_index_lower95_train, C_index_upper95_train, TrainCindexCI = paste0(round(C_index_lower95_train,2),"-",round(C_index_upper95_train,2)), C_index, C_index_lower95, C_index_upper95, TestCindexCI = paste0(round(C_index_lower95,2),"-",round(C_index_upper95,2)))
        print("RNA over")
        
        ## 02 Single Metab
        source("SingleModal/02 Single Metab.R")
        Singlemodal_Metab <- c(Train = Train, Test = Test, Num_Index = length(Active_index), Feature = Feat_Name, C_index_train, C_index_lower95_train, C_index_upper95_train, TrainCindexCI = paste0(round(C_index_lower95_train,2),"-",round(C_index_upper95_train,2)), C_index, C_index_lower95, C_index_upper95, TestCindexCI = paste0(round(C_index_lower95,2),"-",round(C_index_upper95,2)))
        print("Meta over")
        
        ## 03 Single Path
        source("SingleModal/03 Single Path.R")
        Singlemodal_Path  <- c(Train = Train, Test = Test, Num_Index = length(Active_index), Feature = Feat_Name, C_index_train, C_index_lower95_train, C_index_upper95_train, TrainCindexCI = paste0(round(C_index_lower95_train,2),"-",round(C_index_upper95_train,2)), C_index, C_index_lower95, C_index_upper95, TestCindexCI = paste0(round(C_index_lower95,2),"-",round(C_index_upper95,2)))
        print("Path over")
        
        ## 04 Single Rad
        source("SingleModal/04 Single Rad.R")
        Singlemodal_Rad   <- c(Train = Train, Test = Test, Num_Index = length(Active_index), Feature = Feat_Name, C_index_train, C_index_lower95_train, C_index_upper95_train, TrainCindexCI = paste0(round(C_index_lower95_train,2),"-",round(C_index_upper95_train,2)), C_index, C_index_lower95, C_index_upper95, TestCindexCI = paste0(round(C_index_lower95,2),"-",round(C_index_upper95,2)))
        print("Rad over")
        
        ## 05 Single Clin
        source("SingleModal/05 Single Clin.R")
        Singlemodal_Clin  <- c(Train = Train, Test = Test, Num_Index = 3, Feature = c("pT, pN, Subtype"), C_index_train, C_index_lower95_train, C_index_upper95_train, TrainCindexCI = paste0(round(C_index_lower95_train,2),"-",round(C_index_upper95_train,2)), C_index, C_index_lower95, C_index_upper95, TestCindexCI = paste0(round(C_index_lower95,2),"-",round(C_index_upper95,2)))
        print("Clin over")
        
        ## 06 Single IHC
        source("SingleModal/06 Single IHC.R")
        Singlemodal_IHC   <- c(Train = Train, Test = Test, Num_Index = 1, Feature = c("IHC"), C_index_train, C_index_lower95_train, C_index_upper95_train, TrainCindexCI = paste0(round(C_index_lower95_train,2),"-",round(C_index_upper95_train,2)), C_index, C_index_lower95, C_index_upper95, TestCindexCI = paste0(round(C_index_lower95,2),"-",round(C_index_upper95,2)))
        print("IHC over")
        

        Res_Singlemodal <- rbind(Singlemodal_RNA, Singlemodal_Metab, Singlemodal_Path, Singlemodal_Rad, Singlemodal_Clin, Singlemodal_IHC)

        #################################################################################################
        ## Multi modal
        #################################################################################################
        Test_Files     <-  list.files("MultiModal", pattern = "Clin.R")
        Res_Multimodal <- matrix(rep(0, length(Test_Files) * 12), ncol = 12)
        for (i in 1:length(Test_Files)) 
        {
          Test_File <- Test_Files[i]
          source(paste("MultiModal/", Test_File, sep = ""))
          Res_Multimodal[i, ] <- c(Train, Test, NA, NA, C_index_train, C_index_lower95_train, C_index_upper95_train, TrainCindexCI = paste0(round(C_index_lower95_train,2),"-",round(C_index_upper95_train,2)), C_index, C_index_lower95, C_index_upper95, TestCindexCI = paste0(round(C_index_lower95,2),"-",round(C_index_upper95,2)))
          print(i)
        }
        row.names(Res_Multimodal) <- substr(Test_Files, 5, 50)
        row.names(Res_Multimodal) <- gsub(" ", "_", row.names(Res_Multimodal), fixed = T)
        row.names(Res_Multimodal) <- paste("Multimodal", gsub(".R", "", row.names(Res_Multimodal), fixed = T), sep = "_")
        
        Res_Merge <- rbind(Res_Singlemodal, Res_Multimodal)
        colnames(Res_Merge)  <- c("Train", "Test", "Num_Index", "Feature", "C_index_train", "C_index_lower95_train", "C_index_upper95_train", "Cindex_train_95CI", "C_index_test", "C_index_lower95_test", "C_index_upper95_test", "Cindex_test_95CI")
        write.csv(Res_Merge, file = "Result/CBCGA_Multimodal_Cindex.csv")
        

        
      }
      
    }
  }
} 

