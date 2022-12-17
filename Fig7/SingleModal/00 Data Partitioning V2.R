
# rm(list = ls())
# 
# load("D:/ShaoLab/Data Resource/CBCGA多组学数�?/CBCGA.Extended_MergedData_V2.6_220802.Rdata")
# CBCGA_PathAI <- read.csv("Data/AllSampleFeatsDf_New.csv",row.names = 1)

CBCGA_Cohort.Info$Pathomics <- ifelse(CBCGA_Cohort.Info$PatientCode%in%rownames(CBCGA_PathAI),"Yes","No")

ID_All_Diemension <- CBCGA_Cohort.Info$PatientCode[CBCGA_Cohort.Info$`RNA sequencing`=="Yes" &
                                                     CBCGA_Cohort.Info$Metabolomics=="Yes" & 
                                                     CBCGA_Cohort.Info$Radiomics=="Yes"& 
                                                     CBCGA_Cohort.Info$Pathomics=="Yes"] # 201

###### Total ID

ID_Clin <- CBCGA_Cohort.Info$PatientCode # 773
ID_RNA <- CBCGA_Cohort.Info$PatientCode[CBCGA_Cohort.Info$`RNA sequencing`=="Yes"] # 752
ID_Met <- CBCGA_Cohort.Info$PatientCode[CBCGA_Cohort.Info$Metabolomics=="Yes"] # 453
ID_Path <- CBCGA_Cohort.Info$PatientCode[CBCGA_Cohort.Info$Pathomics=="Yes"] # 626
ID_Rad <- CBCGA_Cohort.Info$PatientCode[CBCGA_Cohort.Info$Radiomics=="Yes"] # 424

ID_Clin_RNA <- intersect(ID_Clin,ID_RNA) # 752
ID_Clin_Met <- intersect(ID_Clin,ID_Met) # 453
ID_Clin_Path <- intersect(ID_Clin,ID_Path) # 626
ID_Clin_Rad <- intersect(ID_Clin,ID_Rad) # 424

ID_Clin_RNA_Met <- intersect(ID_Clin_RNA,ID_Met) # 443
ID_Clin_RNA_Path <- intersect(ID_Clin_RNA,ID_Path) # 612
ID_Clin_RNA_Rad <- intersect(ID_Clin_RNA,ID_Rad) # 416

ID_Clin_Met_Path <- intersect(ID_Clin_Met,ID_Path) # 382
ID_Clin_Met_Rad <- intersect(ID_Clin_Met,ID_Rad) # 233
ID_Clin_Path_Rad <- intersect(ID_Clin_Path,ID_Rad) # 356

ID_Clin_RNA_Met_Path <- intersect(ID_Clin_RNA_Met,ID_Path) # 375
ID_Clin_RNA_Met_Rad <- intersect(ID_Clin_RNA_Met,ID_Rad) # 230
ID_Clin_RNA_Path_Rad <- intersect(ID_Clin_RNA_Path,ID_Rad) # 350
ID_Clin_Met_Path_Rad <- intersect(ID_Clin_Met_Path,ID_Rad) # 203

ID_Clin_RNA_Met_Path_Rad <- intersect(ID_Clin_RNA_Met_Path,ID_Rad) #201

###### Testset ID

set.seed(Test_Seed)
ID_Testset_All_Dimension <- ID_All_Diemension[sample(1:length(ID_All_Diemension),SampleSize,replace = F)] # 100

##### Trainset ID

ID_Clin_Train <- setdiff(ID_Clin,ID_Testset_All_Dimension) # 673
ID_RNA_Train <- setdiff(ID_RNA,ID_Testset_All_Dimension) # 652
ID_Met_Train <- setdiff(ID_Met,ID_Testset_All_Dimension) # 353
ID_Path_Train <- setdiff(ID_Path,ID_Testset_All_Dimension) # 526
ID_Rad_Train <- setdiff(ID_Rad,ID_Testset_All_Dimension) # 324

ID_Clin_RNA_Train <- setdiff(ID_Clin_RNA,ID_Testset_All_Dimension) # 652
ID_Clin_Met_Train <- setdiff(ID_Clin_Met,ID_Testset_All_Dimension) # 353
ID_Clin_Path_Train <- setdiff(ID_Clin_Path,ID_Testset_All_Dimension) # 526
ID_Clin_Rad_Train <- setdiff(ID_Clin_Rad,ID_Testset_All_Dimension) # 324

ID_Clin_RNA_Met_Train <- setdiff(ID_Clin_RNA_Met,ID_Testset_All_Dimension) # 343
ID_Clin_RNA_Path_Train <- setdiff(ID_Clin_RNA_Path,ID_Testset_All_Dimension) # 512
ID_Clin_RNA_Rad_Train <- setdiff(ID_Clin_RNA_Rad,ID_Testset_All_Dimension) # 316

ID_Clin_Met_Path_Train <- setdiff(ID_Clin_Met_Path,ID_Testset_All_Dimension) # 282
ID_Clin_Met_Rad_Train <- setdiff(ID_Clin_Met_Rad,ID_Testset_All_Dimension) # 133
ID_Clin_Path_Rad_Train <- setdiff(ID_Clin_Path_Rad,ID_Testset_All_Dimension) # 256

ID_Clin_RNA_Met_Path_Train <- setdiff(ID_Clin_RNA_Met_Path,ID_Testset_All_Dimension) # 275
ID_Clin_RNA_Met_Rad_Train <- setdiff(ID_Clin_RNA_Met_Rad,ID_Testset_All_Dimension) # 130
ID_Clin_RNA_Path_Rad_Train <- setdiff(ID_Clin_RNA_Path_Rad,ID_Testset_All_Dimension) # 250
ID_Clin_Met_Path_Rad_Train <- setdiff(ID_Clin_Met_Path_Rad,ID_Testset_All_Dimension) # 103

ID_Clin_RNA_Met_Path_Rad_Train <- setdiff(ID_Clin_RNA_Met_Path_Rad,ID_Testset_All_Dimension) # 101






