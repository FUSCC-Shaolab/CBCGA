
# CBCGA_ID <- CBCGA_Clin$PatientCode[CBCGA_Clin$CBCGA=="Yes"]
# CBCGA_hospital_ID <- CBCGA_Clin$Hospital_ID[match(CBCGA_ID,CBCGA_Clin$PatientCode)]


OutCome_IHC <- data.frame(RFS_time = CBCGA_Cohort.Info$`RFS time (month)`,
                          RFS_status = CBCGA_Cohort.Info$`RFS status`,
                          IHC = CBCGA_Cohort.Info$`Clinical subtype`)
rownames(OutCome_IHC) <- CBCGA_Cohort.Info$PatientCode

OutCome_IHC <- na.omit(OutCome_IHC)

table(OutCome_IHC$IHC)



# OutCome_IHC$Stage <- case_when(OutCome_IHC$T_stage=="pT1" & OutCome_IHC$N_stage=="pN0"~"I",
#                                 (OutCome_IHC$T_stage=="pT0" & OutCome_IHC$N_stage=="pN1")|
#                                   (OutCome_IHC$T_stage=="pT1" & OutCome_IHC$N_stage=="pN1")|
#                                   (OutCome_IHC$T_stage=="pT2" & OutCome_IHC$N_stage=="pN0")|
#                                   (OutCome_IHC$T_stage=="pT2" & OutCome_IHC$N_stage=="pN1")|
#                                   (OutCome_IHC$T_stage=="pT3" & OutCome_IHC$N_stage=="pN0")~"II",
#                                 TRUE~"III")
# table(OutCome_IHC$Stage)
# 
# OutCome_IHC$Stage_2 <- case_when(OutCome_IHC$T_stage=="pT1" & OutCome_IHC$N_stage=="pN0"~"IA",
#                                   (OutCome_IHC$T_stage=="pT0" & OutCome_IHC$N_stage=="pN1")|
#                                     (OutCome_IHC$T_stage=="pT1" & OutCome_IHC$N_stage=="pN1")|
#                                     (OutCome_IHC$T_stage=="pT2" & OutCome_IHC$N_stage=="pN0")~"IIA",
#                                     (OutCome_IHC$T_stage=="pT2" & OutCome_IHC$N_stage=="pN1")|
#                                     (OutCome_IHC$T_stage=="pT3" & OutCome_IHC$N_stage=="pN0")~"IIB",
#                                   OutCome_IHC$N_stage=="pN3"~"IIIC",
#                                   TRUE~"IIIA")
# table(OutCome_IHC$Stage_2)
# 
# 
# OutCome_IHC$Stage <- factor(OutCome_IHC$Stage,levels = c("I", "II", "III"))
# OutCome_IHC$Stage_2 <- factor(OutCome_IHC$Stage_2,levels = c("IA", "IIA", "IIB", "IIIA", "IIIC"))

OutCome_IHC$IHC <- factor(OutCome_IHC$IHC,levels = c("HR+HER2-", "HR+HER2+", "HR-HER2+","TNBC"))
colnames(OutCome_IHC)
# OutCome_IHC <- OutCome_IHC[,c("RFS_time","RFS_status","IHC_subtype","Stage")]
# OutCome_IHC <- OutCome_IHC[,c("RFS_time","RFS_status","IHC_subtype","Stage_2")]

# OutCome_IHC_1 <- OutCome_IHC[,c("RFS_time","RFS_status","T_stage","N_stage")]
# OutCome_IHC_2 <- OutCome_IHC[,c("RFS_time","RFS_status","IHC_subtype")]



#### C1

IHC_trainset <- OutCome_IHC[rownames(OutCome_IHC) %in% ID_Clin_Train,]
IHC_testset  <- OutCome_IHC[rownames(OutCome_IHC) %in% ID_Testset_All_Dimension,]
Train <- nrow(IHC_trainset)
Test  <- nrow(IHC_testset)

Trainset_candidate_IHC <- IHC_trainset
Testset_candidate_IHC  <- IHC_testset

# Multivariable Cox

multi_formulas <- as.formula(paste('Surv(RFS_time, RFS_status)~',paste(colnames(IHC_trainset)[-(1:2)],collapse = '+')))
res.cox <- coxph(multi_formulas, data = IHC_trainset)

# C-index in Testset

# C_index <- coxph(Surv(RFS_time, RFS_status)~predict(res.cox,IHC_testset), data = IHC_testset)$concordance[6]
# C_index  # CI 0.6417051

C_index <- 1-rcorr.cens(predict(res.cox,IHC_testset),Surv(IHC_testset$RFS_time,IHC_testset$RFS_status)) [[1]]
C_index  # CI 0.6417051


C_index_train <- coxph(Surv(RFS_time, RFS_status)~predict(res.cox,IHC_trainset), data = IHC_trainset)$concordance[6]
C_index_train  # CI 0.5730727


CINDEX <- function(data, indices){
  d <- data[indices,]
  cindex <- as.data.frame(1-rcorr.cens(predict(res.cox,d),Surv(d$RFS_time,d$RFS_status)))[1,1]
  return(cindex)
}

set.seed(1)

result <- boot(data = IHC_testset, statistic = CINDEX, R = 100)
result

a <- boot.ci(result, type = c("norm", "basic", "perc", "stud"))
C_index_lower95 <- a$normal[2]
C_index_upper95 <- a$normal[3]


set.seed(1)

result <- boot(data = IHC_trainset, statistic = CINDEX, R = 100)
result

b <- boot.ci(result, type = c("norm", "basic", "perc", "stud"))
C_index_lower95_train <- b$normal[2]
C_index_upper95_train <- b$normal[3]



