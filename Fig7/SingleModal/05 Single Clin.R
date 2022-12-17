
CBCGA_ID <- CBCGA_Clin$PatientCode[CBCGA_Clin$CBCGA=="Yes"]
CBCGA_hospital_ID <- CBCGA_Clin$Hospital_ID[match(CBCGA_ID,CBCGA_Clin$PatientCode)]


OutCome_Clin <- data.frame(RFS_time = CBCGA_Cohort.Info$`RFS time (month)`,
                           RFS_status = CBCGA_Cohort.Info$`RFS status`,
                           T_stage = CBCGA_Cohort.Info$`Stage (pT)`,
                           N_stage = CBCGA_Cohort.Info$`Stage (pN)`,
                           IHC_subtype = CBCGA_Cohort.Info$`Clinical subtype`)
rownames(OutCome_Clin) <- CBCGA_Cohort.Info$PatientCode

table(OutCome_Clin$T_stage)
table(OutCome_Clin$N_stage)
table(OutCome_Clin$IHC_subtype)

OutCome_Clin$T_stage <- case_when(OutCome_Clin$T_stage=="pT1a"|OutCome_Clin$T_stage=="pT1b"|OutCome_Clin$T_stage=="pT1c"~"pT1",
                                  OutCome_Clin$T_stage=="pT2"~"pT2",
                                  OutCome_Clin$T_stage=="pT3"~"pT3",
                                  TRUE~"Unknown")

OutCome_Clin$T_stage <- gsub("Unknown",NA,OutCome_Clin$T_stage)


# OutCome_Clin$Stage <- case_when(OutCome_Clin$T_stage=="pT1" & OutCome_Clin$N_stage=="pN0"~"I",
#                                 (OutCome_Clin$T_stage=="pT0" & OutCome_Clin$N_stage=="pN1")|
#                                   (OutCome_Clin$T_stage=="pT1" & OutCome_Clin$N_stage=="pN1")|
#                                   (OutCome_Clin$T_stage=="pT2" & OutCome_Clin$N_stage=="pN0")|
#                                   (OutCome_Clin$T_stage=="pT2" & OutCome_Clin$N_stage=="pN1")|
#                                   (OutCome_Clin$T_stage=="pT3" & OutCome_Clin$N_stage=="pN0")~"II",
#                                 TRUE~"III")
# table(OutCome_Clin$Stage)
# 
# OutCome_Clin$Stage_2 <- case_when(OutCome_Clin$T_stage=="pT1" & OutCome_Clin$N_stage=="pN0"~"IA",
#                                   (OutCome_Clin$T_stage=="pT0" & OutCome_Clin$N_stage=="pN1")|
#                                     (OutCome_Clin$T_stage=="pT1" & OutCome_Clin$N_stage=="pN1")|
#                                     (OutCome_Clin$T_stage=="pT2" & OutCome_Clin$N_stage=="pN0")~"IIA",
#                                     (OutCome_Clin$T_stage=="pT2" & OutCome_Clin$N_stage=="pN1")|
#                                     (OutCome_Clin$T_stage=="pT3" & OutCome_Clin$N_stage=="pN0")~"IIB",
#                                   OutCome_Clin$N_stage=="pN3"~"IIIC",
#                                   TRUE~"IIIA")
# table(OutCome_Clin$Stage_2)
# 
# 
# OutCome_Clin$Stage <- factor(OutCome_Clin$Stage,levels = c("I", "II", "III"))
# OutCome_Clin$Stage_2 <- factor(OutCome_Clin$Stage_2,levels = c("IA", "IIA", "IIB", "IIIA", "IIIC"))

OutCome_Clin$IHC_subtype <- factor(OutCome_Clin$IHC_subtype,levels = c("HR+HER2-", "HR+HER2+", "HR-HER2+","TNBC"))
colnames(OutCome_Clin)
# OutCome_Clin <- OutCome_Clin[,c("RFS_time","RFS_status","IHC_subtype","Stage")]
# OutCome_Clin <- OutCome_Clin[,c("RFS_time","RFS_status","IHC_subtype","Stage_2")]

OutCome_Clin_1 <- OutCome_Clin[,c("RFS_time","RFS_status","T_stage","N_stage")]
OutCome_Clin_2 <- OutCome_Clin[,c("RFS_time","RFS_status","IHC_subtype")]



#### C1

Clin_trainset <- OutCome_Clin_1[rownames(OutCome_Clin_1) %in% ID_Clin_Train,]
Clin_testset  <- OutCome_Clin_1[rownames(OutCome_Clin_1) %in% ID_Testset_All_Dimension,]
Train <- nrow(Clin_trainset)
Test  <- nrow(Clin_testset)

Trainset_candidate_Clin <- Clin_trainset
Testset_candidate_Clin  <- Clin_testset

# Multivariable Cox

multi_formulas <- as.formula(paste('Surv(RFS_time, RFS_status)~',paste(colnames(Clin_trainset)[-(1:2)],collapse = '+')))
res.cox <- coxph(multi_formulas, data = Clin_trainset)

# C-index in Testset

C_index <- coxph(Surv(RFS_time, RFS_status)~predict(res.cox,Clin_testset), data = Clin_testset)$concordance[6]
C_index  # CI 0.6733871


C_index_train <- coxph(Surv(RFS_time, RFS_status)~predict(res.cox,Clin_trainset), data = Clin_trainset)$concordance[6]
C_index_train  # CI 0.5730727

# #### C2
# 
# Clin_trainset <- OutCome_Clin_2[rownames(OutCome_Clin) %in% ID_Clin_Train,]
# Clin_testset  <- OutCome_Clin_2[rownames(OutCome_Clin) %in% ID_Testset_All_Dimension,]
# Train <- nrow(Clin_trainset)
# Test  <- nrow(Clin_testset)
# 
# Trainset_candidate_Clin <- Clin_trainset
# Testset_candidate_Clin  <- Clin_testset
# 
# # Multivariable Cox
# 
# multi_formulas <- as.formula(paste('Surv(RFS_time, RFS_status)~',paste(colnames(Clin_trainset)[-(1:2)],collapse = '+')))
# res.cox <- coxph(multi_formulas, data = Clin_trainset)
# 
# # C-index in Testset
# 
# C_index <- coxph(Surv(RFS_time, RFS_status)~predict(res.cox,Clin_testset), data = Clin_testset)$concordance[6]
# C_index  # CI 0.6417051

CINDEX <- function(data, indices){
  d <- data[indices,]
  cindex <- as.data.frame(1-rcorr.cens(predict(res.cox,d),Surv(d$RFS_time,d$RFS_status)))[1,1]
  return(cindex)
}

set.seed(1)

result <- boot(data = Clin_testset, statistic = CINDEX, R = 100)
result

a <- boot.ci(result, type = c("norm", "basic", "perc", "stud"))
C_index_lower95 <- a$normal[2]
C_index_upper95 <- a$normal[3]


set.seed(1)

result <- boot(data = Clin_trainset, statistic = CINDEX, R = 100)
result

b <- boot.ci(result, type = c("norm", "basic", "perc", "stud"))
C_index_lower95_train <- b$normal[2]
C_index_upper95_train <- b$normal[3]
