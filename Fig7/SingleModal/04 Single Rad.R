CBCGA_ID <- CBCGA_Clin$PatientCode[CBCGA_Clin$CBCGA=="Yes"]
CBCGA_hospital_ID <- CBCGA_Clin$Hospital_ID[match(CBCGA_ID,CBCGA_Clin$PatientCode)]
Radiomic_data_tmp <- radiomic_data
Radiomic_data_tmp <- Radiomic_data_tmp[rownames(Radiomic_data_tmp) %in% CBCGA_hospital_ID,]

CBCGA_Radiomics <- Radiomic_data_tmp
rownames(CBCGA_Radiomics) <- CBCGA_Clin$PatientCode[match(rownames(CBCGA_Radiomics),CBCGA_Clin$Hospital_ID)]
CBCGA_Radiomics <- as.data.frame(apply(CBCGA_Radiomics, 2, function(x){(x-min(x))/(max(x)-min(x))}))

OutCome_Rad <- cbind(CBCGA_Cohort.Info$`RFS time (month)`[match(rownames(CBCGA_Radiomics),CBCGA_Cohort.Info$PatientCode)],
                      CBCGA_Cohort.Info$`RFS status`[match(rownames(CBCGA_Radiomics),CBCGA_Cohort.Info$PatientCode)],
                      CBCGA_Radiomics)
colnames(OutCome_Rad)[c(1,2)] <- c("RFS_time", "RFS_status")
colnames(OutCome_Rad) <- gsub("-","_",colnames(OutCome_Rad))



Rad_trainset <- OutCome_Rad[rownames(OutCome_Rad) %in% ID_Rad_Train,] # 278
Rad_testset  <- OutCome_Rad[rownames(OutCome_Rad) %in% ID_Testset_All_Dimension,] # 146
Train <- nrow(Rad_trainset)
Test  <- nrow(Rad_testset)

# LASSO

x <- as.matrix(Rad_trainset[,3:ncol(Rad_trainset)])
y <- Surv(Rad_trainset$RFS_time,Rad_trainset$RFS_status)


set.seed(Test_Seed)

# fit <- glmnet(x, y, family = "cox", alpha = 1)
cv.fit <- cv.glmnet(x, y, family = "cox", type.measure = "C", nfolds = Test_nfolds, alpha = 1)

# LASSO Result Output
if(lambda_Para == "1se")
{
  coefficients <- coef(cv.fit,s=cv.fit$lambda.1se)
} else {
  coefficients <- coef(cv.fit,s=cv.fit$lambda.min)
}


Active_index <- which(coefficients!=0)
Active_index 


if(length(Active_index) > 0) {
  Feat_Name <- paste(colnames(Rad_trainset)[Active_index+2], collapse = ", ") # 33
  Feat_Name
  
  # Matrix Containing Features Selected by LASSO
  
  Trainset_MultiCox <- cbind(Rad_trainset[,1:2],Rad_trainset[,-(1:2)][,Active_index])
  Testset_MultiCox  <- cbind(Rad_testset[,1:2],Rad_testset[,-(1:2)][,Active_index])
  
  colnames(Trainset_MultiCox)[3:ncol(Trainset_MultiCox)] <- colnames(Rad_trainset)[Active_index+2]
  colnames(Testset_MultiCox)[3:ncol(Trainset_MultiCox)] <- colnames(Rad_trainset)[Active_index+2]
  
  
  Trainset_candidate_Rad <- Trainset_MultiCox
  Testset_candidate_Rad  <- Testset_MultiCox
  
  # Multivariable Cox
  
  multi_formulas <- as.formula(paste('Surv(RFS_time, RFS_status)~',paste(colnames(Trainset_MultiCox)[-(1:2)],collapse = '+')))
  res.cox <- coxph(multi_formulas, data = Trainset_MultiCox)
  res.cox.whole.BC.RFS <- res.cox
  # save(res.cox.whole.BC.RFS,file = "Cox_model_for_RFS_Whole_BC.Rdata")
  
  # C-index in Testset
  
  C_index <- coxph(Surv(RFS_time, RFS_status)~predict(res.cox,Testset_MultiCox), data = Testset_MultiCox)$concordance[6]
  C_index  # CI 0.5730727
  
  C_index_train <- coxph(Surv(RFS_time, RFS_status)~predict(res.cox,Trainset_MultiCox), data = Trainset_MultiCox)$concordance[6]
  C_index_train  # CI 0.5730727
  
  date()} else {
    C_index <- NA; Active_index <- NA; Feat_Name <- NA
  }


CINDEX <- function(data, indices){
  d <- data[indices,]
  cindex <- as.data.frame(1-rcorr.cens(predict(res.cox,d),Surv(d$RFS_time,d$RFS_status)))[1,1]
  return(cindex)
}

set.seed(1)

result <- boot(data = Testset_MultiCox, statistic = CINDEX, R = 100)
result

a <- boot.ci(result, type = c("norm", "basic", "perc", "stud"))
C_index_lower95 <- a$normal[2]
C_index_upper95 <- a$normal[3]


set.seed(1)

result <- boot(data = Trainset_MultiCox, statistic = CINDEX, R = 100)
result

b <- boot.ci(result, type = c("norm", "basic", "perc", "stud"))
C_index_lower95_train <- b$normal[2]
C_index_upper95_train <- b$normal[3]
