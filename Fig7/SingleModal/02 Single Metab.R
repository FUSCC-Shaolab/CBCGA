CBCGA_metabolomics <- rbind(CBCGA.Extended_pol,CBCGA.Extended_lip)
CBCGA_metabolomics_T <- CBCGA_metabolomics[,str_detect(colnames(CBCGA_metabolomics),"_T")]
colnames(CBCGA_metabolomics_T) <- substring(colnames(CBCGA_metabolomics_T),1,4)


CBCGA_Metabomics_tmp <- as.data.frame(t(CBCGA_metabolomics_T))
CBCGA_Metabomics_tmp <- as.data.frame(apply(CBCGA_Metabomics_tmp, 2, function(x){(x-min(x))/(max(x)-min(x))}))
OutCome_Met <- cbind(CBCGA_Cohort.Info$`RFS time (month)`[match(rownames(CBCGA_Metabomics_tmp),CBCGA_Cohort.Info$PatientCode)],
                     CBCGA_Cohort.Info$`RFS status`[match(rownames(CBCGA_Metabomics_tmp),CBCGA_Cohort.Info$PatientCode)],
                     CBCGA_Metabomics_tmp)
colnames(OutCome_Met)[c(1,2)] <- c("RFS_time", "RFS_status")

Met_trainset <- OutCome_Met[rownames(OutCome_Met) %in% ID_Met_Train,] 
Met_testset  <- OutCome_Met[rownames(OutCome_Met) %in% ID_Testset_All_Dimension,] 
Train <- nrow(Met_trainset)
Test  <- nrow(Met_testset)

# LASSO

x <- as.matrix(Met_trainset[,3:ncol(Met_trainset)])
y <- Surv(Met_trainset$RFS_time,Met_trainset$RFS_status)


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
  Feat_Name <- paste(colnames(Met_trainset)[Active_index+2], collapse = ", ") # 33
  Feat_Name
  
  # Matrix Containing Features Selected by LASSO
  
  Trainset_MultiCox <- cbind(Met_trainset[,1:2],Met_trainset[,-(1:2)][,Active_index])
  Testset_MultiCox  <- cbind(Met_testset[,1:2],Met_testset[,-(1:2)][,Active_index])
  
  colnames(Trainset_MultiCox)[3:ncol(Trainset_MultiCox)] <- colnames(Met_trainset)[Active_index+2]
  colnames(Testset_MultiCox)[3:ncol(Trainset_MultiCox)] <- colnames(Met_trainset)[Active_index+2]
  
  
  Trainset_candidate_Met <- Trainset_MultiCox
  Testset_candidate_Met  <- Testset_MultiCox
  
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


