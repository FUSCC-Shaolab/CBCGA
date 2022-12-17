#load("CBCGA_Multiomic_Data_Partitioning.Rdata")
#load("CBCGA.Extended_MergedData_V2.6_220802.Rdata")
date()

# CBCGA_Pathway_ssGSEA <- read.csv("hallmark_ssgsea_CBCGA_Multimodal_20221118.csv")
CBCGA_Pathway_ssGSEA <- as.data.frame(apply(CBCGA_Pathway_ssGSEA, 2, function(x){(x-min(x))/(max(x)-min(x))}))
OutCome_RNA <- cbind(CBCGA_Cohort.Info$`RFS time (month)`[match(rownames(CBCGA_Pathway_ssGSEA),CBCGA_Cohort.Info$PatientCode)],
                     CBCGA_Cohort.Info$`RFS status`[match(rownames(CBCGA_Pathway_ssGSEA),CBCGA_Cohort.Info$PatientCode)],
                     CBCGA_Pathway_ssGSEA)
colnames(OutCome_RNA)[c(1,2)] <- c("RFS_time", "RFS_status")
colnames(OutCome_RNA) <- gsub("-","_",colnames(OutCome_RNA))

RNA_trainset <- OutCome_RNA[rownames(OutCome_RNA) %in% ID_RNA_Train, ] 
RNA_testset  <- OutCome_RNA[rownames(OutCome_RNA) %in% ID_Testset_All_Dimension, ] 

Train <- nrow(RNA_trainset)
Test  <- nrow(RNA_testset)

# LASSO

x <- as.matrix(RNA_trainset[,3:ncol(RNA_trainset)])
y <- Surv(RNA_trainset$RFS_time, RNA_trainset$RFS_status)

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
Feat_Name <- paste(colnames(RNA_trainset)[Active_index+2], collapse = ", ") # 33
Feat_Name

# Matrix Containing Features Selected by LASSO

Trainset_MultiCox <- cbind(RNA_trainset[,1:2],RNA_trainset[,-(1:2)][,Active_index])
Testset_MultiCox  <- cbind(RNA_testset[,1:2],RNA_testset[,-(1:2)][,Active_index])

colnames(Trainset_MultiCox)[3:ncol(Trainset_MultiCox)] <- colnames(RNA_trainset)[Active_index+2]
colnames(Testset_MultiCox)[3:ncol(Trainset_MultiCox)] <- colnames(RNA_trainset)[Active_index+2]


Trainset_candidate_RNA <- Trainset_MultiCox
Testset_candidate_RNA  <- Testset_MultiCox

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


