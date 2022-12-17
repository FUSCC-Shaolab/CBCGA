
# Data Integration

Trainset_candidate_Met$ID  <- rownames(Trainset_candidate_Met)
Trainset_candidate_Rad$ID  <- rownames(Trainset_candidate_Rad)
Trainset_candidate_Clin$ID  <- rownames(Trainset_candidate_Clin)

Testset_candidate_Met$ID  <- rownames(Testset_candidate_Met)
Testset_candidate_Rad$ID  <- rownames(Testset_candidate_Rad)
Testset_candidate_Clin$ID  <- rownames(Testset_candidate_Clin)


tmp <- merge(Trainset_candidate_Rad,Trainset_candidate_Clin[,-c(1,2)],by = "ID")
Trainset_Met_Rad_Clin <- data.frame(merge(Trainset_candidate_Met,tmp[,-c(2,3)],by = "ID"),row.names = 1)

tmp <- merge(Testset_candidate_Rad,Testset_candidate_Clin[,-c(1,2)],by = "ID")
Testset_Met_Rad_Clin <- data.frame(merge(Testset_candidate_Met,tmp[,-c(2,3)],by = "ID"),row.names = 1)

Train <- nrow(Trainset_Met_Rad_Clin)
Test  <- nrow(Testset_Met_Rad_Clin)

# Multivariable Cox

multi_formulas <- as.formula(paste('Surv(RFS_time, RFS_status)~',paste(colnames(Trainset_Met_Rad_Clin)[-(1:2)],collapse = '+')))
res.cox <- coxph(multi_formulas, data = Trainset_Met_Rad_Clin)
res.cox.whole.BC.RFS <- res.cox
# save(res.cox.whole.BC.RFS,file = "Cox_model_for_RFS_Whole_BC.Rdata")


# C-index in Testset

C_index <- coxph(Surv(RFS_time, RFS_status)~predict(res.cox,Testset_Met_Rad_Clin), data = Testset_Met_Rad_Clin)$concordance[6]
C_index  # CI 0.6510067

C_index_train <- coxph(Surv(RFS_time, RFS_status)~predict(res.cox,Trainset_Met_Rad_Clin), data = Trainset_Met_Rad_Clin)$concordance[6]
C_index_train  # CI 0.5730727


CINDEX <- function(data, indices){
  d <- data[indices,]
  cindex <- as.data.frame(1-rcorr.cens(predict(res.cox,d),Surv(d$RFS_time,d$RFS_status)))[1,1]
  return(cindex)
}

set.seed(1)

result <- boot(data = Testset_Met_Rad_Clin, statistic = CINDEX, R = 100)
result

a <- boot.ci(result, type = c("norm", "basic", "perc", "stud"))
C_index_lower95 <- a$normal[2]
C_index_upper95 <- a$normal[3]


set.seed(1)

result <- boot(data = Trainset_Met_Rad_Clin, statistic = CINDEX, R = 100)
result

b <- boot.ci(result, type = c("norm", "basic", "perc", "stud"))
C_index_lower95_train <- b$normal[2]
C_index_upper95_train <- b$normal[3]
