setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Drosophila paper/Gemma.gwas/Data")
system("ls")
Pheno_data <- read.csv(paste0(getwd(),"/Phenotype/","dgrp.allpheno.csv"))
head(Pheno_data)
SNP_markers <- readRDS(paste0(getwd(),"/geno/","ch2L.pheno.MAF10.miss20.filtered.txt.RDS"))
SNP_markers[1:5,1:5]
x <- names(SNP_markers)
x <- unlist(strsplit(x, split='X', fixed=TRUE))[2]
library(stringr)
x <- x %>% str_replace("X","")
names(SNP_markers) <- x
SNP_markers[1:5,1:5]
names(Pheno_data)
Pheno <- Pheno_data %>% dplyr::select(BW.F_F) 
Pheno


set.seed(123)
training_entries <- as.matrix(sample(1:93, 55))
testing_entries <- setdiff(1:93, training_entries)

Pheno_training_data <- as.matrix(Pheno_data[training_entries, ]) #%>% na.omit()
SNP_training_data= as.matrix(SNP_markers[training_entries,], K = NULL) 
Pheno_testing_data <- as.matrix(Pheno_data[testing_entries, ])# %>% na.omit
SNP_testing_data= as.matrix(SNP_markers[testing_entries,], K = NULL)


Pheno_training_data <- as.matrix(Pheno[c(1,3,16,21,23,25,27,32,35,37,39,43,46,49,51,58,61,70,72,76,78,81,86,90,92), ]) #%>% na.omit()
SNP_training_data= as.matrix(SNP_markers[c(1,3,16,21,23,25,27,32,35,37,39,43,46,49,51,58,61,70,72,76,78,81,86,90,92),], K = NULL) 
Pheno_testing_data <- as.matrix(Pheno[c(2,4,17,22,24,26,34,36,38,41,45,47,50,53,60,69,71,73,77,79,84,87), ])# %>% na.omit
SNP_testing_data= as.matrix(SNP_markers[c(2,4,17,22,24,26,34,36,38,41,45,47,50,53,60,69,71,73,77,79,84,87),], K = NULL)


trained_model <- mixed.solve(y = Pheno_training_data, Z=SNP_training_data)

marker_effects <- as.matrix(trained_model$u)

#run 500 times
BLUE <- as.vector(trained_model$beta)
BLUE

predicted_train <- as.matrix(SNP_training_data) %*% marker_effects
predicted_test <- as.matrix(SNP_testing_data) %*% marker_effects
predicted_train_result <- as.vector((predicted_train[,1])+BLUE)
predicted_test_result <- as.vector((predicted_test[,1])+BLUE)

summary(as.vector(predicted_train_result))
summary(predicted_test_result)

cor(as.vector(Pheno_testing_data), predicted_test_result, use = "complete")
cor.test(as.vector(Pheno_training_data), predicted_train_result, use = "spearman")



quartz()
plot(Pheno_testing_data,predicted_test_result)

