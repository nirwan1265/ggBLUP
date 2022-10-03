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
training_entries <- as.matrix(sample(1:93, 47))
testing_entries <- setdiff(1:93, training_entries)

Pheno_training_data <- as.matrix(Pheno[training_entries, ]) #%>% na.omit()
SNP_training_data= as.matrix(SNP_markers[training_entries,], K = NULL) 
Pheno_testing_data <- as.matrix(Pheno[testing_entries, ])# %>% na.omit
SNP_testing_data= as.matrix(SNP_markers[testing_entries,], K = NULL)


Pheno_training_data <- as.matrix(Pheno[c(1, 2, 3, 4, 5, 8, 9, 11, 15, 17, 19, 21, 23, 25, 26, 27, 29, 31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 57, 58, 59), ]) #%>% na.omit()
SNP_training_data= as.matrix(SNP_markers[c(1, 2, 3, 4, 5, 8, 9, 11, 15, 17, 19, 21, 23, 25, 26, 27, 29, 31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 57, 58, 59),], K = NULL) 
Pheno_testing_data <- as.matrix(Pheno[c(61, 62, 66, 67, 69, 70, 71, 72, 74, 75, 76, 77, 78, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92), ])# %>% na.omit
SNP_testing_data= as.matrix(SNP_markers[c(61, 62, 66, 67, 69, 70, 71, 72, 74, 75, 76, 77, 78, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92),], K = NULL)


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
cor(as.vector(Pheno_training_data), predicted_train_result, use = "complete")

quartz()
plot(Pheno_testing_data,predicted_test_result)

