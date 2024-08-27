library(NMF)
library(tidyverse)
library(stringr)
library(dplyr)
library(gen3DNet)

### Let's define some functions ######################################################

##### Function 1: find all values that are NaN in a dataframe 
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

##### Function 2: round all values in a dataframe to the specified number of digits
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}

#######Preprocess Histone File##################################################
his_gct <- read.delim(file="/Users/tinatang/Desktop/Mollah_Lab/gen3DNet-data/LINCS_GCP_Plate23_annotated_minimized_2018-05-07_14-31-53_LVL3.gct", skip=2)
his_gct[his_gct == ""] <- NA
#get histone names
his_name <- his_gct[,"pr_gcp_histone_mark"]
his_name <- na.omit(his_name)
#get drug names
hist_drug_name <- his_gct[his_gct$id == "pert_iname",]
hist_drug_name <- t(hist_drug_name)
hist_drug_name <- na.omit(hist_drug_name)
hist_drug_name <- hist_drug_name[-1]
#subset and transpose data
histone_data <- data.frame(t(his_gct[!is.na(his_gct$pr_gene_id),!is.na(his_gct[1,])]))
histone_data <- histone_data[-1,]
#set row/column names
colnames(histone_data) <- his_name
histone_data <- cbind(X = hist_drug_name, histone_data)
names(histone_data) <- make.names(names(histone_data), unique=TRUE) # Create Unique Column Names in case there are repeat histones
#Turn NaNs to 0's and take average of replicates 
n=ncol(histone_data)
histone_data <- histone_data %>% mutate_at(c(2:n), as.numeric)
histone_final <- histone_data |> group_by(X) |> 
  summarize(across(c(1:(n-1)), .fns = mean, na.rm = TRUE))
histone_final[is.nan(histone_final)] <- 0

#######Preprocess Phosphoprotein File##############################################
phos_gct <- read.delim(file="/Users/tinatang/Desktop/Mollah_Lab/gen3DNet-data/LINCS_P100_DIA_Plate23_annotated_minimized_2018-05-02_18-13-55_LVL3.gct", skip=2)
phos_gct[phos_gct ==""] <- NA
#get phosphoprotein names
phos_name <- phos_gct[,"pr_gene_symbol"]
phos_name <- na.omit(phos_name)
phos_code <- phos_gct[,"pr_p100_modified_peptide_code"]
phos_code <- na.omit(phos_code)
phos_id <- paste(phos_name, phos_code, sep="_")
#get drug names
phos_drug_name <- phos_gct[phos_gct$id == "pert_iname",]
phos_drug_name <- t(phos_drug_name)
phos_drug_name <- na.omit(phos_drug_name)
phos_drug_name <- phos_drug_name[-1]
#subset and transpose data
phospho_data <- data.frame(t(phos_gct[!is.na(phos_gct$pr_gene_id),!is.na(phos_gct[1,])]))
phospho_data <- phospho_data[-1,]
#set row/column names
colnames(phospho_data) <- phos_id
phospho_data <- cbind(X = phos_drug_name, phospho_data)
#Turn NaNs to 0's and take average of replicates 
phospho_data <- phospho_data[, colSums(is.na(phospho_data)) <= 31] #remove any proteins with more than 1/3 missing values
n=ncol(phospho_data)
phospho_data <- phospho_data %>% mutate_at(c(2:n), as.numeric)
phospho_final <- phospho_data |> group_by(X) |> 
  summarize(across(c(1:(n-1)), .fns = mean, na.rm = TRUE))
phospho_final[is.nan(phospho_final)] <- 0

#remove any missing drugs between histone and phosphoprotein files
histone_final1 <- subset(histone_final, X %in% phospho_final$X)
phospho_final1 <- subset(phospho_final, X %in% histone_final$X)

index1=which(histone_final1$X=="DMSO")
index2=which(phospho_final1$X=="DMSO")

#set row names 
histone_final1<- as.data.frame(histone_final1)
histone_final2 <- histone_final1[,-1]
rownames(histone_final2) <- histone_final1[,1]
phospho_final1<- as.data.frame(phospho_final1)
phospho_final2 <- phospho_final1[,-1]
rownames(phospho_final2) <- phospho_final1[,1]

#calculate fold change 
m <- data.matrix(histone_final2)
histone_fc <- sweep(m[-index1,], 2, m[index1, ])
histone_fc <- as.data.frame(histone_fc)

m <- data.matrix(phospho_final2)
phospho_fc <- sweep(m[-index2,], 2, m[index2, ])
phospho_fc <- as.data.frame(phospho_fc)

histone_fc <- round_df(histone_fc, digits=5)
phospho_fc <- round_df(phospho_fc, digits=5)

#save preprocessed files
write.csv(histone_fc, "/Users/tinatang/Desktop/Mollah_Lab/gen3DNet-data/histone.csv", row.names=TRUE)
write.csv(phospho_fc, "/Users/tinatang/Desktop/Mollah_Lab/gen3DNet-data/phospho.csv", row.names=TRUE)

################################################################################

##### Use preprocessed files in Gen3DNet analysis package ######################

left <- "/Users/tinatang/Desktop/Mollah_Lab/gen3DNet-data/histone.csv"
right <- "/Users/tinatang/Desktop/Mollah_Lab/gen3DNet-data/phospho.csv"

result <- gen3DNet(
  left,
  right,
  nmf_nrun = 10,
  p_val_threshold = 0.01,
  k_picker = max_ward_kl
)
