library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(caret)
library(mccr)
library(stats)

# Note: prediction data for models here do not have sample IDs
#However, data was predicted from respective test data in the same order as that of 
# the sample IDs in test data

# Read predictions of test data by PAM model

setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\ML\\PAM\\R1\\Prediction")

R1.top400.PAM=read.table("top400.test1.pam.predictions.txt", strip.white = T)

# Read predictions by SVM for test data

setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\ML\\SVM\\R1\\Radial\\Prediction")

R1.top400.SVM=read.table("test1.400DEGs.radial.predictions.txt", strip.white = T)

# Read testdata (for sample IDs)

setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\ML\\Data scaling\\Batch correction\\Batch corrected files")

test1.bcr.data=readRDS("test1.batch-corrected")

# Subset 

# Read pheno data

setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\Pre-processing")

phdata=read.delim("Phenodata_Metafile_7 Batches.txt", strip.white=TRUE)

phdata.df=phdata[c(2,5,6,9)]# Keep relevant info
phdata.df=column_to_rownames(phdata.df,"Sample_ID") 

# subset pheno-info for training and testing data

test1.phdata=phdata.df[colnames(test1.bcr.data),, drop=F]

# Join predictions by ML models and true lables
R1.400.PAM=as.character(R1.top400.PAM$x)
R1.400.SVM=as.character(R1.top400.SVM$x)

R1.ensemble.df=as.data.frame(cbind( R1.400.PAM, R1.400.SVM))

# Apply boolean operator "AND"

R1.ensemble=R1.ensemble.df %>% mutate(Ensemble =
                                        case_when(R1.400.PAM == "SCZ" & R1.400.SVM=="SCZ" ~ "SCZ", 
                                                  R1.400.PAM == "SCZ" & R1.400.SVM=="CNT" ~ "CNT",
                                                  R1.400.PAM == "CNT" & R1.400.SVM=="SCZ" ~ "CNT",
                                                  R1.400.PAM == "CNT" & R1.400.SVM=="CNT" ~ "CNT"))
# add sample IDs to the ensemble data
# Save this data to observe the effect of dominant data on predictions

R1.phdata.ensemble=cbind(R1.ensemble,test1.phdata)


# Predict ensemble results
R1.test.labels=as.factor(test1.phdata$Source)
R1.ensemble.predicted=factor(as.character(R1.ensemble$Ensemble))


# Predict output
R1.ensemble.cmtable=table(true = R1.test.labels, pred = R1.ensemble.predicted)
confusionMatrix(R1.ensemble.cmtable)

# MCC

R1.ensemble.MCC=mcc(TP=R1.ensemble.cmtable[1,1],
                    FN=R1.ensemble.cmtable[1,2],
                    TN=R1.ensemble.cmtable[2,2],
                    FP=R1.ensemble.cmtable[2,1])
R1.ensemble.MCC

# Note: Repeat the process for all the iterations