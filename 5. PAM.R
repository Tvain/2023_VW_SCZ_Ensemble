library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(pamr)
library(caret)
library(mccr)
library(mltools)

#set directory

# Read differentially expressed genes from train data

train1.DEGs=readRDS("train1.tab")

# Read Quantile normalized and Batch corrected data

train1.bcr.data=readRDS("train1.batch-corrected")
test1.bcr.data=readRDS("test1.batch-corrected")

# Read pheno data

setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\Pre-processing")

phdata=read.delim("Phenodata_Metafile_7 Batches.txt", strip.white=TRUE)

phdata.df=phdata[c(2,5,6,9)]# Keep relevant info
phdata.df=column_to_rownames(phdata.df,"Sample_ID") 

# PAM for top DEGs

abc=train1.DEGs[order(train1.DEGs$adj.P.Val, decreasing =F),]
top400DEGs=abc[1:400,]

# Subset top400 DEGs from training and testing dataset

top400.train1=train1.bcr.data[rownames((top400DEGs)),]
top400.test1=test1.bcr.data[rownames((top400DEGs)),]

# subset pheno-info for training and testing data

top400.train1.phdata=phdata.df[colnames(top400.train1),, drop=F]
top400.test1.phdata=phdata.df[colnames(top400.test1),, drop=F]

# Factorize labels (SCZ/CNT) for train and test data

top400.train1.flabels=factor(as.character(top400.train1.phdata$Source))
top400.test1.flabels=factor(as.character(top400.test1.phdata$Source))

# Data file for PAM

SCZ_CNT=top400.train1.flabels
top400.genenames=rownames(top400.train1)
train1.sampleid=colnames(top400.train1)

top400.training.data <- list(x = top400.train1, y = SCZ_CNT, genenames = top400.genenames, 
                           geneid = top400.genenames, sampleid=train1.sampleid)

# Model with training dataset
top400.trained.model=pamr.train(top400.training.data)

# 10 Cross validation (CV). Select delta value based on CV
top400.model.cv <- pamr.cv(top400.trained.model, top400.training.data, nfold = 10)
top400.model.cv

# Delta (threshold) is chosen based on the minimum errors observed in 
# 10 fold CV
Delta = 

# Accuracy for trained model
top400.trained.model.cm=pamr.confusion(top400.trained.model, Delta, extra = F)
confusionMatrix(top400.trained.model.cm)

# Predict test samples
top400.test1.predict=pamr.predict(top400.trained.model, top400.test1, Delta)

top400.test1.predict.table=table(true = top400.test1.flabels, pred = top400.test1.predict)

confusionMatrix(top400.test1.predict.table)

MCC.top400.test1=mcc(TN=top400.test1.predict.table[1,1],
                    FP=top400.test1.predict.table[1,2],
                    TP=top400.test1.predict.table[2,2],
                    FN=top400.test1.predict.table[2,1])

MCC.top400.test1

########################
## Save data

# save trained model

saveRDS(top400.trained.model, "top400.train1.model.rds")

# Save test data predictions for ensemble

write.table(top400.test1.predict, "top400.test1.pam.predictions.txt", sep="\t")

# Save survived genes

# Save list of genes survived
setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\ML\\PAM\\R10\\Genes")

top400.train1.genes.survived=pamr.listgenes(top400.trained.model,top400.training.data, Delta)

write.table(top400.train1.genes.survived, "top400.train1.genes.survived.txt",
            sep="\t")

# Save output

sink("R10.top400.pam.txt")

a="train1.top400.pam"
a

Delta

confusionMatrix(top400.trained.model.cm)

a="test1.top400.pam"
a

confusionMatrix(top400.test1.predict.table)

MCC.top400.test1

sink()

# Note: Repeat this process for all the iterations of train and test data with
# desired number of deferentially expressed genes