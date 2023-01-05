library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(e1071)
library(caret)
library(mccr)
library(mltools)


# Set directory

# Read differentialy expressed genes file from train data

train1.DEGs=readRDS("train1.tab")

# Read Quantile and Batch corrected train and test data 

train1.bcr.data=readRDS("train1.batch-corrected")

test1.bcr.data=readRDS("test1.batch-corrected")

# Read pheno data

setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\Pre-processing")

phdata=read.delim("Phenodata_Metafile_7 Batches.txt", strip.white=TRUE)

phdata.df=phdata[c(2,5,6,9)]# Keep relevan info
phdata.df=column_to_rownames(phdata.df,"Sample_ID") 

# SVM for top DEGs
# Select the DEGs of interest w.r.t to the adj.P.value

abc=train1.DEGs[order(train1.DEGs$adj.P.Val, decreasing =F),]
top400DEGs=abc[1:400,]

# Subset training and testing dataset from top400 DEGs

top400.DEGs.testD=test1.bcr.data[rownames((top400DEGs)),]
top400.DEGs.testD.final=as.matrix(t(top400.DEGs.testD))

top400.DEGs.trainD=train1.bcr.data[rownames((top400DEGs)),]
top400.DEGs.trainD.final=as.matrix(t(top400.DEGs.trainD))

# subset phenoinfo

train1.400DEGs.phdata=phdata.df[colnames(top400.DEGs.trainD),, drop=F]
test1.400DEGs.phdata=phdata.df[colnames(top400.DEGs.testD),, drop=F]

# Factorize labels (SCZ/CNT) for train and test data

train1.400DEGs.flabels=factor(as.character(train1.400DEGs.phdata$Source))
test1.400DEGs.flabels=factor(as.character(test1.400DEGs.phdata$Source))

# SVM for top400

train1.400DEGs.model =best.svm(top400.DEGs.trainD.final, train1.400DEGs.flabels, 
                               cost = c(0.1,1,10), kernel = "radial", 
                               tunecontrol = tune.control(cross = 10), probability = T)

confusionMatrix(train1.400DEGs.flabels, predict(train1.400DEGs.model))

# Prediction test data samples 
test1.400DEGs.predict=predict(train1.400DEGs.model, top400.DEGs.testD.final)

test1.400DEGs.cm=table(true = test1.400DEGs.flabels, pred = test1.400DEGs.predict)
confusionMatrix(test1.400DEGs.cm)

# Save results
setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\ML\\SVM\\For kernel")

sink("R5.top400.DEGs.radial.txt")

a= "train1.top400.radial"
a

confusionMatrix(train1.400DEGs.flabels, predict(train1.400DEGs.model))

a= "test1.top400.radial"
a

confusionMatrix(test1.400DEGs.cm)

sink()
##############################################################

# Save test predicitions and trained model for ensemble learning

# Save test data predictions for ensemble 

write.table(test1.400DEGs.predict, "test1.400DEGs.predict.txt", sep="\t")

# Save trained model 

saveRDS(train1.400DEGs.model, file="train1.400DEGs.model.rds")

# Note: Repeat the process for all the iterations of training and testing datasets
# Similarly, the process can be repeated using different kernels