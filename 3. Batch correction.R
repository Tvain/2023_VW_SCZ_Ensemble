library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(sva)
library(ggfortify)

# Set directory

# Read Quantile normalized train and test data

train1.qnorm=readRDS("train1.qnorm")
test1.qnorm=readRDS("test1.qnorm")

# Read pheno data

setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\Pre-processing")

phdata=read.delim("Phenodata_Metafile_7 Batches.txt", strip.white=TRUE)

phdata.df=phdata[c(2,5,6,9)]# Keep relevan info
phdata.df=column_to_rownames(phdata,"Sample_ID") 

#############################################

## Train

# Subset phenoinfo for train 

train1.phdata=phdata.df[colnames(train1.qnorm),, drop=F]

# Inverse log2 before batch correction

train1.qnorm.df=2^(train1.qnorm)
train1.qnorm.matrix=data.matrix(train1.qnorm.df)# combat needs matrix

# Combat for train
batch.train1=train1.phdata$Batch
modcombat.train1 = model.matrix(~1, data=train1.phdata)

combat.train1 = ComBat(dat=train1.qnorm.matrix, batch=batch.train1, mod=modcombat.train1, par.prior=TRUE, prior.plots=FALSE)

# PCA plot

# Raw train data
abc=prcomp(t(train1.qnorm.df), scale.=TRUE)
autoplot(abc, data=train1.phdata, colour= "Dataset")+ theme_bw()

# Batch corrected train data
abc=prcomp(t(combat.train1), scale.=TRUE)
autoplot(abc, data=train1.phdata, colour= "Dataset")+ theme_bw()

# Save batch corrected file

saveRDS(combat.train1, "train1.batch-corrected")

#############################################################

##Test

# Batch correct data with train data as reference

# Subset pheno info
test1.phdata=phdata.df[colnames(test1.qnorm),, drop=F]

# train data ph info 
train1.phdata.ref=train1.phdata
train1.phdata.ref$Batch=8

# test data ph info
test1.phdata.ref=test1.phdata

# Bind train and test phdata
test1.combat.phdata= rbind(train1.phdata.ref,test1.phdata.ref)
View(test1.combat.phdata)# batch corrected train and raw test data

# Make test data combat ready

# Inverse log2
test1.qnorm.df=2^(test1.qnorm)
test1.qnorm.matrix=data.matrix(test1.qnorm.df)# combat requires matrix

# Bind batch corrected train and raw test data
# batch corrected train and raw test data
test1.combat.data=cbind(combat.train1,test1.qnorm.matrix)


# combat
batch.test1=test1.combat.phdata$Batch
modcombat.test1 = model.matrix(~1, data=test1.combat.phdata)

combat.test1 = ComBat(dat=test1.combat.data, batch=batch.test1, mod=modcombat.test1, 
                      par.prior=TRUE, prior.plots=FALSE,  ref.batch = 8 )



# Split batch corrected test data from train data 

combat.test1.df=combat.test1[,colnames(test1.qnorm)]# test.qnorm has sample ids
dim(combat.test1.df)# Batch corrected test data

# PCA plot

# Raw test data
abc=prcomp(t(test1.qnorm.df), scale.=TRUE)
autoplot(abc, data=test1.phdata, colour= "Dataset")+ theme_bw()

# Batch corrected test data
abc=prcomp(t(combat.test1.df), scale.=TRUE)
autoplot(abc, data=test1.phdata, colour= "Dataset")+ theme_bw()

# Save files
saveRDS(combat.test1.df, "test1.batch-corrected")
#########################################################################

