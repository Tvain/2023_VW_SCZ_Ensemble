
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(preprocessCore)

# Set directory

# Read Train1 and Test1 data
train1=readRDS("train1")

test1=readRDS("test1")

# Quantile normalize train data
# Samples from each batch needs to normalize separately in the train data
# Read pheno info 

setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\Pre-processing")

phdata=read.delim("Phenodata_Metafile_7 Batches.txt", strip.white=TRUE)

phdata.df=phdata[c(2,5,6,9)]# Keep relevant info
phdata.df=column_to_rownames(phdata,"Sample_ID") 


#######################################

## Quantile normalization 

# Each training and testing data had samples from all the 7 datasets
# To perform Q-normalization independently on each dataset. Samples from train and test data
# were split into their respective datasets.

GSE18312.Train1.raw=train1[train1.phdata$Dataset=="GSE18312",]
GSE27383.Train1.raw=train1[train1.phdata$Dataset=="GSE27383",]
GSE38481.Train1.raw=train1[train1.phdata$Dataset=="GSE38481",]
GSE38484.Train1.raw=train1[train1.phdata$Dataset=="GSE38484",]
GSE48072.Train1.raw=train1[train1.phdata$Dataset=="GSE48072",]
GSE54913.Train1.raw=train1[train1.phdata$Dataset=="GSE54913",]
Kum.Train1.raw=train1[train1.phdata$Dataset=="Kumarasinghe et al.",]

GSE18312.test1.raw=test1[test1.phdata$Dataset=="GSE18312",]
GSE27383.test1.raw=test1[test1.phdata$Dataset=="GSE27383",]
GSE38481.test1.raw=test1[test1.phdata$Dataset=="GSE38481",]
GSE38484.test1.raw=test1[test1.phdata$Dataset=="GSE38484",]
GSE48072.test1.raw=test1[test1.phdata$Dataset=="GSE48072",]
GSE54913.test1.raw=test1[test1.phdata$Dataset=="GSE54913",]
Kum.test1.raw=test1[test1.phdata$Dataset=="Kumarasinghe et al.",]


#########################################################################

# Train1 GSE18312

GSE18312.Train1.raw.df=t(GSE18312.Train1.raw)
class(GSE18312.Train1.raw.df)="numeric"
GSE18312.Train1.raw.df=log2(GSE18312.Train1.raw.df)

# Normalize
GSE18312.Train1.qnorm=normalize.quantiles(GSE18312.Train1.raw.df, copy = T)
boxplot(GSE18312.Train1.qnorm)

# Save targets
GSE18312.targets.Train1 <- normalize.quantiles.determine.target(GSE18312.Train1.qnorm,target.length=NULL,subset=NULL)

# Add col and row names from original data

dimnames(GSE18312.Train1.qnorm)=dimnames(GSE18312.Train1.raw.df)

# Normalize test data

GSE18312.test1.raw.df=t(GSE18312.test1.raw)
class(GSE18312.test1.raw.df)="numeric"
GSE18312.test1.raw.df=log2(GSE18312.test1.raw.df)

boxplot(GSE18312.test1.qnorm)

GSE18312.test1.qnorm <- normalize.quantiles.use.target(GSE18312.test1.raw.df,GSE18312.targets.Train1,copy=TRUE,subset=NULL)

# Add col and row names from original data
dimnames(GSE18312.test1.qnorm)=dimnames(GSE18312.test1.raw.df)

# (Optional) Confirm if the test data has been normalize according to the train data

#Alldata.GSE18312=cbind(GSE18312.Train1.qnorm, GSE18312.test1.qnorm)
#boxplot(Alldata.GSE18312)

###################################################################

# Train 1 GSE27383

dim(GSE27383.Train1.raw)

GSE27383.Train1.raw.df=t(GSE27383.Train1.raw)
class(GSE27383.Train1.raw.df)="numeric"
GSE27383.Train1.raw.df=log2(GSE27383.Train1.raw.df)

# Normalize
GSE27383.Train1.qnorm=normalize.quantiles(GSE27383.Train1.raw.df, copy = T)
boxplot(GSE27383.Train1.qnorm)

# Save targets
GSE27383.targets.Train1 <- normalize.quantiles.determine.target(GSE27383.Train1.qnorm,target.length=NULL,subset=NULL)

# Add col and row names from original data

dimnames(GSE27383.Train1.qnorm)=dimnames(GSE27383.Train1.raw.df)

# Normalize test data

GSE27383.test1.raw.df=t(GSE27383.test1.raw)
class(GSE27383.test1.raw.df)="numeric"
GSE27383.test1.raw.df=log2(GSE27383.test1.raw.df)

GSE27383.test1.qnorm <- normalize.quantiles.use.target(GSE27383.test1.raw.df,GSE27383.targets.Train1,copy=TRUE,subset=NULL)
boxplot(GSE27383.test1.qnorm)

# Add col and row names from original data
dimnames(GSE27383.test1.qnorm)=dimnames(GSE27383.test1.raw.df)

################################################################

# Train1 and Test1 GSE38481

dim(GSE38481.Train1.raw)

GSE38481.Train1.raw.df=t(GSE38481.Train1.raw)
class(GSE38481.Train1.raw.df)="numeric"
GSE38481.Train1.raw.df=log2(GSE38481.Train1.raw.df)

# Normalize
GSE38481.Train1.qnorm=normalize.quantiles(GSE38481.Train1.raw.df, copy = T)
boxplot(GSE38481.Train1.qnorm)

# Save targets
GSE38481.targets.Train1 <- normalize.quantiles.determine.target(GSE38481.Train1.qnorm,target.length=NULL,subset=NULL)

# Add col and row names from original data

dimnames(GSE38481.Train1.qnorm)=dimnames(GSE38481.Train1.raw.df)

# Normalize test data

GSE38481.test1.raw.df=t(GSE38481.test1.raw)
class(GSE38481.test1.raw.df)="numeric"
GSE38481.test1.raw.df=log2(GSE38481.test1.raw.df)

GSE38481.test1.qnorm <- normalize.quantiles.use.target(GSE38481.test1.raw.df,GSE38481.targets.Train1,copy=TRUE,subset=NULL)
boxplot(GSE38481.test1.qnorm)

# Add col and row names from original data
dimnames(GSE38481.test1.qnorm)=dimnames(GSE38481.test1.raw.df)

########################################################

# Train1 and Test1 GSE38484

dim(GSE38484.Train1.raw)

GSE38484.Train1.raw.df=t(GSE38484.Train1.raw)
class(GSE38484.Train1.raw.df)="numeric"
GSE38484.Train1.raw.df=log2(GSE38484.Train1.raw.df)

# Normalize
GSE38484.Train1.qnorm=normalize.quantiles(GSE38484.Train1.raw.df, copy = T)
boxplot(GSE38484.Train1.qnorm)

# Save targets
GSE38484.targets.Train1 <- normalize.quantiles.determine.target(GSE38484.Train1.qnorm,target.length=NULL,subset=NULL)

# Add col and row names from original data

dimnames(GSE38484.Train1.qnorm)=dimnames(GSE38484.Train1.raw.df)

# Normalize test data

GSE38484.test1.raw.df=t(GSE38484.test1.raw)
class(GSE38484.test1.raw.df)="numeric"
GSE38484.test1.raw.df=log2(GSE38484.test1.raw.df)

GSE38484.test1.qnorm <- normalize.quantiles.use.target(GSE38484.test1.raw.df,GSE38484.targets.Train1,copy=TRUE,subset=NULL)
boxplot(GSE38484.test1.qnorm)

# Add col and row names from original data
dimnames(GSE38484.test1.qnorm)=dimnames(GSE38484.test1.raw.df)

#######################################################

# Train1 and test1 GSE48072

dim(GSE48072.Train1.raw)

GSE48072.Train1.raw.df=t(GSE48072.Train1.raw)
class(GSE48072.Train1.raw.df)="numeric"
GSE48072.Train1.raw.df=log2(GSE48072.Train1.raw.df)

# Normalize
GSE48072.Train1.qnorm=normalize.quantiles(GSE48072.Train1.raw.df, copy = T)
boxplot(GSE48072.Train1.qnorm)

# Save targets
GSE48072.targets.Train1 <- normalize.quantiles.determine.target(GSE48072.Train1.qnorm,target.length=NULL,subset=NULL)

# Add col and row names from original data

dimnames(GSE48072.Train1.qnorm)=dimnames(GSE48072.Train1.raw.df)

# Normalize test data

GSE48072.test1.raw.df=t(GSE48072.test1.raw)
class(GSE48072.test1.raw.df)="numeric"
GSE48072.test1.raw.df=log2(GSE48072.test1.raw.df)

GSE48072.test1.qnorm <- normalize.quantiles.use.target(GSE48072.test1.raw.df,GSE48072.targets.Train1,copy=TRUE,subset=NULL)
boxplot(GSE48072.test1.qnorm)

# Add col and row names from original data
dimnames(GSE48072.test1.qnorm)=dimnames(GSE48072.test1.raw.df)

############################################################################

# Train1 and Test1 GSE54913

dim(GSE54913.Train1.raw)

GSE54913.Train1.raw.df=t(GSE54913.Train1.raw)
class(GSE54913.Train1.raw.df)="numeric"
GSE54913.Train1.raw.df=log2(GSE54913.Train1.raw.df)

# Normalize
GSE54913.Train1.qnorm=normalize.quantiles(GSE54913.Train1.raw.df, copy = T)
boxplot(GSE54913.Train1.qnorm)

# Save targets
GSE54913.targets.Train1 <- normalize.quantiles.determine.target(GSE54913.Train1.qnorm,target.length=NULL,subset=NULL)

# Add col and row names from original data

dimnames(GSE54913.Train1.qnorm)=dimnames(GSE54913.Train1.raw.df)

# Normalize test data

GSE54913.test1.raw.df=t(GSE54913.test1.raw)
class(GSE54913.test1.raw.df)="numeric"
GSE54913.test1.raw.df=log2(GSE54913.test1.raw.df)

GSE54913.test1.qnorm <- normalize.quantiles.use.target(GSE54913.test1.raw.df,GSE54913.targets.Train1,copy=TRUE,subset=NULL)
boxplot(GSE54913.test1.qnorm)

# Add col and row names from original data
dimnames(GSE54913.test1.qnorm)=dimnames(GSE54913.test1.raw.df)

########################################################################

# Train1 and Test1 Kumarasinghe et al.

dim(Kum.Train1.raw)

Kum.Train1.raw.df=t(Kum.Train1.raw)
class(Kum.Train1.raw.df)="numeric"
Kum.Train1.raw.df=log2(Kum.Train1.raw.df)

# Normalize
Kum.Train1.qnorm=normalize.quantiles(Kum.Train1.raw.df, copy = T)
boxplot(Kum.Train1.qnorm)

# Save targets
Kum.targets.Train1 <- normalize.quantiles.determine.target(Kum.Train1.qnorm,target.length=NULL,subset=NULL)

# Add col and row names from original data

dimnames(Kum.Train1.qnorm)=dimnames(Kum.Train1.raw.df)

# Normalize test data

Kum.test1.raw.df=t(Kum.test1.raw)
class(Kum.test1.raw.df)="numeric"
Kum.test1.raw.df=log2(Kum.test1.raw.df)

Kum.test1.qnorm <- normalize.quantiles.use.target(Kum.test1.raw.df,Kum.targets.Train1,copy=TRUE,subset=NULL)
boxplot(Kum.test1.qnorm)

# Add col and row names from original data
dimnames(Kum.test1.qnorm)=dimnames(Kum.test1.raw.df)

######################################################################################################

# Join all q-normalize data-sets for train and test 

train1.qnorm=cbind(GSE18312.Train1.qnorm, GSE27383.Train1.qnorm, GSE38481.Train1.qnorm, 
                   GSE38484.Train1.qnorm, GSE48072.Train1.qnorm, GSE54913.Train1.qnorm,Kum.Train1.qnorm)

test1.qnorm=cbind(GSE18312.test1.qnorm, GSE27383.test1.qnorm, GSE38481.test1.qnorm, 
                  GSE38484.test1.qnorm, GSE48072.test1.qnorm, GSE54913.test1.qnorm,Kum.test1.qnorm)

# save q-normalize train and test data

train1.qnorm.df=as.data.frame(train1.qnorm)
saveRDS(train1.qnorm.df, "train1.qnorm")

test1.qnorm.df=as.data.frame(test1.qnorm)
saveRDS(test1.qnorm.df, "test1.qnorm")

###########################################################################

# Note repeat this analysis for all the iteration of training an testing datasets.