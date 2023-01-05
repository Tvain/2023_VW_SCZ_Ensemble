
library(ggplot2)
library(ggfortify) 
library(limma)
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(EnhancedVolcano)

# Set directory

# Read Quantile normalized and Batch corrected train data

train1.Bcr=readRDS("train1.batch-corrected")# Matrix array

train1.Bcr.log2=log2(train1.Bcr)# NAN produced
train1.Bcr.log2[is.nan(train1.Bcr.log2)]=0 # Replace NAN with 0

# Read pheno data

setwd("F:\\PhD\\Work\\2022\\Meta-analysis final\\Pre-processing")

phdata=read.delim("Phenodata_Metafile_7 Batches.txt", strip.white=TRUE)

phdata.df=phdata[c(2,5,6,9)]# Keep relevant info
phdata.df=column_to_rownames(phdata,"Sample_ID") 

# Subset pheno-info for train 

train1.phdata=phdata.df[colnames(train1.Bcr),, drop=F]

##Conversion to factors

groups=train1.phdata$Source

f = factor(groups,levels = c("SCZ","CNT"))

##Create design Design matrix

design = model.matrix(~ 0 + f) 
colnames(design)=c("SCZ","CNT")

# Limma
data.fit = lmFit(train1.Bcr.log2,design)

##Comparison between the groups 

contrast.matrix_ = makeContrasts(SCZ-CNT,levels=design)
data.fit.con_ = contrasts.fit(data.fit,contrast.matrix_)
data.fit.eb_ = eBayes(data.fit.con_)

# Top-table for all genes

tab = topTable(data.fit.eb_,number=Inf,adjust.method="BH", confint=0.95)

topgenes = tab[tab[, "adj.P.Val"] < 0.05, ]

# Volcano plot using Enchancedvolcano function with P value


EnhancedVolcano(tab,
                lab = rownames(tab),
                x = 'logFC',
                y = 'adj.P.Val',
                ylab = "-Log10 adj.P.Value",
                pCutoff = 5e-2,
                FCcutoff = 1,
                labSize = 4,
                gridlines.major = FALSE,
                gridlines.minor = FALSE)


# Save DE file


saveRDS(tab, "train1.tab")

# Repeat the process for all the iteration of training datasets
###############################################################################################
