library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(e1071)

# Ten fold random sampling using raw meta-file

# Read files GSE18312

GSE18312=readRDS("GSE18312_raw")
dim(GSE18312) #17131 21
GSE18312$rn <- rownames(GSE18312)

# Read files GSE27383
GSE27383=readRDS("GSE27383_raw")
dim(GSE27383)#21664 72
GSE27383$rn <- rownames(GSE27383)

#Read files GSE38481

GSE38481=readRDS("GSE38481_raw")
dim(GSE38481)# 12647 37
GSE38481$rn <- rownames(GSE38481)

# Read files GSE38484

GSE38484=readRDS("GSE38484_raw")
dim(GSE38484)# 17233 202
GSE38484$rn <- rownames(GSE38484)

# Read files GSE48072

GSE48072=readRDS("GSE48072_raw")
dim(GSE48072)# 15155    66
GSE48072$rn <- rownames(GSE48072)

# Read files GSE48072

GSE54913=readRDS("GSE54913_raw")
dim(GSE54913)# 13003  30
GSE54913$rn <- rownames(GSE54913)

# Read files kumarasinghe

Kum2013= readRDS("kum2013_raw")
dim(Kum2013)# 10544    21
Kum2013$rn <- rownames(Kum2013)


# Meta-file 

Meta.file_raw <- join_all(list(GSE18312,GSE27383, GSE38481, GSE38484, GSE48072, GSE54913, Kum2013), 
                          by = 'rn', type = 'inner')

# Remove "rn" (rowname) column

Meta.file_raw=column_to_rownames(Meta.file_raw, "rn")

# check dimensions
dim(Meta.file_raw)# 6775 449

# Have numbers in the df in same format (optional)
Meta.file_raw= format(Meta.file_raw, scientific = FALSE)

# Shuffle the samples (Columns)

set.seed(123)# Will ensure reproducibility

Meta.file_shf <- Meta.file_raw[,sample(ncol(Meta.file_raw))]

# Divide the df in 10 equal parts
a=t(Meta.file_shf)

nfolds = 10

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 0  # which ever fold you want to test with
train1 = a[folds != fold,]

test1 = a[folds == fold,]

# T2

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 1  # which ever fold you want to test with
train2 = a[folds != fold,]

test2 = a[folds == fold,]

# T3

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 2  # which ever fold you want to test with
train3 = a[folds != fold,]

test3 = a[folds == fold,]

# T4

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 3  # which ever fold you want to test with
train4 = a[folds != fold,]

test4 = a[folds == fold,]

# T5

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 4  # which ever fold you want to test with
train5 = a[folds != fold,]

test5 = a[folds == fold,]

# T6

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 5  # which ever fold you want to test with
train6 = a[folds != fold,]

test6 = a[folds == fold,]

# T7

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 6  # which ever fold you want to test with
train7 = a[folds != fold,]

test7 = a[folds == fold,]

# T8

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 7  # which ever fold you want to test with
train8 = a[folds != fold,]

test8 = a[folds == fold,]

# T9

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 8  # which ever fold you want to test with
train9 = a[folds != fold,]

test9 = a[folds == fold,]

# T10

set.seed(123)

folds = sample(1:nrow(a)%%nfolds) #randomize
fold = 9  # which ever fold you want to test with
train10 = a[folds != fold,]

test10 = a[folds == fold,]

# Save train and test files

saveRDS(train1, "train1")
saveRDS(test1, "test1")

