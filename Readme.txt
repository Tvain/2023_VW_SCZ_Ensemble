Read me file for schizophrenia (SCZ) project

Aim of the project: To achieve higher machine learning based classification precision for psychiatric disorders using ensemble learning.

Data type: Microarray Gene Expression Datasets

Packages used: Refer scripts

About the scripts: The gene expression data was processed using the scripts shared. The numbers before the title of the script indicates the sequence in which analysis was performed. Below are the details of the analysis performed using the R scripts.

1. Training and Testing data: The script includes creation of a meta-file from the gene expression datasets followed by random selection of training and testing datasets.

2. Quantile normalization: Training data was independently normalized and the quantile targets generated for each dataset/batch was used for normalizing samples of the same dataset/batch in testing data.

3. Batch correction: Different batches in train data were batch corrected using comBat. Test was batch corrected with batch corrected train data as reference. 

4. Feature selection: Differential gene expression analysis as feature selection method. Feature genes identified from each train data was used to build machine learning models.

5. SVM and 6. PAM: Models were built using different number of feature genes. Models were tested using test datasets. Models and test data predictions were saved for further analysis. 

7. Ensemble learning: Test data prediction of SVM and PAM models were ensemble using Boolean operator “AND”.  Ensemble models were evaluated based on precision.

Contact:
Vipul : vipulwagh31@gmail.com
