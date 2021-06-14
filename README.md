# Core_Prediction

Did my best to include everything necessary to reproduce the results from "Machine learning approaches to identify core and dispensable genes in pangenomes"

bin/
-- Directory containing python and shell scripts

data/
-- Directory containing data to be used


## Gene feature matrix creation

We collected several features for each gene using publically available data.
These features were collected using the scripts `annotate_core_genes_osat_nested.py` and `annotate_core_genes_bdis.py`.
These scripts have several customizable command line arguments. The exact arguments used for our manuscript are listed in the shell scripts in the bin/ directory

## Model testing

We coallated model metric calculation into a single script `model_test_all_metrics_21_04_26.py`
This script calculates AUC-ROC, accuracy, MCC, and the confusion matrix for three machine learning models (support vector classifier, gaussian naive bayes, and random forest), and all pairwise combinations of training and testing data for *Oryza sativa*, and three *Brachypodium distachyon* genotypes

## AUC-ROC graph generation

We generated AUC-ROC graphs using the script `auc_curves.py`
This script uses the matplotlib python library to generate these graphs

## Recursive feature elimination

We performed recursive feature elimination using the script `osat_bdis_rfe.py`

## Confusion matrix graphs

Confusion matrix graphs were generated using the output from `annotate_core_genes_bdis.py`
Code to generate these graphs (supplementary figures 4-6) can be found in the file `bin/core_pred_cm.Rmd`

## Figure generation

Code for generating all figures can be found in the included `.Rmd` scripts
