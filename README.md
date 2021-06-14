# Core_Prediction

Did my best to include everything necessary to reproduce the results from "Machine learning approaches to identify core and dispensable genes in pangenomes"

bin/

-- Directory containing python and shell/submission scripts

data/

-- Directory containing data to be used

I included `environment.yml` which lists all (probably also more than is needed) python packages necessary to execute the accompanying scripts.

To get started, download and unzip this repository:

`$ unzip Core_Prediction.zip`

Now, change into this directory:

`$ cd Core_Prediction`

We need to create a conda evironment that loads in the necessary dependencies to execute our python scripts:

`$ conda create -n core_prediction --file=environment.yml`

Make sure to activate this environment before starting.

`$ conda activate core_prediction`

We are now ready to run our analyses!

## Gene feature matrix creation

We collected several features for each gene using publically available data.
These features were collected using the scripts `annotate_core_genes_osat_nested.py` and `annotate_core_genes_bdis.py`.
These scripts have several customizable command line arguments. The exact arguments used for our manuscript are listed in the shell scripts in the bin/ directory

Let's create the gene feature matricies for each genotype used in this study:

`$ bash bin/ann_core_osat.sh`

`$ bash bin/ann_core_bdis_ref.sh`

`$ bash bin/ann_core_bdis_ABR2.sh`

`$ bash bin/ann_core_bdis_Tek2.sh`

## Model testing

We coallated model metric calculation into a single script `model_test_all_metrics_21_04_26.py`
This script calculates AUC-ROC, accuracy, MCC, and the confusion matrix for three machine learning models (support vector classifier, gaussian naive bayes, and random forest), and all pairwise combinations of training and testing data for *Oryza sativa*, and three *Brachypodium distachyon* genotypes.

To run the model test, use the following command:

`$ python bin/model_test_all_metrics_21_04_26.py --slurm_array_id 0`

**WARNING** This script take a long time to run (>100 hours??). Therefore, you can alter the `--slurm_array_id` flag to run a single comparison and parallelize through a method of your choosing.
The flag is called slurm_array_id since I parallelized this script using the SLURM job management system.

## AUC-ROC graph generation

We generated AUC-ROC graphs using the script `auc_curves.py`
This script uses the matplotlib python library to generate these graphs

`$ python bin/auc_curves.py --slurm_array_id 0`

## Recursive feature elimination

We performed recursive feature elimination using the script `osat_bdis_rfe.py`
To execute this script for all datasets, run the following script:

`$ bash bin/rfe.sh`

## Feature importance

The feature importance values were calculated using the script `osat_bdis_feat_imp.py`
This script takes three command line arguments: the gene feature file (--meta_file), whether or not to balance the dataset before training (--balance True/False), and whether or not duplication type is one-hot encoded (--ohe True/False)
To execute this script for all datasets, run the following script:

`$ bash bin/feat_imp.sh`

## Confusion matrix graphs

Confusion matrix graphs were generated using the output from `model_test_all_metrics_21_04_26.py` (also included in the data/ directory)
Code to generate these graphs (supplementary figures 4-6) can be found in the file `bin/core_pred_cm.Rmd`


## Figure generation

Code for generating all figures can be found in the included `.Rmd` scripts
