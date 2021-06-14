#!/bin/python
#Alan E. Yocca
#04-13-21
#osat_bdis_rfe.py
#going to try and pick best features

import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.datasets import make_classification
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--balance", default = False, required=False, help = "if true, balance, will change output name also")
parser.add_argument("--meta", required=True, help = "meta file to use")
parser.add_argument("--ohe", default = False, required=False, help = "one hot encoding for duplication type")
args = parser.parse_args()

print("Running recursive feature elimination on the file: %s" % (args.meta))
print("Balance == %s" % (str(args.balance)))
print("One Hot Encoding of dup_type == %s" % (str(args.ohe)))

def balance_df(ml_df = []):
	core_sub = ml_df[ml_df.Membership == 1]
	ncore = len(core_sub)
	disp_sub = ml_df[ml_df.Membership == 0]
	ndisp = len(disp_sub)
	
	minority_class = ""
	if ncore < ndisp:
		minority_class = 1
	else:
		minority_class = 0
	
	minority_subset = ml_df[ml_df.Membership == minority_class]
	nmin = len(minority_subset)
	majority_subset = ml_df[ml_df.Membership == abs(minority_class - 1)].sample(n=nmin, random_state=1)
	ml_df_balanced = minority_subset.append(majority_subset)
	return(ml_df_balanced)


#load in df
ml_df = pd.read_csv(args.meta, delimiter = "\t")
tag = args.meta.replace(".txt","")

feature_list = [ "GC_Per", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks", "Para_Ka_Ks", "Para_Ka",
				   "Para_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", 
				   "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", 
				   "TA", "TC", "TG", "TT", "Dup_Type"]

if args.ohe == "True":
	feature_list = [ "GC_Per", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks", "Para_Ka_Ks", "Para_Ka",
				   "Para_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", 
				   "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", 
				   "TA", "TC", "TG", "TT", "Singleton", "Dispersed", "Proximal", "Tandem", "WGD"]

ml_df = ml_df.dropna()
ml_df['Membership'] = ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})

if args.balance == "True":
	ml_df = balance_df(ml_df = ml_df)
	tag = tag + "_bal"

target=np.ravel(ml_df['Membership'])
features_scale = preprocessing.scale(ml_df[feature_list])


# Create the RFE object and compute a cross-validated score.
model = SVC(kernel='linear')
min_features_to_select = 1  # Minimum number of features to consider
rfecv = RFECV(estimator=model, step=1, cv=10,
              scoring='accuracy',
              min_features_to_select=min_features_to_select)
rfecv.fit(features_scale, target)

print("Optimal number of features SVC: %d" % rfecv.n_features_)

excluded = []
for i in range(len(list(rfecv.support_))):
	if not list(rfecv.support_)[i]:
		excluded.append(feature_list[i])

n_excluded = len(excluded)
if n_excluded > 0:
	with open(tag + '_svc_rfe_excluded.txt', 'w') as out:
		if n_excluded == 1:
			tmp = out.write(str(excluded[0]) + "\n")
		for i in range(n_excluded):
			tmp = out.write(str(excluded[i]) + "\n")

# Plot number of features VS. cross-validation scores
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Cross validation score (nb of correct classifications)")
plt.plot(range(min_features_to_select,
               len(rfecv.grid_scores_) + min_features_to_select),
         rfecv.grid_scores_)
plt.savefig(tag + '_svc_rfe.pdf')
plt.close('all')

#### Random Forest 
model = RandomForestClassifier(n_estimators=100)
min_features_to_select = 1  # Minimum number of features to consider
rfecv = RFECV(estimator=model, step=1, cv=10,
              scoring='accuracy',
              min_features_to_select=min_features_to_select)
rfecv.fit(features_scale, target)

print("Optimal number of features RFT: %d" % rfecv.n_features_)

excluded = []
for i in range(len(list(rfecv.support_))):
	if not list(rfecv.support_)[i]:
		excluded.append(feature_list[i])

n_excluded = len(excluded)
if n_excluded > 0:
	with open(tag + '_rft_rfe_excluded.txt', 'w') as out:
		if n_excluded == 1:
			tmp = out.write(str(excluded[0]) + "\n")
		for i in range(n_excluded):
			tmp = out.write(str(excluded[i]) + "\n")

# Plot number of features VS. cross-validation scores
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Cross validation score (nb of correct classifications)")
plt.plot(range(min_features_to_select,
               len(rfecv.grid_scores_) + min_features_to_select),
         rfecv.grid_scores_)
plt.savefig(tag + '_rft_rfe.pdf')
plt.close('all')

