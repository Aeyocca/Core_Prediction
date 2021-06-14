#!/bin/python
#Alan E. Yocca
#11-16-20
#osat_bdis_feat_imp.py
#updated 02-09-21 to include beter ka/ks estimates
#and remove amino acid proportions
#updated 02-25-21 to include ortho/para ka/ks and dup type
#sooooo if testing feature importance, why not just train with all the data?
#04-13-2021
#updated to read in cmd line args, runs svc/rft feature importances
#output file names based on input meta files

import csv
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from sklearn.inspection import permutation_importance
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import sklearn
from sklearn.model_selection import cross_val_score
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--balance", default = False, required=False, help = "if true, balance, will change output name also")
parser.add_argument("--meta", required=True, help = "meta file to use")
parser.add_argument("--ohe", default = False, required=False, help = "one hot encoding for duplication type")
args = parser.parse_args()

print("Running feature importance calculation on the file: %s" % (args.meta))
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

ml_df = ml_df.dropna()
ml_df['Membership'] = ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})

if args.balance == "True":
	ml_df = balance_df(ml_df = ml_df)
	tag = tag + "_bal"

target=np.ravel(ml_df['Membership'])

features_scale = []
#Not needed, defined above
#feature_list = []
if args.ohe == "True":
	#need to not scale ohe
	feature_list = [ "GC_Per", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks", "Para_Ka_Ks", "Para_Ka",
				   "Para_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", 
				   "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", 
				   "TA", "TC", "TG", "TT"]
	features_scale = preprocessing.scale(ml_df[feature_list])
	#append ohe_list
	#hmmm long append though..
	ohe_list = ["Singleton", "Dispersed", "Proximal", "Tandem", "WGD"]
	feature_list = feature_list + ohe_list
	features_scale = [np.append(features_scale[i],np.array(ml_df[ohe_list])[i]) for i in range(len(features_scale))]
	tag = tag + "_ohe"
else:
	features_scale = preprocessing.scale(ml_df[feature_list])


clf=RandomForestClassifier(n_estimators=100)
svc = SVC(kernel='linear', C=1)

rft_fi = clf.fit(features_scale,target).feature_importances_
svc_fi = svc.fit(features_scale,target).coef_[0]

rft_fi_dict = {"Features" : feature_list,
			   "Importance" : rft_fi}
svc_fi_dict = {"Features" : feature_list,
			   "Importance" : svc_fi}

####Convert to data frame
rft_fi_df = pd.DataFrame (rft_fi_dict, columns = ['Features','Importance'])
svc_fi_df = pd.DataFrame (svc_fi_dict, columns = ['Features','Importance'])

####OUTPUT
rft_fi_df.to_csv('rft_feat_imp_' + tag + '.csv', index=False)
svc_fi_df.to_csv('svc_feat_imp_' + tag + '.csv', index=False)




