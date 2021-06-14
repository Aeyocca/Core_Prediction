#!/bin/python
#Alan E. Yocca
#04-29-2021
#auc_curves.py

import pandas as pd
import numpy as np
from scipy import interp
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn import preprocessing
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score
from statistics import mean
from statistics import stdev
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--slurm_array_id", default = 0, required=False, help = "slurm array id, if zero will run all comps")
args = parser.parse_args()


def preprocessing_ohe(ml_df = []):
	#need to not scale the one hot encoded features, do this several times...
	#umm tempted to hard code in the ohe features and feature list here actually
	#taking it out of the loop function
	feature_list = [ "GC_Per", "Ortho_Ka_Ks", "Ortho_Ka", "Ortho_Ks", "Para_Ka_Ks", "Para_Ka",
				   "Para_Ks", "Length", "Exon_Count", "Exon_Length", "Intron_Length", 
				   "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", 
				   "TA", "TC", "TG", "TT"]
	features_scale = preprocessing.scale(ml_df[feature_list])
	#append ohe_list
	#hmmm long append though..
	ohe_list = ["Singleton", "Dispersed", "Proximal", "Tandem", "WGD"]
	feature_list = feature_list + ohe_list
	features_scale = np.append(features_scale, ml_df[ohe_list], axis = 1)
	
	return features_scale

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


def cv_auc_roc_curve(model = "", ml_df = dict(), tag = ""):
	print("Cross validating %s" % (tag))
	
	target=np.ravel(ml_df['Membership'])
	features_scale = preprocessing_ohe(ml_df = ml_df)
	
	n_splits = 10
	skf = StratifiedKFold(n_splits=n_splits)
	cv = skf.split(features_scale, target)	
	
	mean_tpr = 0.0
	mean_fpr = np.linspace(0, 1, 100)
	all_tpr = []
	all_acc = []
	all_pre = []
	all_rec = []
	
	fold = 0
	for train, test in cv:
		fold += 1
		print(fold)
		model_fit = model.fit(features_scale[train], target[train])
		probas_ = model_fit.predict_proba(features_scale[test])
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(target[test], probas_[:, 1])
		mean_tpr += np.interp(mean_fpr, fpr, tpr)
		mean_tpr[0] = 0.0
		roc_auc = auc(fpr, tpr)
		plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (fold, roc_auc))
		#Also keep track of accuracy, precision, and recall
		pred = model_fit.predict(features_scale[test])
		all_acc.append(accuracy_score(target[test], pred, normalize = True))
		all_rec.append(len([x for x in range(0,len(target[test])) if target[test][x] == 1 and pred[x] == 1])
						/ sum(target[test]))
		all_pre.append(len([x for x in range(0,len(target[test])) if target[test][x] == 1 and pred[x] == 1])
						/ sum(pred))
	
	print("Accuracy: %0.2f (+/- %0.2f)" % (mean(all_acc), stdev(all_acc) * 2))
	print("Precision: %0.2f (+/- %0.2f)" % (mean(all_pre), stdev(all_pre) * 2))
	print("Recall: %0.2f (+/- %0.2f)" % (mean(all_rec), stdev(all_rec) * 2))
	
	plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
	
	mean_tpr /= fold
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	plt.plot(mean_fpr, mean_tpr, 'k--',
	         label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)
	
	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('AUC_ROC %s' % (tag))
	plt.legend(loc="lower right")
	plt.savefig('auc_roc_%s.pdf' % (tag)) 
	plt.close('all')


def loop_model_spec(df_dict = dict(), loop = ""):
	model_list = ["SVC","GNB","RFT"]
	species_list = ["Osat_Bal","Osat","Bdis_Bal","Bdis"]
				   
	#These are all just intra-specific comps
	combo_array = []
	for model_tag in model_list:
		for species in species_list:
			combo_array.append([model_tag, species])
	#depending on slurm array, do diff combo
	
	model = ""
	if loop == 0:
		#run them all
		combo = combo_array
	else:
		#subset
		combo_array = [combo_array[loop - 1]]
	
	#initialize
	loop_dict = dict()
	
	for combo in combo_array:
		model_tag = combo[0]
		if model_tag == "SVC":
			model = SVC(kernel='linear', C=1, probability = True)
		elif model_tag == "GNB":
			model = GaussianNB()
		else:
			model = RandomForestClassifier(n_estimators=100)
	
		cv_auc_roc_curve(model = model, ml_df = df_dict[combo[1]], tag = '_'.join(combo))



if __name__ == "__main__":
	#load data
	wkdir = "."
	
	osat_ml_df = pd.read_csv(wkdir + "osat_meta_21_04_27.txt", delimiter = "\t")
	osat_ml_df = osat_ml_df.dropna()
	osat_ml_df['Membership'] = osat_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	bdis_ml_df = pd.read_csv(wkdir + "bdis_meta_21_04_30.txt", delimiter = "\t")
	bdis_ml_df = bdis_ml_df.dropna()
	bdis_ml_df['Membership'] = bdis_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	osat_ml_df_bal = balance_df(ml_df = osat_ml_df)
	bdis_ml_df_bal = balance_df(ml_df = bdis_ml_df)
	
	df_dict = {"Osat_Bal" : osat_ml_df_bal,
			   "Osat" : osat_ml_df,
			   "Bdis_Bal" : bdis_ml_df_bal,
			   "Bdis" : bdis_ml_df}
	
	loop_model_spec(df_dict = df_dict, loop = int(args.slurm_array_id))
	
