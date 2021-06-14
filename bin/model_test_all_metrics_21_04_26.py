#!/bin/python
#04-26-2021
#model_test_all_metrics_21_04_26.py
#Alan E. Yocca

import csv
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from sklearn.inspection import permutation_importance
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import sklearn
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, auc
import statistics
import argparse
from sklearn.metrics import confusion_matrix

#not sure if all the above are necessary but hey, we got em

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

#Lets simplify this, since adding accuracy, just make three separate function
def calc_cv(ml_df = "", cv = 10, model = "", out_dict = dict()):
	#extract features and scale
	target=np.ravel(ml_df['Membership'])
	features_scale = preprocessing_ohe(ml_df = ml_df)
		
	n_splits = cv
	skf = StratifiedKFold(n_splits=n_splits)
	cv = skf.split(features_scale, target)
	
	roc_auc = []
	acc = []
	mcc = []
	cm = {"TN" : [],
		  "FP" : [],
		  "FN" : [],
		  "TP" : []}
	for train, test in cv:
		fit_model = model.fit(features_scale[train], target[train])
		probas_ = fit_model.predict_proba(features_scale[test])
		pred = fit_model.predict(features_scale[test])
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(target[test], probas_[:, 1])
		roc_auc.append(auc(fpr, tpr))
		acc.append(sklearn.metrics.accuracy_score(target[test], pred))
		mcc.append(sklearn.metrics.matthews_corrcoef(target[test], pred))
		tn, fp, fn, tp = confusion_matrix(target[test], pred).ravel()
		
		cm["TN"].append(tn)
		cm["FP"].append(fp)
		cm["FN"].append(fn)
		cm["TP"].append(tp)
	
	cm_stdev =  {"TN" : [statistics.mean(cm["TN"]), statistics.stdev(cm["TN"])],
		  		 "FP" : [statistics.mean(cm["FP"]), statistics.stdev(cm["FP"])],
		  		 "FN" : [statistics.mean(cm["FN"]), statistics.stdev(cm["FN"])],
		 		 "TP" : [statistics.mean(cm["TP"]), statistics.stdev(cm["TP"])]}
	
	loop_dict = {"ACC" : [statistics.mean(acc), statistics.stdev(acc)],
				 "AUC" : [statistics.mean(roc_auc), statistics.stdev(roc_auc)],
				 "MCC" : [statistics.mean(mcc), statistics.stdev(mcc)],
				 "CM" : cm_stdev}
	return loop_dict


#keep these for now incase we messed something up and need to reference
def calc_cv_auc(ml_df = "", cv = 10, model = ""):
	#extract features and scale
	target=np.ravel(ml_df['Membership'])
	features_scale = preprocessing_ohe(ml_df = ml_df)
		
	n_splits = cv
	skf = StratifiedKFold(n_splits=n_splits)
	cv = skf.split(features_scale, target)
	
	roc_auc = []
	for train, test in cv:
		probas_ = model.fit(features_scale[train], target[train]).predict_proba(features_scale[test])
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(target[test], probas_[:, 1])
		roc_auc.append(auc(fpr, tpr))
	
	auc_stdev = [statistics.mean(roc_auc), statistics.stdev(roc_auc)]
	return auc_stdev

def calc_cv_cm(ml_df = "", cv = 10, model = ""):
	#extract features and scale
	target=np.ravel(ml_df['Membership'])
	features_scale = preprocessing_ohe(ml_df = ml_df)
	
	n_splits = cv
	skf = StratifiedKFold(n_splits=n_splits)
	cv = skf.split(features_scale, target)
	
	cm = {"TN" : [],
		  "FP" : [],
		  "FN" : [],
		  "TP" : []}
	for train, test in cv:
		pred = model.fit(features_scale[train], target[train]).predict(features_scale[test])
		# Compute ROC curve and area the curve
		#fpr, tpr, thresholds = roc_curve(target[test], probas_[:, 1])
		#roc_auc.append(auc(fpr, tpr))
		#instead of probas_ only need pred..
		tn, fp, fn, tp = confusion_matrix(target[test], pred).ravel()
		
		cm["TN"].append(tn)
		cm["FP"].append(fp)
		cm["FN"].append(fn)
		cm["TP"].append(tp)
	
	#hmmm loop these somehow to reduce, do we want wide or tall??
	#dictionary!
	cm_stdev =  {"TN" : [statistics.mean(cm["TN"]), statistics.stdev(cm["TN"])],
		  		 "FP" : [statistics.mean(cm["FP"]), statistics.stdev(cm["FP"])],
		  		 "FN" : [statistics.mean(cm["FN"]), statistics.stdev(cm["FN"])],
		 		 "TP" : [statistics.mean(cm["TP"]), statistics.stdev(cm["TP"])]}
	return cm_stdev

def calc_cv_mcc(ml_df = "", cv = 10, model = ""):
	#extract features and scale
	target=np.ravel(ml_df['Membership'])
	features_scale = preprocessing_ohe(ml_df = ml_df)
	
	n_splits = cv
	skf = StratifiedKFold(n_splits=n_splits)
	cv = skf.split(features_scale, target)
	
	mcc = []
	for train, test in cv:
		pred = model.fit(features_scale[train], target[train]).predict(features_scale[test])
		# Compute ROC curve and area the curve
		#fpr, tpr, thresholds = roc_curve(target[test], probas_[:, 1])
		#roc_auc.append(auc(fpr, tpr))
		#instead of probas_ only need pred..
		loop_mcc = sklearn.metrics.matthews_corrcoef(target[test], pred)
		mcc.append(loop_mcc)
	
	mcc_stdev = [statistics.mean(mcc), statistics.stdev(mcc)]
	return mcc_stdev


def calc_cross(train_ml_df = "", test_ml_df = "", model = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features_scale = preprocessing_ohe(ml_df = train_ml_df)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features_scale = preprocessing_ohe(ml_df = test_ml_df)
	
	#train
	fit_model = model.fit(train_features_scale, train_target)
	probas_ = fit_model.predict_proba(test_features_scale)
	pred = fit_model.predict(test_features_scale)
	# Compute ROC curve and area the curve
	fpr, tpr, thresholds = roc_curve(test_target, probas_[:, 1])
	roc_auc = auc(fpr, tpr)
	
	acc = sklearn.metrics.accuracy_score(test_target, pred)
	mcc = sklearn.metrics.matthews_corrcoef(test_target, pred)
	
	tn, fp, fn, tp = confusion_matrix(test_target, pred).ravel()
	cm_stdev = {"TN" : [tn, 0],
		  		"FP" : [fp, 0],
		  		"FN" : [fn, 0],
		  		"TP" : [tp, 0]}
	
	loop_dict = {"ACC" : [acc, 0],
				 "AUC" : [roc_auc, 0],
				 "MCC" : [mcc, 0],
				 "CM" : cm_stdev}
	return loop_dict


def calc_auc(train_ml_df = "", test_ml_df = "", model = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features_scale = preprocessing_ohe(ml_df = train_ml_df)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features_scale = preprocessing_ohe(ml_df = test_ml_df)
	
	#train
	probas_ = model.fit(train_features_scale, train_target).predict_proba(test_features_scale)
	# Compute ROC curve and area the curve
	fpr, tpr, thresholds = roc_curve(test_target, probas_[:, 1])
	roc_auc = auc(fpr, tpr)
	
	#compare
	auc_stdev = [roc_auc, 0]
	return auc_stdev

def calc_cm(train_ml_df = "", test_ml_df = "", model = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features_scale = preprocessing_ohe(ml_df = train_ml_df)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features_scale = preprocessing_ohe(ml_df = test_ml_df)
	
	#train
	pred = model.fit(train_features_scale, train_target).predict(test_features_scale)
	# Compute confusion matrix
	tn, fp, fn, tp = confusion_matrix(test_target, pred).ravel()
	
	cm_stdev = {"TN" : [tn, 0],
		  		"FP" : [fp, 0],
		  		"FN" : [fn, 0],
		  		"TP" : [tp, 0]}
	
	return cm_stdev

def calc_mcc(train_ml_df = "", test_ml_df = "", model = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features_scale = preprocessing_ohe(ml_df = train_ml_df)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features_scale = preprocessing_ohe(ml_df = test_ml_df)
	
	#train
	pred = model.fit(train_features_scale, train_target).predict(test_features_scale)
	# Compute MCC
	mcc = sklearn.metrics.matthews_corrcoef(test_target, pred)
	
	mcc_stdev = [mcc, 0]
	return mcc_stdev


def calc_split(train_ml_df = "", test_ml_df = "", model = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features_scale = preprocessing_ohe(ml_df = train_ml_df)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features_scale = preprocessing_ohe(ml_df = test_ml_df)
	
	#ugh I need to train test split otherwise hella over training on balance v unbal on sames species
	data_train_test, data_test_test, target_train_test, target_test_test = train_test_split(test_features_scale, 
		test_target, test_size = 0.10, random_state = 10)
	
	data_train_train, data_test_train, target_train_train, target_test_train = train_test_split(train_features_scale, 
		train_target, test_size = 0.10, random_state = 10)
	
	fit_model = model.fit(data_train_train, target_train_train)
	probas_ = fit_model.predict_proba(data_test_test)
	pred = fit_model.predict(data_test_test)
	# Compute ROC curve and area the curve
	fpr, tpr, thresholds = roc_curve(target_test_test, probas_[:, 1])
	roc_auc = auc(fpr, tpr)
	
	acc = sklearn.metrics.accuracy_score(target_test_test, pred)
	mcc = sklearn.metrics.matthews_corrcoef(target_test_test, pred)
	tn, fp, fn, tp = confusion_matrix(target_test_test, pred).ravel()
	cm_stdev = {"TN" : [tn, 0],
		  		"FP" : [fp, 0],
		  		"FN" : [fn, 0],
		  		"TP" : [tp, 0]}
	
	loop_dict = {"ACC" : [acc, 0],
				 "AUC" : [roc_auc, 0],
				 "MCC" : [mcc, 0],
				 "CM" : cm_stdev}
	return loop_dict


def calc_auc_split(train_ml_df = "", test_ml_df = "", model = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features_scale = preprocessing_ohe(ml_df = train_ml_df)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features_scale = preprocessing_ohe(ml_df = test_ml_df)
	
	#ugh I need to train test split otherwise hella over training on balance v unbal on sames species
	data_train_test, data_test_test, target_train_test, target_test_test = train_test_split(test_features_scale, 
		test_target, test_size = 0.10, random_state = 10)
	
	data_train_train, data_test_train, target_train_train, target_test_train = train_test_split(train_features_scale, 
		train_target, test_size = 0.10, random_state = 10)
	
	probas_ = model.fit(data_train_train, target_train_train).predict_proba(data_test_test)
	# Compute ROC curve and area the curve
	fpr, tpr, thresholds = roc_curve(target_test_test, probas_[:, 1])
	roc_auc = auc(fpr, tpr)
	
	#compare
	auc_stdev = [roc_auc, 0]
	return auc_stdev

def calc_cm_split(train_ml_df = "", test_ml_df = "", model = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features_scale = preprocessing_ohe(ml_df = train_ml_df)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features_scale = preprocessing_ohe(ml_df = test_ml_df)
	
	#ugh I need to train test split otherwise hella over training on balance v unbal on sames species
	data_train_test, data_test_test, target_train_test, target_test_test = train_test_split(test_features_scale, 
		test_target, test_size = 0.10, random_state = 10)
	
	data_train_train, data_test_train, target_train_train, target_test_train = train_test_split(train_features_scale, 
		train_target, test_size = 0.10, random_state = 10)
	
	pred = model.fit(data_train_train, target_train_train).predict(data_test_test)
	# Compute MCC
	tn, fp, fn, tp = confusion_matrix(target_test_test, pred).ravel()
	cm_stdev = {"TN" : [tn, 0],
		  		"FP" : [fp, 0],
		  		"FN" : [fn, 0],
		  		"TP" : [tp, 0]}
	
	return cm_stdev

def calc_mcc_split(train_ml_df = "", test_ml_df = "", model = ""):
	#extract features and scale
	train_target=np.ravel(train_ml_df['Membership'])
	train_features_scale = preprocessing_ohe(ml_df = train_ml_df)
	
	test_target=np.ravel(test_ml_df['Membership'])
	test_features_scale = preprocessing_ohe(ml_df = test_ml_df)
	
	#ugh I need to train test split otherwise hella over training on balance v unbal on sames species
	data_train_test, data_test_test, target_train_test, target_test_test = train_test_split(test_features_scale, 
		test_target, test_size = 0.10, random_state = 10)
	
	data_train_train, data_test_train, target_train_train, target_test_train = train_test_split(train_features_scale, 
		train_target, test_size = 0.10, random_state = 10)
	
	pred = model.fit(data_train_train, target_train_train).predict(data_test_test)
	# Compute MCC
	mcc = sklearn.metrics.matthews_corrcoef(target_test_test, pred)
	mcc_stdev = [mcc, 0]
	return mcc_stdev

def loop_model_spec(df_dict = "", loop = 0):
	#model_list = ["SVC","GNB","RFT","SVC_wtd","RFT_wtd"]
	model_list = ["SVC","GNB","RFT"]
	species_list = ["Osat_Bal","Osat","Bdis_Bal","Bdis", "BdisA_Bal", "BdisA","BdisT","BdisT_Bal"]
				   
	combo_array = []
	for model_tag in model_list:
		for train_species in species_list:
			for test_species in species_list:
				combo_array.append([model_tag, train_species, test_species])
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
		train_species = combo[1]
		test_species = combo[2]
		print("Starting " + model_tag + " " + train_species + " " + test_species)
		
		if model_tag == "SVC":
			model = SVC(kernel='linear', C=1, probability = True)
		elif model_tag == "GNB":
			model = GaussianNB()
		else:
			model = RandomForestClassifier(n_estimators=100)

		if train_species.split("_")[0] == test_species.split("_")[0]:
			if train_species == test_species:
				loop_dict = calc_cv(ml_df = df_dict[train_species], cv = 10, 
									model = model)
							
			else:
				#train test split to avoid overtraining
				loop_dict = calc_split(train_ml_df = df_dict[train_species],
								 		   test_ml_df = df_dict[test_species],
								 		   model = model)	
				
		else:
			#cross data, can't do cv, just list acc, stdev = 0
			loop_dict = calc_cross(train_ml_df = df_dict[train_species],
								 test_ml_df = df_dict[test_species],
								 model = model)
		
		out_dict = {"ACC" : {model_tag : {train_species : {test_species : loop_dict["ACC"]}}},
					"AUC" : {model_tag : {train_species : {test_species : loop_dict["AUC"]}}},
					"MCC" : {model_tag : {train_species : {test_species : loop_dict["MCC"]}}},
					"CM"  : {model_tag : {train_species : {test_species : loop_dict["CM"]}}}}
	return out_dict
	
def output_table_auc(out_dict = dict(), filename = "", drop_header = False):
	header = ["Model","Train","Test","AUC-ROC"]
	with open(filename, "w") as out_fh:
		if not drop_header:
			for i in range(len(header) - 1):
				out_fh.write(header[i] + "\t")
			out_fh.write(header[-1] + "\n")
		for model in out_dict.keys():
			for train in out_dict[model].keys():
				for test in out_dict[model][train].keys():
					out_fh.write(model + "\t" + train + "\t" + test + "\t")
					out_fh.write(str(out_dict[model][train][test][0]) + " +/- " + 
									 str(out_dict[model][train][test][1]) + "\n")

def output_table_acc(out_dict = dict(), filename = "", drop_header = False):
	header = ["Model","Train","Test","ACC"]
	with open(filename, "w") as out_fh:
		if not drop_header:
			for i in range(len(header) - 1):
				out_fh.write(header[i] + "\t")
			out_fh.write(header[-1] + "\n")
		for model in out_dict.keys():
			for train in out_dict[model].keys():
				for test in out_dict[model][train].keys():
					out_fh.write(model + "\t" + train + "\t" + test + "\t")
					out_fh.write(str(out_dict[model][train][test][0]) + " +/- " + 
									 str(out_dict[model][train][test][1]) + "\n")

def output_table_mcc(out_dict = dict(), filename = "", drop_header = False):
	header = ["Model","Train","Test","MCC"]
	with open(filename, "w") as out_fh:
		if not drop_header:
			for i in range(len(header) - 1):
				out_fh.write(header[i] + "\t")
			out_fh.write(header[-1] + "\n")
		for model in out_dict.keys():
			for train in out_dict[model].keys():
				for test in out_dict[model][train].keys():
					out_fh.write(model + "\t" + train + "\t" + test + "\t")
					out_fh.write(str(out_dict[model][train][test][0]) + " +/- " + 
									 str(out_dict[model][train][test][1]) + "\n")

def output_table_cm(out_dict = dict(), filename = "", drop_header = False):
	header = ["Model", "Train", "Test", "TN_avg", "TN_std", "FP_avg", "FP_std",
				"FN_avg", "FN_std", "TP_avg", "TP_std"]
	#list instead of key loop to make sure order is correct
	cm_list = ["TN", "FP", "FN", "TP"]
	with open(filename, "w") as out_fh:
		if not drop_header:
			for i in range(len(header) - 1):
				out_fh.write(header[i] + "\t")
			out_fh.write(header[-1] + "\n")
		for model in out_dict.keys():
			for train in out_dict[model].keys():
				for test in out_dict[model][train].keys():
					out_fh.write(model + "\t" + train + "\t" + test + "\t")
					for i in range(0,len(cm_list)-1):
						out_fh.write(str(out_dict[model][train][test][cm_list[i]][0]) + "\t" + 
										 str(out_dict[model][train][test][cm_list[i]][1]) + "\t")
					out_fh.write(str(out_dict[model][train][test][cm_list[-1]][0]) + "\t" + 
									 str(out_dict[model][train][test][cm_list[-1]][1]) + "\n")

#want to get a nice output table of train, test, model, accuracy with +/- stdev
if __name__ == "__main__":

	wkdir = "."
	#load in osat and bdis data
	osat_ml_df = pd.read_csv(wkdir + "osat_meta_21_04_27.txt", delimiter = "\t")
	osat_ml_df = osat_ml_df.dropna()
	osat_ml_df['Membership'] = osat_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	bdis_ml_df = pd.read_csv(wkdir + "bdis_meta_21_04_30.txt", delimiter = "\t")
	bdis_ml_df = bdis_ml_df.dropna()
	bdis_ml_df['Membership'] = bdis_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	bdisA_ml_df = pd.read_csv(wkdir + "bdis_ABR2_meta_21_04_27.txt", delimiter = "\t")
	bdisA_ml_df = bdisA_ml_df.dropna()
	bdisA_ml_df['Membership'] = bdisA_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	bdisT_ml_df = pd.read_csv(wkdir + "bdis_Tek2_meta_21_04_27.txt", delimiter = "\t")
	bdisT_ml_df = bdisT_ml_df.dropna()
	bdisT_ml_df['Membership'] = bdisT_ml_df['Membership'].map({'Core': 1, 'Dispensable': 0})
	
	#balance
	osat_ml_df_bal = balance_df(ml_df = osat_ml_df)
	bdis_ml_df_bal = balance_df(ml_df = bdis_ml_df)
	bdisA_ml_df_bal = balance_df(ml_df = bdisA_ml_df)
	bdisT_ml_df_bal = balance_df(ml_df = bdisT_ml_df)
	
	df_dict = {"Osat_Bal" : osat_ml_df_bal,
			   "Osat" : osat_ml_df,
			   "Bdis_Bal" : bdis_ml_df_bal,
			   "Bdis" : bdis_ml_df,
			   "BdisA_Bal" : bdisA_ml_df_bal,
			   "BdisA" : bdisA_ml_df,
			   "BdisT_Bal" : bdisT_ml_df_bal,
			   "BdisT" : bdisT_ml_df}
	out_dict = loop_model_spec(df_dict = df_dict, loop = int(args.slurm_array_id))
	
	#output
	output_table_acc(out_dict = out_dict["ACC"], 
		filename = wkdir + '/09_ACC/osat_bdis_bdisA_bdisT_cv_acc_table_21_04_' + args.slurm_array_id + '.txt',
		drop_header = True)
	output_table_auc(out_dict = out_dict["AUC"], 
		filename = wkdir + '/03_auc_roc/osat_bdis_bdisA_bdisT_cv_auc_table_21_04_' + args.slurm_array_id + '.txt',
		drop_header = True)
	output_table_mcc(out_dict = out_dict["MCC"], 
		filename = wkdir + '/06_MCC/osat_bdis_bdisA_bdisT_cv_mcc_table_21_04_' + args.slurm_array_id + '.txt',
		drop_header = True)
	output_table_cm(out_dict = out_dict["CM"], 
		filename = wkdir + '/07_CM/osat_bdis_bdisA_bdisT_cv_cm_table_21_04_' + args.slurm_array_id + '.txt',
		drop_header = True)


