#!/bin/bash
#Alan E. Yocca
#04-13-21
#feat_imp.sh
#wrapper for python script to calculate feature importance scores

#three meta files, balance and unbalanced
python_dir="bin/"

source activate core_prediction

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta osat_meta_21_04_27.txt \
--ohe True"
echo $cmd
eval $cmd

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta osat_meta_21_04_27.txt \
--ohe True \
--balance True"
echo $cmd
eval $cmd

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta bdis_meta_21_04_30.txt \
--ohe True"
echo $cmd
eval $cmd

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta bdis_meta_21_04_30.txt \
--ohe True \
--balance True"
echo $cmd
eval $cmd

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta bdis_ABR2_meta_21_04_09.txt \
--ohe True"
echo $cmd
eval $cmd

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta bdis_ABR2_meta_21_04_09.txt \
--ohe True \
--balance True"
echo $cmd
eval $cmd

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta bdis_Tek2_meta_21_04_27.txt \
--ohe True"
echo $cmd
eval $cmd

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta bdis_Tek2_meta_21_04_27.txt \
--ohe True \
--balance True"
echo $cmd
eval $cmd