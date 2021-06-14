#!/bin/bash
#06-10-2021
#Alan E. Yocca
#rfe.sh

conda activate core_prediction

cmd="python bin/osat_bdis_rfe.py \
--meta osat_meta_21_04_27.txt --ohe True"
echo $cmd
eval $cmd

cmd="python bin/osat_bdis_rfe.py \
--meta osat_meta_21_04_27.txt --ohe True --balance True"
echo $cmd
eval $cmd

cmd="python bin/osat_bdis_rfe.py \
--meta bdis_meta_21_04_30.txt --ohe True"
echo $cmd
eval $cmd

cmd="python bin/osat_bdis_rfe.py \
--meta bdis_meta_21_04_30.txt --ohe True --balance True"
echo $cmd
eval $cmd

cmd="python bin/osat_bdis_rfe.py \
--meta bdis_ABR2_meta_21_04_09.txt --ohe True"
echo $cmd
eval $cmd

cmd="python bin/osat_bdis_rfe.py \
--meta bdis_ABR2_meta_21_04_09.txt --ohe True --balance True"
echo $cmd
eval $cmd

cmd="python bin/osat_bdis_rfe.py \
--meta bdis_Tek2_meta_21_04_27.txt --ohe True"
echo $cmd
eval $cmd

cmd="python bin/osat_bdis_rfe.py \
--meta bdis_Tek2_meta_21_04_27.txt --ohe True --balance True"
echo $cmd
eval $cmd