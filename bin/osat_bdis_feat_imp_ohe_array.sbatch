#!/bin/bash
#template
#Alan E. Yocca

#SBATCH --job-name=feat_imp_array
#SBATCH --mem=4Gb
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --output=/mnt/research/edgerpat_lab/AlanY/Error_files/%x-%A_%a.SLURMout

#hello there

#/mnt/gs18/scratch/users/yoccaala/05_ka_ks/ka_ks_19_03_08.txt

echo "Hello" 

META=${META:-`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`}

source activate core_pred

cd /mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/07_ms_edits/02_feat_imp

cmd="python /mnt/research/edgerpat_lab/AlanY/python/osat_bdis_feat_imp.py \
--meta ${META} \
--ohe True"

echo $cmd
eval $cmd

cmd="python /mnt/research/edgerpat_lab/AlanY/python/osat_bdis_feat_imp.py \
--meta ${META} \
--ohe True \
--balance True"

echo $cmd
eval $cmd
