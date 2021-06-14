#!/bin/bash
#template
#Alan E. Yocca

#SBATCH --job-name=rfe
#SBATCH --mem=30Gb
#SBATCH --ntasks=1
#SBATCH --time=40:00:00
#SBATCH --output=/mnt/research/edgerpat_lab/AlanY/Error_files/%x-%A_%a.SLURMout

#hello there

#/mnt/gs18/scratch/users/yoccaala/05_ka_ks/ka_ks_19_03_08.txt

echo "Hello" 

source activate core_pred

cd /mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/07_ms_edits/01_rfe

META=${DIR:-`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`}

cmd="python /mnt/research/edgerpat_lab/AlanY/python/osat_bdis_rfe.py \
--meta ${META} --ohe True"
echo $cmd
eval $cmd

cmd="python /mnt/research/edgerpat_lab/AlanY/python/osat_bdis_rfe.py \
--meta ${META} --ohe True --balance True"
echo $cmd
eval $cmd
