#!/bin/bash
#template
#Alan E. Yocca

#SBATCH --job-name=feat_imp_array
#SBATCH --mem=4Gb
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --output=~/Error_files/%x-%A_%a.SLURMout

#hello there

echo "Hello" 

META=${META:-`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`}

source activate core_pred

python_dir=""

cd ~/00_collab/02_core_pred/07_ms_edits/02_feat_imp

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta ${META} \
--ohe True"

echo $cmd
eval $cmd

cmd="python ${python_dir}/osat_bdis_feat_imp.py \
--meta ${META} \
--ohe True \
--balance True"

echo $cmd
eval $cmd
