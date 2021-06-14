#!/bin/bash
#Alan E. Yocca
#04-26-2021
#ann_core_bdis_Tek2.sh
#keeping commands to create meta file here
#because command line args are many and don't want to retype this

#-u for unbuffered so prints things out to command line when running this script
wkdir="."

python -u bin/annotate_core_genes_bdis.py \
--pav_file data/bdis_Tek2_pav_mat.txt \
--ortho_ka_ks data/bdis_Tek2.osat_cds.ka.ks.txt \
--para_ka_ks data/bdis_Tek2.bdis_Tek2.ka.ks.txt \
--cds_file data/BdistachyonTek_2_341_v1.Tek-2.1.cds_primaryTranscriptOnly.fa \
--gff data/BdistachyonTek_2_341_v1.Tek-2.1.gene.gff3 \
--gene_type data/bdis_Tek2.gene_type \
--arb 5 \
--output bdis_Tek2_meta_21_04_27.txt

