#!/bin/bash
#Alan E. Yocca
#04-09-2021
#ann_core_bdis_ABR2.sh
#keeping commands to create meta file here
#because command line args are many and don't want to retype this

#-u for unbuffered so prints things out to command line when running this script
wkdir="."

python -u bin/annotate_core_genes_bdis.py \
--pav_file data/bdis_ABR2_pav_mat.txt \
--ortho_ka_ks data/bdis_ABR2.osat_cds.ka.ks.txt \
--para_ka_ks data/bdis_ABR2.bdis_ABR2.ka.ks.txt \
--cds_file data/BdistachyonABR2_337_v1.ABR2.1.cds_primaryTranscriptOnly.rename.fa \
--gff data/BdistachyonABR2_337_v1.ABR2.1.gene.gff3 \
--gene_type data/bdis_ABR2.gene_type \
--output bdis_ABR2_meta_21_04_09.txt \
--arb 5
