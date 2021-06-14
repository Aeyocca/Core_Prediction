#!/bin/bash
#Alan E. Yocca
#04-09-2021
#ann_core_bdis_ref.sh
#keeping commands to create meta file here
#because command line args are many and don't want to retype this

#-u for unbuffered so prints things out to command line when running this script
wkdir="."

python -u bin/annotate_core_genes_bdis.py \
--pav_file data/bdis_pav_mat_ortho.txt \
--ortho_ka_ks data/bdis_cds.osat_cds.ka.ks.txt \
--para_ka_ks data/bdis_cds.bdis_cds.ka.ks.txt \
--cds_file data/Bdistachyon_283_v2.1.cds_primaryTrancsriptOnly.fa \
--gff data/Bdistachyon_283_v2.1.gene.gff3 \
--gene_type data/bdis.gene_type \
--arb 5 \
--output bdis_meta_21_04_30.txt

