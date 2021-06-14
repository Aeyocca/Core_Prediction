#!/bin/bash
#Alan E. Yocca
#04-09-2021
#ann_core_bdis_ref.sh
#keeping commands to create meta file here
#because command line args are many and don't want to retype this

#-u for unbuffered so prints things out to command line when running this script
wkdir="."

python -u bin/annotate_core_genes_osat_nested.py \
--pav_file data/pav_info_list_osat.txt \
--para_ka_ks data/osat_cds.osat_cds.ka.ks.txt \
--ortho_ka_ks data/osat_cds.bdis_cds_mcscan.ka.ks.txt \
--cds_file data/IRGSP-1.0_cds_2020-06-03.fasta \
--gff data/osat_transcripts_exon.gff \
--gene_type data/osat.gene_type \
--arb 5 \
--output osat_meta_21_04_27.txt
