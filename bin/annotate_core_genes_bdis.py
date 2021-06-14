#!/bin/python
#Alan E. Yocca
#07-20-20
#annotate_core_genes_bdis.py

import csv
import sys
import argparse
import re

#make every file a cmd line arg so can test different accessions
#defaults will be for the reference
parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('--arb', required=False, default = 99, help='arbitrary ka/ks value')
parser.add_argument('--pav_file', required=True,
  help='SPACE SEPARATED! PAV matrix, first column is gene names, additional cols are binary calls for each accession')
parser.add_argument('--core_threshold', type = float, required=False, default = 1, help='proportion of individuals to call gene Core')
parser.add_argument('--ortho_ka_ks', required=True,
  help = 'Ka/ks file, absolute path')
parser.add_argument('--cds_file', required=True,
  help = 'cds file, absolute path')
parser.add_argument('--para_ka_ks', required=True, 
  help = 'Ka/Ks file, self vs. self, absolute path')
parser.add_argument('--gff', required=True,
parser.add_argument('--gene_type', required=True,
  help = 'gene_type file from MCScanX, absolute path')
parser.add_argument('--intron_fasta', required=False, 
  default = "",
  help = 'fasta file of first introns')

parser.add_argument('--output', required=True, help='output tsv file')

#eventually add an argument to define a smaller set of individuals

args = parser.parse_args()

def load_pav_file(filename = "", core_threshold = 1):
  pav_dict = dict()
  with open(filename) as fh:
    #next(fh)
    for line in fh:
      line_array = line.strip().split(" ")
      #number of accessions will be ncol minus first one
      #nacc=len(line_array) - 1
      npresent=sum([int(x) for x in line_array[1:]])
      prop_present=npresent/(len(line_array) - 1)
      if prop_present >= core_threshold:
      	pav_dict[line_array[0]] = {"Membership" : "Core"}
      else:
      	pav_dict[line_array[0]] = {"Membership" : "Dispensable"}
  return pav_dict
  
def add_ka_ks(filename = "", pav_dict = dict(), arb_value = 99, tag = ""):
  print("Adding arbitrary value for genes without ortholog: %s" % (arb_value))
  print("Reading ka/ks file: %s" % (filename))
  #4 is ka/ks 5 is bna gene
  ka_ks_dict = dict()
  ka_dict = dict()
  ks_dict = dict()
  with open(filename) as fh:
    for line in fh:
      next(fh)
      for line in fh:
        line_array = line.strip().split("\t")
        gene_no_trans = line_array[5].split('.')[0]
        ka_ks_dict[gene_no_trans] = line_array[4]
        ka_dict[gene_no_trans] = line_array[2]
        ks_dict[gene_no_trans] = line_array[3]
  for gene in pav_dict.keys():
    #a little difference in the string of gene names, edit here
    #gene_no_trans = gene.split('-')[0].replace("t","g")
    #pretty sure all transcript identifiers just 0.1
    try:
      pav_dict[gene][tag + "Ka_Ks"] = ka_ks_dict[gene]
      pav_dict[gene][tag + "Ka"] = ka_dict[gene]
      pav_dict[gene][tag + "Ks"] = ks_dict[gene]
    except:
      pav_dict[gene][tag + "Ka_Ks"] = arb_value
      pav_dict[gene][tag + "Ka"] = arb_value
      pav_dict[gene][tag + "Ks"] = arb_value
  
  return pav_dict

def calc_ec_el_il(cg_breaks = [], cg_strand = ""):
  #if minus strand, reverse order of list
  if cg_strand == "-":
    cg_breaks = cg_breaks[::-1]
  
  if len(cg_breaks) == 0:
    print("Failed in calc_ec_el_il")
    raise
  exon_length = 0
  intron_length = 0
  exon_count = 0
  #initialize
  previous_stop = 0
  for interval in cg_breaks:
    exon_count +=1
    exon_length += (int(interval[1]) - int(interval[0]) + 1)
    if exon_count > 1:
      intron_length += (int(interval[0]) - int(previous_stop))
    previous_stop = interval[1]
  return (exon_count, exon_length, intron_length)

def parse_gff(filename = "", pav_dict = dict()):
  loop_dict = dict()
  cg_breaks = []
  cg = ""
  transcript = ""
  cg_aed = ""
  cg_length = ""
  cg_strand = ""
  with open(filename) as fh:
    for line in fh:
      if re.match("^#",line):
        continue
      line_array = line.strip().split("\t")
      if re.match("mRNA",line_array[2]):
        #ugh, load in things, when the variables made?
        
        if len(cg_breaks) == 0:
          #print("Passing first")
          pass
        else:
          (cg_exon_count, cg_exon_length, cg_intron_length) = calc_ec_el_il(cg_breaks = cg_breaks, cg_strand = cg_strand)
          
          try:
            if loop_dict[cg][1] < cg_length:
              loop_dict[cg] = (transcript, cg_length, cg_exon_count, cg_exon_length,
            				cg_intron_length)
          except KeyError:
            loop_dict[cg] = (transcript, cg_length, cg_exon_count, cg_exon_length,
            				cg_intron_length)
        
        transcript = re.search("ID=([a-zA-Z0-9\-]*)\.",line_array[8]).group(1)
        cg = transcript
        #cg = cg.replace("t","g")
        cg_breaks = []
        cg_length = int(line_array[4]) - int(line_array[3])
        cg_strand = line_array[6]
        if cg_length == 0:
          sys.exit("Zero length gene, huh?? %s %s" % (cg, cg_length))
      
      if re.match("CDS",line_array[2]):
        cg_breaks.append([line_array[3],line_array[4]])
  
  #A Little clunky, but easiest way to replace info if found a longer isoform
  no_gene_info = 0
  longest_transcript_list = []
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      #(transcript, cg_length, cg_exon_count, cg_exon_length,cg_intron_length)
      pav_dict[gene]["Length"] = loop_dict[gene][1]
      pav_dict[gene]["Exon_Count"] = loop_dict[gene][2]
      pav_dict[gene]["Exon_Length"] = loop_dict[gene][3]
      pav_dict[gene]["Intron_Length"] = loop_dict[gene][4]
      
    except:
      no_gene_info += 1
      membership.append(pav_dict[gene]["Membership"])
      remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No pav info for %s genes in gff file, removing" % (no_gene_info))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  """
  with open("osat_longest_transcript.txt", "w") as output:
    for item in longest_transcript_list:
      output.write(item)
      output.write("\n")
  """
  
  return pav_dict

def add_go_terms(filename = "", pav_dict = dict()):
  sys.exit("You don't want to add go terms, need to make nested version later if you do")
  key_error = 0
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      #add transcript identifier so we can find it
     # gene_trans = line_array[0].replace("D2",".1D2")
      try:
        pav_dict[line_array[0]].extend([line_array[1]])
      except KeyError:
        key_error+=1
  
  print("Key errors adding go terms: %s" % (key_error))
  return pav_dict

def add_aa_table(filename = "", pav_dict = dict()):
  failed = 0
  aa_dict = dict()
  aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X"]
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      gene = line_array[0].split(".")[0]
      if gene in aa_dict:
      	continue
      for i in range(1,len(line_array)):
        try:
          aa_dict[gene].extend([line_array[i]])
        except KeyError:
          aa_dict[gene] = [line_array[i]]
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      for i in range(0,len(aa_dict[gene])):
        pav_dict[gene][aa_list[i]] = aa_dict[gene][i]
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])
      remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No aa info added for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  return pav_dict

def calc_gc_per(filename = "", pav_dict = dict()):
  header = ""
  cds_seq = dict()
  with open(filename) as fh:
    for line in fh:
      if re.match("^>",line):
        line_array = line.strip().split(" ")
        header = line_array[0].strip(">").split(".")[0]
      else:
        #also need to remove stop codons
        line = line.replace("*","")
        try:
          cds_seq[header] += line.strip()
        except KeyError:
          cds_seq[header] = line.strip()
  failed = 0
  membership = []
  for gene in pav_dict.keys():
    gc_per = ""
    try:
       gc_per = float(cds_seq[gene].count("G") + cds_seq[gene].count("C") +
      				 cds_seq[gene].count("g") + cds_seq[gene].count("c")) / float(len(cds_seq[gene]))
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])
    #ugh... going to just take the first transcript
    try:
      pav_dict[gene]["GC_Per"]
    except:
      #set if undefined
      pav_dict[gene]["GC_Per"] = gc_per
  
  print("No GC info for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  
  return pav_dict
  
def calc_dinucleotide_per(filename = "", pav_dict = dict()):
  header = ""
  cds_seq = dict()
  dinuc_pair_list = ["AA","AC","AG","AT",
                     "CA","CC","CG","CT",
                     "GA","GC","GG","GT",
                     "TA","TC","TG","TT"]
  
  with open(filename) as fh:
    for line in fh:
      if re.match("^>",line):
        line_array = line.strip().split(" ")
        header = line_array[0].strip(">").split(".")[0]
      else:
        #also need to remove stop codons
        line = line.replace("*","")
        try:
          cds_seq[header] += line.strip()
        except KeyError:
          cds_seq[header] = line.strip()
  failed = 0
  membership = []
  for gene in pav_dict.keys():
    try:
      gene_length = float(len(cds_seq[gene]))
      for dinuc_pair in dinuc_pair_list:
        #ugh... going to just take the first transcript 
        try:
          pav_dict[gene][dinuc_pair]
        except:
          #set if undefined
          #why do I have this try except here?
          pav_dict[gene][dinuc_pair] = float(cds_seq[gene].upper().count(dinuc_pair)) / gene_length
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])  
  print("No dinucleotide info for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  
  return pav_dict

def add_dup_type(filename = "", pav_dict = dict()):
  failed = 0
  dup_dict = dict()
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      dup_dict[line_array[0]] = line_array[1]
  membership = []
  remove = []
  for gene in pav_dict.keys():
    try:
      pav_dict[gene]["Dup_Type"] = dup_dict[gene]
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])
      remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No Duplication info for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  return pav_dict

def add_dup_type_one_hot_encoding(filename = "", pav_dict = dict()):
  #per reviewer request, try one hot encoding instead of single "dup_type class"
  #Type of dup	Code	Number
  #Singleton	0	
  #Dispersed	1	
  #Proximal	2	
  #Tandem	3	
  #WGD or segmental	4	
  trans_dict = {"0" : "Singleton",
  				"1" : "Dispersed",
  				"2" : "Proximal",
  				"3" : "Tandem",
  				"4" : "WGD"}
  failed = 0
  dup_dict = dict()
  with open(filename) as fh:
    for line in fh:
      line_array = line.strip().split("\t")
      dup_dict[line_array[0]] = line_array[1]
  membership = []
  remove = []
  for gene in pav_dict.keys():
    for dup_type in trans_dict.keys():
      try:
        if dup_dict[gene] == dup_type:
          pav_dict[gene][trans_dict[dup_type]] = 1
        else:
          pav_dict[gene][trans_dict[dup_type]] = 0
      except KeyError:
        failed += 1
        membership.append(pav_dict[gene]["Membership"])
        remove.append(gene)
  for gene in remove:
    del pav_dict[gene]
  print("No Duplication info for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  return pav_dict

def calc_intron_one_stats(filename = "", pav_dict = dict()):
  header = ""
  intron_seq = dict()
  with open(filename) as fh:
    for line in fh:
      if re.match("^>",line):
        line_array = line.strip().split(" ")
        header = line_array[0].strip(">").split(".")[0]
      else:
        #also need to remove stop codons
        line = line.replace("*","")
        try:
          intron_seq[header] += line.strip()
        except KeyError:
          intron_seq[header] = line.strip()
  failed = 0
  membership = []
  
  #print("Length pav_dict: %s" % (len(pav_dict.keys())))
  #print("Length intron_seq: %s" % (len(intron_seq.keys())))
  
  for gene in pav_dict.keys():
    gc_per = ""
    try:
      if len(intron_seq[gene]) == 0:
        gc_per = 0
      else:
        gc_per = float(intron_seq[gene].count("G") + intron_seq[gene].count("C") +
      				 intron_seq[gene].count("g") + intron_seq[gene].count("c")) / float(len(intron_seq[gene]))
    except KeyError:
      failed += 1
      membership.append(pav_dict[gene]["Membership"])
    #ugh... going to just take the first transcript
    try:
      pav_dict[gene]["Intron_one_GC"]
    except:
      #set if undefined
      pav_dict[gene]["Intron_one_GC"] = gc_per
      pav_dict[gene]["Intron_one_Length"] = len(intron_seq[gene])
  
  print("No First intron info for %s genes" % (failed))
  print("%s core, %s dispensable" % (len([x for x in membership if x == "Core"]),
                                     len([x for x in membership if x != "Core"])))
  
  return pav_dict

if __name__ == "__main__":

  #load feature files
  pav_dict = dict()
  #wkdir="/mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/04_brachy"
  pav_dict = load_pav_file(filename = args.pav_file, core_threshold = args.core_threshold)
  pav_dict = calc_gc_per(filename = args.cds_file, pav_dict = pav_dict)
  pav_dict = calc_dinucleotide_per(filename = args.cds_file, pav_dict = pav_dict)
  pav_dict = add_ka_ks(filename = args.ortho_ka_ks, tag = "Ortho_",
                       pav_dict = pav_dict, arb_value = args.arb)
  pav_dict = add_ka_ks(filename = args.para_ka_ks, tag = "Para_",
                       pav_dict = pav_dict, arb_value = args.arb)
  pav_dict = parse_gff(filename = args.gff, pav_dict = pav_dict)
  #pav_dict = add_exp(filename = wkdir + "/rice_exp_table.txt", pav_dict = pav_dict)              
  #pav_dict = add_aa_table(filename = wkdir + "/ref_genome/Bdis_aa_mat.txt", pav_dict = pav_dict)
  #pav_dict = add_dup_type(pav_dict = pav_dict, filename = args.gene_type)
  pav_dict = add_dup_type_one_hot_encoding(pav_dict = pav_dict, filename = args.gene_type)
  if len(args.intron_fasta) > 1:
    pav_dict = calc_intron_one_stats(pav_dict = pav_dict, filename = args.intron_fasta)
  
  with open(args.output, "w") as output:
    key_count = 1
    for key in pav_dict.keys():
      if key_count == 1:
        #add header
        key_count = 0
        output.write("Gene")
        for colname in pav_dict[key].keys():
          output.write("\t")
          output.write(colname)
        output.write("\n")
      output.write(key)
      for sub_key in pav_dict[key].keys():
        output.write("\t")
        output.write(str(pav_dict[key][sub_key]))
      output.write("\n")

  print("Finished")


"""
python /mnt/research/edgerpat_lab/AlanY/python/annotate_core_genes_bdis.py \
--ortho_ka_ks /mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/02_ka_ks/07_osat_brachy/03_codeml/bdis_cds.osat_cds.ka.ks.txt \
--para_ka_ks /mnt/research/edgerpat_lab/AlanY/00_collab/02_core_pred/02_ka_ks/07_osat_brachy/03_codeml/bdis_cds.bdis_cds.ka.ks.txt \
--output ../04_brachy/bdis_meta_21_02.txt
"""

