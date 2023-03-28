#!/usr/bin/env python

import time
timer_start = time.time()

import os
import sys
sys.path.append(os.getcwd())
#from control_file import *
from functions import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import logging


###########argparse#############
parser = argparse.ArgumentParser(description='searches databases for SNV of a gene')
parser.add_argument('uniprot_id', type=str, help='name of gene to search')
parser.add_argument('res_loc', type=str, help='where you want the output to be saved')
diagram = parser.add_subparsers(title="diagram", dest='diagram', help='select type of diagram')


parser_2d = diagram.add_parser('2d', help='shows distribution of SNV along amino acid sequence')

parser_3d = diagram.add_parser('3d', help='shows distribution of SNV as a heatmap on pdb structure')
parser_3d.add_argument('--user_pdb_ID', default='', help='you can specify pdb if convenient. Include the chain(e.g. 2L7B chain A => 2L7BA)')

args = parser.parse_args()
###############################

#make a folder to put all outputs in
res_loc = args.res_loc+"/"+args.uniprot_id+"_"+time.strftime("%Y%m%d_%H%M%S")
os.makedirs(res_loc)

###########logging#############
logger = logging.getLogger("SNVdistro")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')

#log on command line
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)

#log file
file_handler = logging.FileHandler(res_loc+'/'+args.uniprot_id+'.log')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

###############################

if (len(args.user_pdb_ID)>0) & (len(args.user_pdb_ID) != 5):
    logger.error("Please include the chain in your PDB ID(e.g. 2L7B chain A => 2L7BA)")
    exit()

#find CCDS ID on uniprot
CCDS_ID = Up_CCDS_ID(args.uniprot_id)
logger.debug("looking for the gene on CCDS")

#locate exons and the relevent section within GRCh38
CCDS = CCDS(CCDS_ID)
C_num = CCDS[0].at[0,"Chromosome"]
g_start = CCDS[0]["Start"].min()
g_stop = CCDS[0]["Stop"].max()

#extract SNVs from databases based on the location data from CCDS
if database_to_use["ClinVar"]:
    logger.debug("searching ClinVar")
    ClinVar = ClinVar_snv(C_num, g_start, g_stop)
else:
    logger.debug("ClinVar was not searched")
    ClinVar = pd.DataFrame(columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'CLNSIG'])

if database_to_use["dbSNP"]:
    logger.debug("searching dbSNP")
    dbSNP = dbSNP(C_num, g_start, g_stop)
else:
    logger.debug("dbSNP was not searched")
    dbSNP = pd.DataFrame(columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'CLNSIG'])

if database_to_use["gnomAD"]:
    logger.debug("searching gnomAD")
    gnomAD = gnomAD_snv(C_num, g_start, g_stop)
else:
    logger.debug("gnomAD was not searched")
    gnomAD = pd.DataFrame(columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'CLNSIG'])

if database_to_use["usersdata"]:
    logger.debug("searching usersdata")
    usersdata = usersdata_snv(C_num, g_start, g_stop)
else:
    logger.debug("usersdata was not searched")
    usersdata = pd.DataFrame(columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'CLNSIG'])


#combine the data and remove duplicates
SNV = pd.concat([ClinVar, dbSNP, gnomAD, usersdata])
SNV = SNV.drop_duplicates(subset=['POS','ALT'], keep="first")
SNV['POS'] = SNV['POS'].astype(int)
SNV = SNV.sort_values(by=['POS'])
SNV = SNV.reset_index(drop=True)

clnsig_included = [key for key, value in clnsig_to_include.items() if value]
SNV_used = SNV[SNV['CLNSIG'].isin(clnsig_included)]
if "NaN" in clnsig_included:
    SNV_used = pd.concat([SNV_used, SNV[SNV['CLNSIG'].isnull()]])

SNV_used = SNV_used.sort_values(by=['POS'])
SNV_used = SNV_used.reset_index(drop=True)

#workout nuclotide residue number and original/ mutated codon/amino acid
logger.debug("working out amino acid")
SNV_used = mut_amino(CCDS[0],SNV_used,g_start, g_stop,CCDS[1])

#make a label
SNV_used['p_POS'] = np.ceil(SNV_used['nrn']/3).astype(int)
SNV_used['label']="p."+SNV_used["RefAminoAcid"].astype(str)+SNV_used["p_POS"].astype(str)+SNV_used["AltAminoAcid"].astype(str)

#extract SNV in exon
exsnv = SNV_used[SNV_used["exon_intron"]=="exon"]
exsnv = exsnv.reset_index(drop=True)
logger.debug("SNVs in exon")
logger.debug(exsnv)
exsnv.to_csv(res_loc+"/"+args.uniprot_id+"_snv.csv", index=False)

#make a histogram
if args.diagram == "2d":
    diagram2d(exsnv,int(CCDS[0]["rn_len"].max()/3))
    plt.savefig(res_loc+"/"+args.uniprot_id+".png")


#heatmap of frequency of variation in each amino acid
if args.diagram == "3d":
    diagram3d(args.uniprot_id, res_loc, exsnv, args.user_pdb_ID)


#log
files = os.listdir(res_loc)
files_file = [f for f in files if os.path.isfile(os.path.join(res_loc, f))]
ipt = vars(args).items()
run_time = time.time() - timer_start
s = sys.argv[0] \
+ "\nrun time:" + run_time \
+ "\n\ninputs\n" + '\n'.join(map(str, ipt)) \
+ "\n\noutputs\n" + '\n'.join(files_file) \

logger.info(s)
