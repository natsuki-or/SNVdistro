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
import time

###########argparse#############
parser = argparse.ArgumentParser(description='searches databases for SNV of a gene')
parser.add_argument('gene_name', type=str, help='gene to search')
parser.add_argument('res_loc', type=str, help='where you want the output to be saved')
diagram = parser.add_subparsers(title="diagram", dest='diagram', help='select type of diagram')


parser_2d = diagram.add_parser('2d', help='shows distribution of SNV along amino acid sequence')

parser_3d = diagram.add_parser('3d', help='shows distribution of SNV as a heatmap on pdb structure')
parser_3d.add_argument('uniprot_name', help="uniprot name of the gene is needed to search for pdb structure")
parser_3d.add_argument('--pdb_code', default='', help='you can specify pdb if you want')

args = parser.parse_args()
###############################

#make a folder to put all outputs in
res_loc = args.res_loc+"/"+args.gene_name+time.strftime("%Y%m%d_%H%M%S")
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
file_handler = logging.FileHandler(res_loc+'/'+args.gene_name+'.log')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

###############################

#locate exons and the relevent section within GRCh38
CCDS = CCDS(args.gene_name)
C_num = CCDS[0].at[0,"Chromosome"]
g_start = CCDS[0]["Start"].min()
g_stop = CCDS[0]["Stop"].max()

#extract SNVs from databases based on the location data from CCDS
logger.debug("searching ClinVar")
ClinVar = ClinVar_snv(C_num, g_start, g_stop)
logger.debug("searching dbSNP")
dbSNP = dbSNP(C_num, g_start, g_stop)
logger.debug("searching gnomAD")
gnomAD = gnomAD_snv(C_num, g_start, g_stop)

#combine the data and remove duplicates
SNV = pd.concat([ClinVar, dbSNP, gnomAD])
SNV = SNV.drop_duplicates(subset=['POS','ALT'], keep="first")
SNV['POS'] = SNV['POS'].astype(int)
SNV = SNV.sort_values(by=['POS'])
SNV = SNV.reset_index(drop=True)


#workout nuclotide residue number and original/ mutated codon/amino acid
logger.debug("working out amino acid")
SNV = mut_amino(CCDS[0],SNV,g_start, g_stop,CCDS[1])

#make a label
SNV['p_POS'] = np.ceil(SNV['nrn']/3).astype(int)
SNV['label']="p."+SNV["RefAminoAcid"].astype(str)+SNV["p_POS"].astype(str)+SNV["AltAminoAcid"].astype(str)

#extract SNV in exon
exsnv = SNV[SNV["exon_intron"]=="exon"]
exsnv = exsnv.reset_index(drop=True)
logger.debug("SNVs in exon")
logger.debug(exsnv)

#make a histogram
if args.diagram == "2d":
    diagram2d(exsnv,int(CCDS[0]["rn_len"].max()/3))
    plt.savefig(res_loc+"/"+args.gene_name+".png")


#heatmap of frequency of variation in each amino acid
if args.diagram == "3d":
    #from subprocess import *
    diagram3d(args.uniprot_name, res_loc, exsnv, args.pdb_code)


#log
files = os.listdir(res_loc)
files_file = [f for f in files if os.path.isfile(os.path.join(res_loc, f))]
ipt = vars(args).items()
s = sys.argv[0] \
+ "\n\ninputs\n" + '\n'.join(map(str, ipt)) \
+ "\n\noutputs\n" + '\n'.join(files_file) \

logger.info(s)
