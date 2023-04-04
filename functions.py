import os
import sys
sys.path.append(os.getcwd())
from control_file import *

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import re
import datetime
from urllib import request
from Bio import SeqIO
from Bio.Seq import MutableSeq
import logging
module_logger = logging.getLogger("SNVdistro.mod")



#obtains uniprot data based on id
def getuniprot(uniprot_id):
    URL = 'https://rest.uniprot.org/uniprotkb/' + uniprot_id + '.txt'
    response = request.urlopen(URL)
    content = response.read()
    response.close()
    txt = content.decode().split('\n')
    return txt


#converts uniprot id to ccds id
def Up_CCDS_ID(uniprot_id):
    Up_CCDS_ID.logger = logging.getLogger("SNVdistro.mod.Up_CCDS_ID")
    list_uniprot = getuniprot(uniprot_id)
    ccds_id_loc = [n for n, l in enumerate(list_uniprot) if l.startswith('DR   CCDS;')]
    ccds_id = re.findall('; (.*);', list_uniprot[ccds_id_loc[0]])
    if not ccds_id:
        Up_CCDS_ID.logger.error("CCDS ID not found in UniProt")
        exit()
    return ccds_id


#searches through CCDS to obtain exon locations and the relevant section of GRCh38
def CCDS(CCDS_ID):
    CCDS.logger = logging.getLogger("SNVdistro.mod.CCDS")
    CCDS_doc = pd.read_table(database_loc["CCDS"])
    CCDS_gene = CCDS_doc.loc[CCDS_doc['ccds_id'] == CCDS_ID[0]]
    if len(CCDS_gene) == 0:
        CCDS.logger.error("gene not found in CCDS")
        exit()
    CCDS_gene = CCDS_gene.reset_index(drop=True)
    C_num = CCDS_gene.iat[0,0]
    nc = CCDS_gene.iat[0,1]
    strand = CCDS_gene.iat[0,6]
    #
    exsome_loc = CCDS_gene["cds_locations"].str.extract('\[(.+)\]')
    exsome_loc = exsome_loc.iat[0,0].split(", ")
    #
    list = [[''for i in range(5)] for j in range(len(exsome_loc))]
    for k in range(len(exsome_loc)):
        list[k][0] = C_num
        list[k][1] = nc
        list[k][2] = strand
        list[k][3] = int(exsome_loc[k].split("-")[0])+1
        list[k][4] = int(exsome_loc[k].split("-")[1])+1
    df = pd.DataFrame(list, columns= ['Chromosome', 'nc_accession','Strand',"Start","Stop"])
    ###########start/stop location###############
    g_start = int(df["Start"].min())
    g_stop = int(df["Stop"].max())
    df["rn_len"] = df["Stop"].astype(int)-df["Start"].astype(int)+1
    df["rn_len"] = df["rn_len"].cumsum()
    for Chr in SeqIO.parse(database_loc["GRCh38"][str(C_num)], "fasta"):
        id_part = Chr.id
        desc_part = Chr.description
        seq = Chr.seq
    g = Chr.seq[g_start-1:g_stop]
    return df,g


#search ClinVar for snv
def ClinVar_snv(C_num, g_start, g_stop):
    ClinVar_snv.logger = logging.getLogger("SNVdistro.mod.ClinVar_snv")
    chr = []
    nt = []
    pt = []
    with open (hashtb_loc["ClinVar"],"r") as f:
        for line in f:
            tmp = line.strip().split()
            chr.append(tmp[0])
            nt.append(int(tmp[1]))
            pt.append(int(tmp[2]))
    total = len(chr)
    t_chr = str(C_num)
    st    = int(g_start)
    fn    = int(g_stop)
    lines_t=[]
    for i in range(total):
        if t_chr == chr[i]:
            if nt[i]+hashtb_DIV["ClinVar"] > st:
                pointer = pt[i]
                break
    with open(database_loc["ClinVar"],"r") as f:
        f.seek(pointer)
        for line in f:
            if not str(line.split()[0]) == t_chr:
                break
            if int(line.split()[1]) > fn:
                break
            if int(line.split()[1]) >= st:
                #print(line[0:50])
                lines_t.append(re.sub('\s+',' ',line))
    v = pd.DataFrame(lines_t)
    if v.empty:
        ClinVar.logger.warning("No varient found on ClinVar")
        cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'CLNSIG']
        g_snv = pd.DataFrame(columns=cols)
        return g_snv
    v = v[0].str.split(expand=True)
    v.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    v["CLNSIG"] = v['INFO'].str.extract('(?:.*CLNSIG=)(.*?)(?:;)')
    v["CLNVC"] = v['INFO'].str.extract('(?:.*CLNVC=)(.*?)(?:;)')
    #delete "INFO" column
    v.drop('INFO', inplace=True, axis=1)
    #extract single nucleotide variants
    g_snv = v[v["CLNVC"]=="single_nucleotide_variant"]
    g_snv.drop('CLNVC', inplace=True, axis=1)
    g_snv = g_snv.reset_index(drop=True)
    clnsig_included = [key for key, value in clnsig_to_include.items() if value]
    g_snv = g_snv[g_snv['CLNSIG'].isin(clnsig_included)]
    if "NaN" in clnsig_included:
        g_snv = pd.concat([g_snv, g_snv[g_snv['CLNSIG'].isnull()]])
    return g_snv


# for dbSNP entry with multiple clnsig
def mean_clnsig(string):
    numbers = string.replace('/', '|').split('|')
    numbers = ["0" if x == '' or x == '.' else x for x in numbers]
    num_in_str =[clin_sig_dict.get(item,item)  for item in numbers]
    scaled_num =[clnsig_scale.get(item,item)  for item in num_in_str]
    scaled_mean=sum(scaled_num) / len(scaled_num)
    if scaled_mean != 0:
        return round(math.sqrt(abs(scaled_mean))*(scaled_mean/abs(scaled_mean)),3)
    else:
        return 0.000


#search dbSNP for snv
def dbSNP(C_num, g_start, g_stop):
    dbSNP.logger = logging.getLogger("SNVdistro.mod.dbSNP")
    chr_an = {"NC_000001.11":1, "NC_000002.12":2, "NC_000003.12":3, "NC_000004.12":4, "NC_000005.10":5, \
            "NC_000006.12":6, "NC_000007.14":7, "NC_000008.11":8, "NC_000009.12":9, "NC_000010.11":10, \
            "NC_000011.10":11, "NC_000012.12":12, "NC_000013.11":13, "NC_000014.9":14, "NC_000015.10":15, \
            "NC_000016.10":16, "NC_000017.11":17, "NC_000018.10":18, "NC_000019.10":19, "NC_000020.11":20, \
            "NC_000021.9":21, "NC_000022.11":22, "NC_000023.11":"X", "NC_000024.10":"Y", "NC_012920.1":"MT"}
    clin_sig_dict = {"0":"Uncertain significance",\
                ".": "Uncertain significance",\
                "1": "not provided",\
                "2": "Benign",\
                "3": "Likely benign",\
                "4": "Likely pathogenic",\
                "5": "Pathogenic",\
                "6": "drug response",\
                "8": "confers sensitivity",\
                "9": "risk-factor",\
                "10": "association",\
                "11": "protective",\
                "12": "conflict",\
                "13": "affects",\
                "255": "other"}
    chr = []
    nt = []
    pt = []
    with open (hashtb_loc["dbSNP"],"r") as f:
        for line in f:
            tmp = line.strip().split()
            chr.append(tmp[0])
            nt.append(int(tmp[1]))
            pt.append(int(tmp[2]))
    total = len(chr)
    t_chr = str([i for i, j in chr_an.items() if j == int(C_num)][0])
    st    = int(g_start)
    fn    = int(g_stop)
    lines_t=[]
    for k in range(total):
        if t_chr == chr[k]:
            if nt[k]+hashtb_DIV["dbSNP"] > st:
                pointer = pt[k]
                break
    with open(database_loc["dbSNP"],"r") as f:
        f.seek(pointer)
        for line in f:
            if not str(line.split()[0]) == t_chr:
                break
            if int(line.split()[1]) > fn:
                break
            if int(line.split()[1]) >= st:
                #print(line[0:50])
                lines_t.append(re.sub('\s+',' ',line))
    v = pd.DataFrame(lines_t)
    if v.empty:
        dbSNP.logger.warning("No varient found on dbSNP")
        cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'CLNSIG']
        g_snv = pd.DataFrame(columns=cols)
        return g_snv
    v = v[0].str.split(expand=True)
    v.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    v["CLNVC"] = v['INFO'].str.extract('.+[;]VC[=](.*?)[;]')
    v["CLNSIG"] = v['INFO'].str.extract('(?:.*CLNSIG=.,)(.*?)(?:;)')
    v.drop('INFO', inplace=True, axis=1)
    v["CHROM"] = chr_an[v["CHROM"][0]]
    g_snv = v[v["CLNVC"]=="SNV"]
    g_snv.drop('CLNVC', inplace=True, axis=1)
    g_snv = g_snv.reset_index(drop=True)
    g_snv["CLNSIG"] = g_snv["CLNSIG"].fillna('')
    for ls in range(len(g_snv)):
        g_snv["ALT"][ls]= g_snv["ALT"][ls].split(",")
        g_snv["CLNSIG"][ls]= g_snv["CLNSIG"][ls].split(",")
        if len(g_snv["ALT"][ls]) > len(g_snv["CLNSIG"][ls]):
            for i in range(len(g_snv["ALT"][ls])-len(g_snv["CLNSIG"][ls])):
                g_snv["CLNSIG"][ls].append('')
    g_snv = g_snv.explode(["ALT","CLNSIG"])
    g_snv = g_snv.reset_index(drop=True)
    for mc in range(len(g_snv)):
        g_snv["CLNSIG"][mc] = mean_clnsig(g_snv["CLNSIG"][mc])
    g_snv.loc[(g_snv['CLNSIG'] >= clnsig_threshold["lower_limit"]) & (g_snv['CLNSIG'] <= clnsig_threshold["upper_limit"])]
    return g_snv


#search gnomAD for snv
def gnomAD_snv(C_num, g_start, g_stop):
    gnomAD_snv.logger = logging.getLogger("SNVdistro.mod.gnomAD_snv")
    nt = []
    pt = []
    with open (hashtb_loc["gnomAD"]+"gnomAD.chr"+str(C_num)+".hashtb","r") as f:
        for line in f:
            tmp = line.strip().split()
            nt.append(int(tmp[0]))
            pt.append(int(tmp[1]))
    total = len(nt)
    st    = int(g_start)
    fn    = int(g_stop)
    lines_t=[]
    for i in range(total):
        if nt[i]+hashtb_DIV["gnomAD"] > st:
            pointer = pt[i]
            break
    with open(database_loc["gnomAD"]+"gnomad.genomes.v3.1.2.sites.chr"+str(C_num)+".vcf","r") as f:
        f.seek(pointer)
        for line in f:
            if int(line.split()[1]) > fn:
                break
            if int(line.split()[1]) >= st:
                #print(line[0:50])
                lines_t.append(re.sub('\s+',' ',line))
    v = pd.DataFrame(lines_t)
    if v.empty:
        gnomAD_snv.logger.warning("No varient found on gnomAD")
        cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'CLNSIG']
        g_snv = pd.DataFrame(columns=cols)
        return g_snv
    v = v[0].str.split(expand=True)
    v.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    v.drop('INFO', inplace=True, axis=1)
    v["CLNVC"] = ""
    for j in range(len(v)):
        if (len(v["REF"][j])==1) & (len(v["ALT"][j])==1):
            v["CLNVC"][j] = "SNV"
    #extract single nucleotide variants
    g_snv = v[v["CLNVC"]=="SNV"]
    g_snv.drop('CLNVC', inplace=True, axis=1)
    g_snv = g_snv.reset_index(drop=True)
    g_snv["CHROM"] = C_num
    return g_snv


#search user's data for snv
def usersdata_snv(C_num, g_start, g_stop):
    usersdata_snv.logger = logging.getLogger("SNVdistro.mod.usersdata_snv")
    t_chr = str(C_num)
    st    = int(g_start)
    fn    = int(g_stop)
    lines_t=[]
    with open(database_loc["usersdata"],"r") as f:
        lines = [l for l in f if not l.startswith('#')]
        for line in lines:
            if (str(line.split()[0]) == t_chr) & (int(line.split()[1]) >= st):
                if (int(line.split()[1]) > fn):
                    break
                else:
                    lines_t.append(re.sub('\s+',' ',line))
    v = pd.DataFrame(lines_t)
    v = v[0].str.split(expand=True)
    v.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    v["CLNSIG"] = v['INFO'].str.extract('(?:.*CLNSIG=)(.*?)(?:;)')
    v["CLNVC"] = v['INFO'].str.extract('(?:.*CLNVC=)(.*?)(?:;)')
    v.drop('INFO', inplace=True, axis=1) #delete "INFO" column
    g_snv = v[v["CLNVC"]=="single_nucleotide_variant"] #extract single nucleotide variants
    g_snv.drop('CLNVC', inplace=True, axis=1)
    g_snv = g_snv.reset_index(drop=True)
    return g_snv


#checks if snv causes amino acid change
def mut_amino(df,g_snv,g_start, g_stop,g):
    g_snv["exon_intron"] = '' #8
    g_snv['nrn'] = '' #9
    #
    g_start = int(g_start)
    g_stop = int(g_stop)
    df[["Start", "Stop", "rn_len"]] = df[["Start", "Stop", "rn_len"]].apply(pd.to_numeric)
    for l in range(len(df)):
        g_snv.loc[g_snv["POS"].between(df["Start"][l], df["Stop"][l]), 'exon_intron'] = 'exon'
        if l == 0:
            g_snv.loc[g_snv["POS"].between(df["Start"][l], df["Stop"][l]), 'nrn'] = g_snv["POS"]-df["Start"][l]+1
        else:
            g_snv.loc[g_snv["POS"].between(df["Start"][l], df["Stop"][l]), 'nrn'] = g_snv["POS"]-df["Start"][l]+1+df["rn_len"][l-1]
    g_snv["exon_intron"]= g_snv["exon_intron"].replace("","intron")
    g_snv["nrn"]= g_snv["nrn"].replace("",0)
    #######
    if df["Strand"][0] == "-":
        g_snv["nrn"] = df["rn_len"].max()-g_snv[g_snv["exon_intron"]=="exon"]["nrn"]+1
        g_snv["nrn"]= g_snv["nrn"].replace(np.NaN,0).astype(int)
        comp = {"A":"T","T":"A","G":"C","C":"G"}
        g_snv.replace({"REF": comp}, inplace=True)
        g_snv.replace({"ALT": comp}, inplace=True)
    ########
    modulo = divmod(g_snv["nrn"].astype(int), 3)
    g_snv['RefCodon'] = ''      #col_num 10
    g_snv['AltCodon'] = ''      #col_num 11
    g_snv['RefAminoAcid'] = ''  #col_num 12
    g_snv['AltAminoAcid'] = ''  #col_num 13
    #
    if df["Strand"][0] == "+":
        for o in range(len(g_snv)):
            if g_snv["nrn"][o] > 0:
                try:
                    x = modulo[1][o]
                    y = g_snv["POS"].astype(int)-g_start+1
                    z = int(y[o]-3/2*x**2+7/2*x-2)
                    g_snv.iloc[o,10] = g[z-1:z+2]
                    g_snv.iloc[o,11] = MutableSeq(g_snv.iloc[o,10])
                    g_snv.iloc[o,11][int(modulo[1][o]-1)] = g_snv.iloc[o,4]
                    g_snv.iloc[o,12] = g_snv.iloc[o,10].translate()
                    g_snv.iloc[o,13] = g_snv.iloc[o,11].translate()
                except Exception:
                    pass
    ###
    if df["Strand"][0] == "-":
        for o in range(len(g_snv)):
            if g_snv["nrn"][o] > 0:
                try:
                    x = modulo[1][o]
                    y = len(g)-(g_snv["POS"].astype(int)-g_start)
                    z = int(y[o]-3/2*x**2+7/2*x-2)
                    g_snv.iloc[o,10] = g.reverse_complement()[z-1:z+2]
                    g_snv.iloc[o,11] = MutableSeq(g_snv.iloc[o,10])
                    g_snv.iloc[o,11][int(modulo[1][o]-1)] = g_snv.iloc[o,4]
                    g_snv.iloc[o,12] = g_snv.iloc[o,10].translate()
                    g_snv.iloc[o,13] = g_snv.iloc[o,11].translate()
                except Exception:
                    pass
    return g_snv


#blast search on pdb
def pdbblast(seq):
    pdbblast.logger = logging.getLogger("SNVdistro.mod.pdbblast\n")
    import requests
    s = {
      "query": {
      "type": "terminal",
      "service": "sequence",
      "parameters": {
      "evalue_cutoff": 1,
      "identity_cutoff": 0.9,
      "sequence_type": "protein",
      "value": seq
      }
      },
      "request_options": {"scoring_strategy": "sequence"},
      "return_type": "polymer_instance"
      }
    r = requests.post('https://search.rcsb.org/rcsbsearch/v2/query', json=s).json()
    ID = list()
    evalue = list()
    for result in r['result_set']:
        ID.append([result][0]['identifier'])
        evalue.append(float([result][0]["score"]))
    pdbblast.logger.info(pd.DataFrame(list(zip(ID, evalue)),columns =['ID', 'evalue']))
    return ID,evalue




#obtains amino acid sequence baced on mmcif

RESDAT = "./residue.dat"
ATOM = ['_atom_site.group_PDB',
        '_atom_site.id',
        '_atom_site.type_symbol',
        '_atom_site.label_atom_id',
        '_atom_site.label_alt_id',
        '_atom_site.label_comp_id',
        '_atom_site.label_asym_id',
        '_atom_site.label_entity_id',
        '_atom_site.label_seq_id',
        '_atom_site.pdbx_PDB_ins_code',
        '_atom_site.Cartn_x',
        '_atom_site.Cartn_y',
        '_atom_site.Cartn_z',
        '_atom_site.occupancy',
        '_atom_site.B_iso_or_equiv',
        '_atom_site.pdbx_formal_charge',
        '_atom_site.auth_seq_id',
        '_atom_site.auth_comp_id',
        '_atom_site.auth_asym_id',
        '_atom_site.auth_atom_id',
        '_atom_site.pdbx_PDB_model_num']

def readmmcif(file):
    pdb = list()
    with open(file,"r") as f:
        for line in f:
            pdb.append(line.replace('\n',''))
            if "_entry.id" in line:
                pdbid = line.strip().split()[1]

    return pdbid, pdb

def ext_atom(pdb):
    atom = list()
    for i in range(len(pdb)):
        if ('ATOM'in pdb[i] or 'HETATM' in pdb[i]) and \
            (pdb[i].split()[ATOM.index('_atom_site.label_alt_id')] == '.' or \
             pdb[i].split()[ATOM.index('_atom_site.label_alt_id')] == 'A' or \
             pdb[i].split()[ATOM.index('_atom_site.label_alt_id')] == '1') and \
            pdb[i].split()[ATOM.index('_atom_site.pdbx_PDB_model_num')] == '1':
            atom.append(pdb[i].split())
    return(atom)

def newnum(atom):
    sernum = list()
    num = 1
    sernum.append(num)
    for i in range(1,len(atom)):
        if  atom[i][ATOM.index('_atom_site.label_asym_id')] != \
            atom[i-1][ATOM.index('_atom_site.label_asym_id')]:
            num = 1
        elif  atom[i][ATOM.index('_atom_site.label_seq_id')] != \
              atom[i-1][ATOM.index('_atom_site.label_seq_id')]:
            num += 1
        sernum.append(num)
    return sernum

def het2atom(atom,sernum):
    total = len(atom)
    st = -1
    alt = 0
    for i in range(total):
        if st == -1 and \
           atom[i][ATOM.index('_atom_site.group_PDB')] == 'HETATM':
            st = i-1
        if st != -1 and \
           atom[i][ATOM.index('_atom_site.group_PDB')] == 'ATOM':
            fn = i
            if atom[st][ATOM.index('_atom_site.label_asym_id')] == \
               atom[fn][ATOM.index('_atom_site.label_asym_id')] and \
               sernum[st] < sernum[fn] and \
               int(atom[st][ATOM.index('_atom_site.label_seq_id')]) <= \
               int(atom[fn][ATOM.index('_atom_site.label_seq_id')]):
                for j in range(st+1,fn):
                    atom[j][ATOM.index('_atom_site.group_PDB')] = 'ATOM'
                alt = 1
                print("REMARK   HETATM -> ATOM [%s]" % atom[st+1][ATOM.index('_atom_site.label_comp_id')],file=sys.stderr)
            st = -1
    for i in range(total):
        if atom[i][ATOM.index('_atom_site.group_PDB')] == 'ATOM' and \
           sernum[i] != int(atom[i][ATOM.index('_atom_site.label_seq_id')]):
            print("REMARK   RENUMBERING",file=sys.stderr)
            break
    return atom

def replace_res(atom,sernum):
    total = len(atom)
    for i in range(total):
        atom[i][ATOM.index('_atom_site.label_seq_id')] = sernum[i]
    return atom

def conv(atom):
    ot1 = list()
    ot3 = list()
    with open(RESDAT,"r") as f:
        for line in f:
            tmp = line.strip().split()
            ot3.append(tmp[0])
            ot1.append(tmp[1])

    chain = list()
    ch    = atom[0][ATOM.index('_atom_site.label_asym_id')]
    chain.append(atom[0][ATOM.index('_atom_site.label_asym_id')])

    seq   = list()
    sq    = ''

    total = len(atom)
    for i in range(total):
        if i != 0 and ch != atom[i][ATOM.index('_atom_site.label_asym_id')]:
            if len(sq) != 0:
                seq.append(sq)
            else:
                chain.pop()
            sq = ''

            ch = atom[i][ATOM.index('_atom_site.label_asym_id')]
            chain.append(ch)

        if atom[i][ATOM.index('_atom_site.group_PDB')] == "ATOM" and \
           atom[i][ATOM.index('_atom_site.label_atom_id')] == "CA":
            sq += ot1[ot3.index(atom[i][ATOM.index('_atom_site.label_comp_id')])]
    if len(sq) != 0:
        seq.append(sq)
    else:
        chain.pop()
    return(chain,seq)

def atom2seq(id):
    pdbid,pdb = readmmcif(id)
    atom = ext_atom(pdb)
    sernum = newnum(atom)
    atom = het2atom(atom,sernum)
    atom = replace_res(atom,sernum)
    chain,seq = conv(atom)
    return dict(zip(chain, seq))




#make a histogram of snv
def diagram2d(exsnv,b):
    import matplotlib.pyplot as plt
    #
    a = exsnv["p_POS"].to_numpy()
    cm = plt.cm.get_cmap('hot_r')
    n, bins, patches = plt.hist(a, bins=b, color='green')
    col = (n-n.min())/(n.max()-n.min())
    for c, p in zip(col, patches):
        plt.setp(p, 'facecolor', cm(c))
    return print("diagram made")



def rgb(minimum, maximum, value):
    rgb.logger = logging.getLogger("SNVdistro.mod.rgb")
    minimum, maximum = float(minimum), float(maximum)
    try:
        ratio = 2 * (value-minimum) / (maximum - minimum)
    except ZeroDivisionError:
        rgb.logger.exception('Tried to divide by zero(i.e. min and max number of snv was the same). Check freq_z')
    else:
        b = int(max(0, 256*(1 - ratio)))
        r = int(max(0, 256*(ratio - 1)))
        g = 256 - b - r
        b = int(min(256, 256*(1 - ratio)+256))
        r = int(min(256, 256*(ratio - 1)+256))
        return [r, g, b]


#make a 3d heatmap of snv
def diagram3d(uniprot_id,res_loc,exsnv,user_pdb_ID):
    diagram3d.logger = logging.getLogger("SNVdistro.mod.diagram3d")
    from scipy import stats
    from pymol import cmd
    ###############  getuniprot  ###############
    list_uniprot = getuniprot(uniprot_id)
    for seq_loc in range(len(list_uniprot)):
        if list_uniprot[seq_loc].startswith("SQ"):
            break
    refseq =  re.findall("[A-Z]+", "".join(list_uniprot[seq_loc+1:]).replace(" ", ""))[0]
    diagram3d.logger.info("\nreference sequence:\n" +refseq)
    ###############  get pdb sequence ID  ###############
    if user_pdb_ID=="":
        ################  pdbblast  ###############
        pdbblast_result = pdbblast(refseq)
        best_match_pdb = pdbblast_result[0][0]
        pdb_ID = best_match_pdb[:4]    ##to put it in pymol
        pdb_chain = best_match_pdb[5:]  ##to use it in getpdbsw -s ID
    else:
        pdb_ID = user_pdb_ID[:4] #requires user to type in the chain
        pdb_chain = user_pdb_ID[4:]
    #################　fetch mmcif file from PDB　###############
    from Bio.PDB import PDBList
    pdbl=PDBList()
    cif_file_path = pdbl.retrieve_pdb_file(pdb_ID, file_format="mmCif", pdir = res_loc) #get pdb
    ################  get pdb sequence  ###############
    pdbseq = atom2seq(cif_file_path)[pdb_chain]
    #################  get alignment to reference sequence  ###############
    from Bio import Align
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(refseq, pdbseq)
    alignment = alignments[0]
    aliloc = alignment.aligned
    ################# alignment log #################
    diagram3d.logger.info("Alignment Score = " + str(alignment.score))
    ali_for = alignment.format().split('\n')
    split_num = 60
    ali_ref = [ali_for[0][i:i+split_num] for i in range(0, len(ali_for[0]), split_num)]
    ali_ali = [ali_for[1][i:i+split_num] for i in range(0, len(ali_for[1]), split_num)]
    ali_tar = [ali_for[2][i:i+split_num] for i in range(0, len(ali_for[2]), split_num)]
    log_ali = ""
    for i in range(len(ali_ref)):
        log_ali = log_ali + str(ali_ref[i])+"\n"+str(ali_ali[i])+"\n"+str(ali_tar[i])+"\n\n"
    diagram3d.logger.info("Alignment\n"+log_ali)
    #################　make dict to convert　###############
    aliloc_dict = {}
    for i, tup in enumerate(aliloc[0]):
        for j in range(tup[0], tup[1]):
            aliloc_dict[j+1] = aliloc[1][i][0] + j - tup[0] + 1
    aliloc_dict = {v: k for k, v in aliloc_dict.items()}
    #################　reindex cif file　###############
    cif_file = pd.read_table(cif_file_path, header=None)
    modi_cif_file = cif_file[0].str.split(expand=True)
    #################　first, reset residue position so that it starts from one　###############
    def rank_list(nums):
        ranks = {}
        rank = 1
        for num in set(nums):
            ranks[num] = rank
            rank += 1
        return [ranks[num] for num in nums]
    cha = pdb_chain.upper()
    mask = ((modi_cif_file.iloc[:, 0] == "ATOM") | (modi_cif_file.iloc[:, 0] == "HETATM")) & (modi_cif_file.iloc[:, 18] == cha)
    modi_cif_file.loc[mask, 16] = rank_list(modi_cif_file[((modi_cif_file[0]==("ATOM"))|(modi_cif_file[0]==("HETATM"))) & (modi_cif_file[18] == cha)][16].astype(int))
    #################　second, reindex residue position so that the amino acid position is aligned to the reference sequence　###############
    modi_cif_file.loc[mask] = modi_cif_file.loc[mask].replace({16:aliloc_dict})
    modi_cif_file.to_csv(res_loc+"/"+pdb_ID+'_reindexed.cif',sep="\t", header=False, index=False) #save reindexed_cif file
    ################# workout variations within the range of pdb file  ###############
    exsnv["pdb_POS"] = ""
    for i in range(len(exsnv)):
        for j in range(len(aliloc[0])):
            if aliloc[0][j][0] <= exsnv["p_POS"][i] < aliloc[0][j][1]:
                exsnv["pdb_POS"][i] = exsnv["p_POS"][i]
    ############### match frequency to heatmap ###############
    freq = exsnv["pdb_POS"].mask(exsnv["pdb_POS"]=="").value_counts().to_dict()
    for i in range(1,aliloc[1][-1][1]-1):
        if i not in freq:
            freq[i] = 0
    npfreq = np.array(list(freq.items()), dtype=object)
    z = stats.zscore(npfreq[:, 1].astype(int))
    freq_z = np.append(npfreq, z.reshape(len(z),1),axis=1)
    freq_z = pd.DataFrame(freq_z)
    freq_z["rgb"] = ""
    for j in range(len(freq_z.loc[freq_z[1] != 0])):
        freq_z["rgb"][j] = rgb(-abs(freq_z[2]).max(), abs(freq_z[2]).max(), freq_z[2][j])
    ############### pymol ###############
    cmd.load(res_loc+"/"+pdb_ID+'_reindexed.cif')
    cmd.color("green",  pdb_ID)
    f = freq_z["rgb"].drop_duplicates().reset_index(drop="True")
    f = f[f.astype(bool)]
    for colour in range(len(f)):
        cmd.set_color(str(f[colour]), f[colour])
    for k in range(len(freq_z.loc[freq_z[1] != 0])):
        cmd.color(str(freq_z["rgb"][k]),"chain "+ cha +" & resi "+str(freq_z[0][k]))
    cmd.color("blue",  "color green & chain "+ cha)
    cmd.color("grey80",  "color green")
    cmd.set("cartoon_transparency", 0.5, "color grey80")
    cmd.set("stick_transparency", 0.5, "color grey80")
    cmd.show("wire", "all")
    cmd.save(filename=res_loc+"/"+pdb_ID+"_mutfreq.pse", selection =pdb_ID , state="-1")
    return print("pymol file saved!")
