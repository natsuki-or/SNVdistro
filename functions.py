import os
import sys
sys.path.append(os.getcwd())
from control_file import *

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import re
from Bio import SeqIO
from Bio.Seq import MutableSeq
import subprocess as sp
import logging
module_logger = logging.getLogger("SNVdistro.mod")


# function: selectfile()
def selectfile(id,header,idfile,idlist,number,maximum):
    idlen = len(id)
    idx = 0
    for i in range(0,idlen):
        if header.find(id[i]) >= 0:
            idx += number[header.index(id[i])]
    idx = idx % maximum
    if idx not in idfile:
        idfile.append(idx)
        idlist.append([id])
    else:
        idlist[idfile.index(idx)].append(id)
    return idfile,idlist

#
# function: retrieve()
def retrieve(db,idfile,idlist,path,header,output):
    """ Retrieve Uniprot Data """
    g = open(db,'r')
    length = len(idfile)
    for num in range(0,length):
        if len(idlist[num]) == 0:
            continue

        f = open(path+str(idfile[num]),'r')
        dic = {}
        for l in f:
            tmp = l.strip().split()
            dic[tmp[0]] = int(tmp[1])
        f.close()

        tmp = idlist[num][:]
        for id in idlist[num]:
            if id in dic:
                tmp.remove(id)
                pt = dic[id]
                g.seek(pt)
                for l in g:
                    output.append(l)
                    if l[0:2] == '//':
                        break
        idlist[num] = tmp
    g.close()

    return idlist, output

#
# function: init()
def init():
    header = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
    number = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157]
    maximum    = 2047
    return header,number,maximum

#
# main function
def getuniprot(ids):
    header, number, maximum = init()
    idfile = []
    idlist = []
    output = []
    for id in ids:
        idfile,idlist = selectfile(id,header,idfile,idlist,number,maximum)

    db   = database_loc["swissprot"]
    path = '/db/swissprot/key/swissprot/'
    idlist,output = retrieve(db,idfile,idlist,path,header,output)

    empty = 0
    for ids in idlist:
        if len(ids) != 0:
            empty += len(ids)

    if empty == 0:
        return output

    db   = database_loc["trembl"]
    path = '/db/swissprot/key/trembl/'
    idlist,output = retrieve(db,idfile,idlist,path,header,output)

    for ids in idlist:
        if len(ids) != 0:
            for e in ids:
                print ('%s does not exist' % (e), file=stderr)

    return output
#EOF



def Up_CCDS_ID(uniprot_name):
    Up_CCDS_ID.logger = logging.getLogger("SNVdistro.mod.Up_CCDS_ID")
    list_uniprot = getuniprot([uniprot_name])
    ccds_id_loc = [n for n, l in enumerate(list_uniprot) if l.startswith('DR   CCDS;')]
    ccds_id = re.findall('; (.*);', list_uniprot[ccds_id_loc[0]])
    return ccds_id


def Up_PDB_ID(uniprot_name):
    list_uniprot = getuniprot([uniprot_name])
    PDB_id_loc = [n for n, l in enumerate(list_uniprot) if l.startswith('DR   PDB;')]
    PDB_id = re.findall('; (.*);', list_uniprot[PDB_id_loc[0]])
    return PDB_id


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
    return g_snv


def dbSNP(C_num, g_start, g_stop):
    # import pandas as pd
    # pd.options.mode.chained_assignment = None
    # import re
    dbSNP.logger = logging.getLogger("SNVdistro.mod.dbSNP")
    chr_an = {"NC_000001.11":1, "NC_000002.12":2, "NC_000003.12":3, "NC_000004.12":4, "NC_000005.10":5, \
            "NC_000006.12":6, "NC_000007.14":7, "NC_000008.11":8, "NC_000009.12":9, "NC_000010.11":10, \
            "NC_000011.10":11, "NC_000012.12":12, "NC_000013.11":13, "NC_000014.9":14, "NC_000015.10":15, \
            "NC_000016.10":16, "NC_000017.11":17, "NC_000018.10":18, "NC_000019.10":19, "NC_000020.11":20, \
            "NC_000021.9":21, "NC_000022.11":22, "NC_000023.11":"X", "NC_000024.10":"Y", "NC_012920.1":"MT"}
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
    v.drop('INFO', inplace=True, axis=1)
    v["CHROM"] = chr_an[v["CHROM"][0]]
    g_snv = v[v["CLNVC"]=="SNV"]
    g_snv.drop('CLNVC', inplace=True, axis=1)
    g_snv = g_snv.reset_index(drop=True)
    for ls in range(len(g_snv)):
        g_snv["ALT"][ls]= g_snv["ALT"][ls].split(",")
    g_snv = g_snv.explode("ALT")
    return g_snv


def gnomAD_snv(C_num, g_start, g_stop):
    #from sys import *
    # import pandas as pd
    # pd.options.mode.chained_assignment = None
    # import re
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


def mut_amino(df,g_snv,g_start, g_stop,g):
    # import pandas as pd
    # pd.options.mode.chained_assignment = None
    # from Bio import SeqIO
    # from Bio.Seq import MutableSeq
    # import numpy as np
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


def diagram2d(exsnv,b):
    # import pandas as pd
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


def diagram3d(upname,res_loc,exsnv,pdb_ID_chain):
    from scipy import stats
    from pymol import cmd
    ###############  getuniprot  ###############
    list_uniprot = getuniprot([upname])
    for s in range(len(list_uniprot)):
        if list_uniprot[s].startswith("SQ"):
            break
    refseq = list_uniprot[s:-1]
    with open(res_loc+"/"+upname+'.faa', 'w') as file:
        for row in refseq:
            file.write(''.join([str(item) for item in row]))
            file.write('\n')
    if pdb_ID_chain=="":
        ################  pdbblast  ###############
        sp.run(['pdbblast', res_loc+"/"+upname+'.faa', res_loc+"/"+upname+'.blast'])
        ################  open blast result and get pdb code ###############
        with open(res_loc+"/"+upname+'.blast', 'r') as f:
            lines = f.readlines()
            #print(lines)
        for i in range(len(lines)):
            if lines[i].startswith(">"):
                break
        for j in range(i+1, len(lines)):
            if lines[j].startswith(">"):
                break
        first_pdb = lines[i:j-1]
        pdb_ID_chain = first_pdb[0][5:10]  ##to use it in getpdbsw -s ID
        pdb_ID = first_pdb[0][5:9]    ##to put it in pymol
    else:
        pdb_ID = pdb_ID_chain[0:4]    #requires user to type in the chain
    ################  get pdb sequence  ###############
    cmd2 = ['getpdbsw', '-s', pdb_ID_chain]
    p2 = sp.Popen(cmd2, stdout=sp.PIPE)
    output2 = p2.stdout.read().decode('utf-8')
    list2 = output2.split("\n")
    for s in range(len(list2)):
        if list2[s].startswith("SQ"):
            break
    pdbseq = "".join(list2[s+1:-2]).replace(" ", "")
    #################  get alignment to reference sequence  ###############
    with open(res_loc+"/"+upname+".faa") as f:
        lines = f.read().splitlines()
    refseq =  "".join(lines[1:-1]).replace(" ", "")
    from Bio import Align
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(refseq, pdbseq)
    alignment = alignments[0]
    aliloc = alignment.aligned
    #################　make dict to convert　###############
    aliloc_dict = {}
    for i, tup in enumerate(aliloc[0]):
        for j in range(tup[0], tup[1]):
            aliloc_dict[j+1] = aliloc[1][i][0] + j - tup[0] + 1
    aliloc_dict = {v: k for k, v in aliloc_dict.items()}
    #################　fetch mmcif file from PDB　###############
    from Bio.PDB import PDBList
    pdbl=PDBList()
    cif_file_path = pdbl.retrieve_pdb_file(pdb_ID, file_format="mmCif", pdir = res_loc) #get pdb
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
    cha = pdb_ID_chain[-1].upper()
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
    cmd.show("wire", "all")
    cmd.save(filename=res_loc+"/"+pdb_ID+"_mutfreq.pse", selection =pdb_ID , state="-1")
    return print("pymol file saved!")
