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

def Up_CCDS_ID(uniprot_name):
    Up_CCDS_ID.logger = logging.getLogger("SNVdistro.mod.Up_CCDS_ID")
    cmd1 = ['getuniprot', '-s', uniprot_name]
    p = sp.Popen(cmd1, stdout=sp.PIPE)
    output = p.stdout.read().decode('utf-8')
    list_uniprot = output.split("\n")
    ccds_id_loc = [n for n, l in enumerate(list_uniprot) if l.startswith('DR   CCDS;')]
    ccds_id = re.findall('; (.*);', list_uniprot[ccds_id_loc[0]])
    return ccds_id


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


def COSMIC_snv(gene,C_num, g_start, g_stop):
    """
    This is incomplete.
    The original COSMIC file is not sorted in anyway. (or maybe I just couldn't figure out how it was sorted)
    There is an alphabetically sorted COMIC file available at /work13/natsuki/Cosmic_sort.tsv
    However, it lacked consistency over how mutation location were stored in the database
    and it was quite challenging to use. As I felt the other 3 was sufficient to capture most
    SNVs, COSMIC was not included in the search in the end.
    """
    # import pandas as pd
    # pd.options.mode.chained_assignment = None
    # import numpy as np
    ftls = []
    pt_l = []
    pt_b = []
    with open ("/work13/natsuki/hashtb/COSMIC.hashtb","r") as f:
        for line in f:
            tmp = line.strip().split(",")
            ftls.append(tmp[1])
            pt_l.append(int(tmp[2]))
            pt_b.append(int(tmp[3]))
    total = len(ftls)
    t_ftls = gene[0:2]
    for i in range(total):
        if t_ftls.casefold() == "A1".casefold():
            start_l = 0
            stop_l = pt_l[1]
            start_b = 0
            stop_b = pt_b[2]
            break
        if t_ftls.casefold() == ftls[i].casefold():
            start_l = pt_l[i-1]
            stop_l = pt_l[i]
            start_b = pt_b[i-1]
            stop_b = pt_b[i]
            break
    with open("/work13/natsuki/Cosmic_sort.tsv","r", encoding="latin") as f:
        f.seek(start_b)
        doc = f.read(stop_b-start_b)
    doc = doc.split("\n")
    results = []
    for j in range(len(doc)):
        if doc[j].casefold().startswith(gene.casefold()):
            results.append(doc[j])
    #################
    df = pd.DataFrame(results)
    df = df[0].str.split("\t", expand=True)
    df.columns = ["Gene name", "Accession Number", "Gene CDS length", "HGNC ID", \
                "Sample name", "Sample id", "ID_tumour", "Primary site", "Site subtype 1", \
                "Site subtype 2", "Site subtype 3", "Primary histology", "Histology subtype 1", \
                "Histology subtype 2", "Histology subtype 3", "Genome-wide screen", \
                "GENOMIC_MUTATION_ID", "LEGACY_MUTATION_ID", "MUTATION_ID", "Mutation CDS", \
                "Mutation AA", "Mutation Description", "Mutation zygosity", "LOH", "GRCh", \
                "Mutation genome position", "Mutation strand", "Resistance Mutation", \
                "Mutation somatic status", "Pubmed_PMID", "ID_STUDY", "Sample Type", \
                "Tumour origin", "Age", "HGVSP", "HGVSC", "HGVSG"]
    #print(df)
    df_loc_data = df.loc[:,["Gene name", "Accession Number", "Mutation CDS", "Mutation AA", \
                            "Mutation Description", "GRCh", "Mutation genome position", \
                            "Mutation strand","HGVSG"]]
    df_loc_data.replace('', np.nan, inplace=True)
    v = df_loc_data.dropna(subset=["Mutation CDS", "Mutation AA", "Mutation Description", \
                        "GRCh", "Mutation genome position", "Mutation strand","HGVSG"], how="all")
    #######################
    ################
    ###########
    #######
    ####
    ##
    #
    g_snv = v[v["Mutation Description"]=="single_nucleotide_variant"]
    g_snv.drop('CLNVC', inplace=True, axis=1)
    g_snv = g_snv.reset_index(drop=True)
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


def diagram3d(upname,res_loc,exsnv,pdb_code):
    from scipy import stats
    #import subprocess as sp
    # import pandas as pd
    # import numpy as np
    from pymol import cmd
    ###############  getuniprot  ###############
    #upname = args.uniprot_name
    cmd1 = ['getuniprot', '-s', upname]
    p = sp.Popen(cmd1, stdout=sp.PIPE)
    output = p.stdout.read().decode('utf-8')
    list_uniprot = output.split("\n")
    for s in range(len(list_uniprot)):
        if list_uniprot[s].startswith("SQ"):
            break
    refseq = list_uniprot[s:-1]
    with open(res_loc+"/"+upname+'.faa', 'w') as file:
        for row in refseq:
            file.write(''.join([str(item) for item in row]))
            file.write('\n')
    if pdb_code=="":
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
        pdb_codeA = first_pdb[0][5:10]  ##to use it in getpdbsw -s ID
        pdb_code = first_pdb[0][5:9]    ##to put it in pymol
    else:
        pdb_codeA = pdb_code + "A"
    ################  get pdb sequence  ###############
    cmd2 = ['getpdbsw', '-s', pdb_codeA] #2L7BA
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
    # from Bio import SeqIO
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(refseq, pdbseq)
    alignment = alignments[0]
    aliloc = alignment.aligned
    ################# convert residue position to that of pdb  ###############
    exsnv["pdb_POS"] = ""
    for i in range(len(exsnv)):
        for j in range(len(aliloc[0])):
            if aliloc[0][j][0] <= exsnv["p_POS"][i] < aliloc[0][j][1]:
                dif = aliloc[0][0][0]
                while j>0:
                    dif += aliloc[0][j][0]-aliloc[0][j-1][1]
                    j -= 1
                exsnv["pdb_POS"][i] = exsnv["p_POS"][i] - dif
    ###############  get pdb structure file  ###############
    from Bio.PDB import PDBList
    pdbl=PDBList()
    cif_file_path = pdbl.retrieve_pdb_file(pdb_code, file_format="mmCif", pdir = res_loc) #get pdb
    cif_file = pd.read_table(cif_file_path, header=None)
    subtitle = ["data", "_entry.id", "_cell.", "_symmetry.", "loop_", "_atom_site.", "ATOM"] #extracts nessesary info to feed it to pymol
    str_file = cif_file[cif_file[0].str.startswith(tuple(subtitle))]
    str_file = str_file[0].str.split(expand=True).reset_index(drop=True)
    for n in range(len(str_file)):
        if str_file[0][n].startswith("ATOM"):
            break
    for m in range(n,len(str_file)):
        if str_file[0][m].startswith("loop_"):
            break
    str_file.drop(str_file.tail(len(str_file)-m).index, inplace=True) #gets rid of extra "loop_" at the end
    for t in range(n,len(str_file)):
        if str_file[20][t]>"1":
            break
    str_file.drop(str_file.tail(len(str_file)-t).index, inplace=True) #gets rid of other model numbers
    ###############  re-index amino acid number  ###############
    str_file = str_file[(str_file[7].isnull()) | (str_file[7]=="1")]
    idx = 0
    for i in range(n-1,len(str_file)):
        try:
            if str_file[16][i] == str_file[16][i+1]:
                str_file[16][i]= idx
            else:
                idx += 1
                str_file[16][i]= idx
        except Exception:
            pass
    str_file[16] = str_file[16].shift(1, fill_value=None)
    str_file.fillna("",inplace=True)
    s = str_file.values.tolist()
    with open(res_loc+"/"+pdb_code+'_reindexed.cif', 'w') as file:
        for row in s:
            file.write(' '.join([str(item) for item in row]))
            file.write('\n')
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
    cmd.load(res_loc+"/"+pdb_code+'_reindexed.cif')
    cmd.color("green",  pdb_code)
    f = freq_z["rgb"].drop_duplicates().reset_index(drop="True")
    f = f[f.astype(bool)]
    for colour in range(len(f)):
        cmd.set_color(str(f[colour]), f[colour])
    for k in range(len(freq_z.loc[freq_z[1] != 0])):
        cmd.color(str(freq_z["rgb"][k]),"resi "+str(freq_z[0][k]))
    cmd.color("blue",  "color green")
    cmd.show("wire", "all")
    cmd.save(filename=res_loc+"/"+pdb_code+"_mutfreq.pse", selection =pdb_code , state="-1")
    return print("pymol file saved!")
