'''
please edit the parameters below according to  your environment
'''

#location of files

intern_func_loc = "/work13/natsuki/python"

database_loc = {"CCDS": "/db/NCBI/CCDS/CCDS.20180614.txt", \
            "GRCh38": "/db/ensembl/Homo_sapiens/dna/", \
            "ClinVar": "/db/ClinVar/clinvar.vcf", \
            "dbSNP": "/db/dbSNP/GCF_000001405.39", \
            "gnomAD": "/db/gnomAD3.1.2/", \
            "COSMIC": "/work13/natsuki/Cosmic_sort.tsv"}

hashtb_loc = {"ClinVar": "/db/ClinVar/ClinVar.hashtb", \
            "dbSNP": "/db/dbSNP/dbSNP.hashtb", \
            "gnomAD": "/db/gnomAD3.1.2/TABLE/", \
            "COSMIC": "/work13/natsuki/hashtb/COSMIC.hashtb"}

"""
for GRCh38 and gnomAD, which the user may store them differently, how should I locate them?
e.g. in our cluster, GRCh38 fasta file are stored as:
        Homo_sapiens.GRCh38.dna.chromosome.1.fa
        Homo_sapiens.GRCh38.dna.chromosome.2.fa
        Homo_sapiens.GRCh38.dna.chromosome.3.fa
    etc...
        Homo_sapiens.GRCh38.dna.chromosome.X.fa
        Homo_sapiens.GRCh38.dna.chromosome.Y.fa
        Homo_sapiens.GRCh38.dna.chromosome.MT.fa
    so I can access them by chromosome by
        "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."+str(C_num)+".fa"

    but if in another place it was stored like
        Homo_sapiens.GRCh38.dna.chromosome.01.fa
        Homo_sapiens.GRCh38.dna.chromosome.02.fa
        Homo_sapiens.GRCh38.dna.chromosome.03.fa (to stop chr10 coming before chr2)
    or if it was not separated by chromosome
    or if it was separated by chromosome but in different folder
    etc...
    I wouldn't be able to access it and not sure how I could accommodate for that.
    maybe make an alias folders?

"""


getuniprot_loc = ""
pdbblast_loc = ""
getpdbsw = ""

#This parameter changes the amount of incliments when forming hash table of the database.
#e.g. if set to 10000, it will store the bite location
ClinVar_hashtb_DIV = 100000
dbSNP_hashtb_DIV = 100000
gnomAD_hashtb_DIV = 50000
