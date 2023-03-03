'''
please edit the parameters below according to  your environment
'''

#location of database files
database_loc = {"CCDS": "/db/NCBI/CCDS/CCDS.20180614.txt", \
            "ClinVar": "/db/ClinVar/clinvar.vcf", \
            "dbSNP": "/db/dbSNP/GCF_000001405.39", \
            "gnomAD": "/db/gnomAD3.1.2/", \
            "swissprot": "/db/swissprot/uniprot_sprot.dat", \
            "trembl": "/db/swissprot/uniprot_trembl.dat", \
            "GRCh38":  {"1": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa", \
                        "2": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.2.fa", \
                        "3": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa", \
                        "4": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.4.fa", \
                        "5": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa", \
                        "6": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa", \
                        "7": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa", \
                        "8": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.8.fa", \
                        "9": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.9.fa", \
                        "10": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.10.fa", \
                        "11": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa", \
                        "12": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.12.fa", \
                        "13": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa", \
                        "14": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa", \
                        "15": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.15.fa", \
                        "16": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.16.fa", \
                        "17": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa", \
                        "18": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.18.fa", \
                        "19": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa", \
                        "20": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa", \
                        "21": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa", \
                        "22": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa", \
                        "X": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa", \
                        "Y": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa", \
                        "MT": "/db/ensembl/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa", }}

#select the database
database_to_use = { "ClinVar": False, \
                    "dbSNP": True, \
                    "gnomAD": True}

#select what clinical significnace to be included
clnsig_to_include = {"Benign": True, \
                    "Likely_benign": True, \
                    "Uncertain_significance": True, \
                    "Likely_pathogenic": False, \
                    "Pathogenic": False, \
                    "NaN" : True, \
                    "Conflicting_interpretations_of_pathogenicity": False,}

#maybe a way of putting a weight of clinical significance in ClinVar
#perhaps it might be better to add output file location on here?


#location of hash table to be stored  <= is this necessary? I could preset the program
#to store it under where hashtb.py file is...
hashtb_loc = {"ClinVar": "/db/ClinVar/ClinVar.hashtb", \
            "dbSNP": "/db/dbSNP/dbSNP.hashtb", \
            "gnomAD": "/db/gnomAD3.1.2/TABLE/", \
            "swissprot": "/db/swissprot/key/swissprot/", \
            "trembl": "/db/swissprot/key/trembl/"}

#location of separate functions



getuniprot_loc = ""
pdbblast_loc = ""
getpdbsw = ""

#This parameter changes the amount of incliments when forming hash table of the database (bin/ jump window).
#e.g. if set to 100,000, it will store the bite location of mutations at every 100,000 residues
hashtb_DIV = {"ClinVar": 100000, \
            "dbSNP": 100000, \
            "gnomAD": 50000}
