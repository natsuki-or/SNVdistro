'''
please edit the parameters below according to  your environment
'''
 
#location of database files
database_loc = {"CCDS": "/db/NCBI/CCDS/CCDS.20180614.txt", \
            "ClinVar": "/db/ClinVar/clinvar.vcf", \
            "dbSNP": "/db/dbSNP/GCF_000001405.39", \
            "gnomAD": "/db/gnomAD3.1.2/", \
            "usersdata": "/db/usersdata.vcf", \
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
database_to_use = { "ClinVar": True, \
                    "dbSNP": True, \
                    "gnomAD": True,
                    "usersdata": False}

#select what clinical significnace to be included from ClinVar
clnsig_to_include = {"Benign": False, \
                    "Likely_benign": False, \
                    "Uncertain_significance": False, \
                    "Likely_pathogenic": True, \
                    "Pathogenic": True, \
                    "NaN" : False, \
                    "Conflicting_interpretations_of_pathogenicity": False,}


#set an arbitary unit form -1 to 1 for each clinical significnace labelling on dbSNP
clnsig_scale = {"Uncertain significance": 0,\
                "not provided": 0,\
                "Benign": -1,\
                "Likely benign": -0.9,\
                "Likely pathogenic": 0.9,\
                "Pathogenic": 1,\
                "drug response": 0,\
                "confers sensitivity": 0,\
                "risk-factor": 0.7,\
                "association": 0.5,\
                "protective": -1,\
                "conflict": 0,\
                "affects": 0.3,\
                "other": 0}

#select thresholds for the arbitary unit for dbSNP(inclusive)
clnsig_threshold = {"upper_limit": 1,\
                    "lower_limit": 0}

clnsig_func = {"mean": True,\
                "max": False,\
                "min": False}


#select thresholds for allele frequency for gnomAD(inclusive)
allele_freq_threshold = {"upper_limit": 0.01,\
                        "lower_limit": 0.00}


#location of hash table to be stored  <= is this necessary? I could preset the program
#to store it under where hashtb.py file is...
hashtb_loc = {"ClinVar": "/db/ClinVar/ClinVar.hashtb", \
            "dbSNP": "/db/dbSNP/dbSNP.hashtb", \
            "gnomAD": "/db/gnomAD3.1.2/TABLE/", \
            "swissprot": "/db/swissprot/key/swissprot/", \
            "trembl": "/db/swissprot/key/trembl/"}




#This parameter changes the amount of incliments when forming hash table of the database (bin/ jump window).
#e.g. if set to 100,000, it will store the bite location of mutations at every 100,000 residues
hashtb_DIV = {"ClinVar": 100000, \
            "dbSNP": 100000, \
            "gnomAD": 50000}
