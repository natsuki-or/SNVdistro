# SNVdistro
This program visualises the frequency, distribution and clinical significance of single nucleotide variations(SNV) of a gene of interest on 
a three-dimensional structure. It searches ClinVar, dbSNP and gnomAD for SNVs which are filtered according to users' interests and mapped
to a structural file retrieved from PDB. The resulting heatmap of SNV frequency, distribution and clinical significance can be viewed on
PyMOL. 

To set up
1. Install all files from this repository
2. Download prerequisites.
3. Open control_file.py and edit “database_loc” so that the location of databases in the user's environment is reflected.
4. Set “hashtb_loc” to wherever the user wishes to store the hash table for each database.
5. Either:
      - Move the hash table provided under hashtb to the location specified in the previous step.
      - Create hash tables for databases by running “database_hashtb.py” under hashtb. This may take a while.
6. Choose the database to search SNVs and filter to apply on control_file.py; the program is ready to run.


To run the program use the following syntax
   $ python3 main.py {Uniprot ID} {location of output folder} {type of diagram (2d or 3d)} {pdb_ID (optional)}

- The first argument takes the UniProt ID(e.g. BRAF_HUMAN) or accession ID(e.g.P15056). 
- The second argument takes the location user wishes to store the outputs. 
- The third argument takes the type of diagram the user wishes to create. 
   - Whilst we believe the main strength of the program lies in 3D visualisation, the program can also be used to visualise the frequency and distribution of SNVs as a histogram. Selecting “2d” outputs a histogram and “3d” outputs a pymol heatmap. 
- The fourth argument is optional. If user select 3d, they can choose a structural file to map the variations onto by entering the PDB ID of choice.
