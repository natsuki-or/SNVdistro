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

|Name|URL|Reference|
|----|---|-----------------|
|The Consensus Coding Sequence (CCDS)|https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi|(1)|
|Genome Reference Consortium Human Build 38 (GRCh38)|https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40|(2)|
|ClinVar|https://www.ncbi.nlm.nih.gov/clinvar/|(3)|
|The Single Nucleotide Polymorphism Database (dbSNP)|https://www.ncbi.nlm.nih.gov/snp/|(4)|
|The Genome Aggregation Database (gnomAD)|https://gnomad.broadinstitute.org/|(5)|
|Biopython|https://biopython.org/|(6)|
|Matplotlib|https://matplotlib.org/stable/index.html|(7)|
|NumPy|https://numpy.org/|(8)|
|Pandas|https://pandas.pydata.org/docs/#|(9)|
|SciPy|https://scipy.org/|(10)|
|PyMOL|https://github.com/schrodinger/pymol-open-source|(11)|


To run the program use the following syntax
   $ python3 main.py {Uniprot ID} {location of output folder} {type of diagram (2d or 3d)} {pdb_ID (optional)}

- The first argument takes the UniProt ID(e.g. BRAF_HUMAN) or accession ID(e.g.P15056). 
- The second argument takes the location user wishes to store the outputs. 
- The third argument takes the type of diagram the user wishes to create. 
   - Whilst we believe the main strength of the program lies in 3D visualisation, the program can also be used to visualise the frequency and distribution of SNVs as a histogram. Selecting “2d” outputs a histogram and “3d” outputs a pymol heatmap. 
- The fourth argument is optional. If user select 3d, they can choose a structural file to map the variations onto by entering the PDB ID of choice.


Reference
1. Pujar S, O’Leary NA, Farrell CM, Loveland JE, Mudge JM, Wallin C, et al. Consensus coding sequence (CCDS) database: A standardized set of human and mouse protein-coding regions supported by expert curation. Nucleic Acids Res. 2018;46(D1):D221–8.
2. International Human Genome Sequencing Consortium. Finishing the euchromatic sequence of the human genome. Nature [Internet]. 2004 Oct;431(7011):931–45. Available from: http://www.nature.com/articles/nature03001
3. Landrum MJ, Chitipiralla S, Brown GR, Chen C, Gu B, Hart J, et al. ClinVar: Improvements to accessing data. Nucleic Acids Res. 2020;48(D1):D835–44.
4. Sherry ST, Ward M, Sirotkin K. dbSNP - database for single nucleotide polymorphisms and other classes of minor genetic variation. Genome Res. 1999;9(8):677–9.
5. Chen S, Francioli LC, Goodrich JK, Collins RL, Kanai M, Wang Q, et al. A genome-wide mutational constraint map quantified from variation in 76,156 human genomes. bioRxiv [Internet]. 2022 Jan 1;2022.03.20.485034. Available from: http://biorxiv.org/content/early/2022/10/10/2022.03.20.485034.abstract
6. Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, et al. Biopython: Freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422–3.
7. Hunter JD. Matplotlib: A 2D Graphics Environment. Comput Sci Eng [Internet]. 2007;9(3):90–5. Available from: http://ieeexplore.ieee.org/document/4160265/
8. Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, et al. Array programming with NumPy. Nature [Internet]. 2020;585(7825):357–62. Available from: http://dx.doi.org/10.1038/s41586-020-2649-2
9. McKinney W. Data Structures for Statistical Computing in Python. Proc 9th Python Sci Conf. 2010;1(Scipy):56–61.
10. Virtanen P, Gommers R, Oliphant TE, Haberland M, Reddy T, Cournapeau D, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods. 2020;17(3):261–72.
11. Schrödinger L. The PyMOL Molecular Graphics System, Version 2.5.2 [Internet]. 2022. Available from: http://www.pymol.org/pymol

