from sys import *
import os
import math
sys.path.append(os.path.dirname(os.getcwd()))
from control_file import *

#chr = [ '1', '2', '3', '4', '5', '6', '7', '8', '9','10','11','12',\
       #'13','14','15','16','17','18','19','20','21','22','23', 'X',\
        #'Y','MT']
DIV = hashtb_DIV["gnomAD"]
nt  = 0
pt = 0
l = []
c_num  = str(argv[1])

def myround(x, base=DIV):
    return base * math.floor(x/base)

with open(database_loc["gnomAD"]+"gnomad.genomes.v3.1.2.sites.chr"+c_num+".vcf","r", encoding="latin") as f:
    for line in f:
        if line[0] == '#':
            pt += len(line)
        else:
            tmp = line.split()
            if nt == 0:
                nt = myround(int(tmp[1]))
            if myround(int(tmp[1])) == nt:
                l.append([str(nt),str(pt)])
                nt += DIV
            elif myround(int(tmp[1])) > nt:
                while myround(int(tmp[1])) > nt:
                    nt += DIV
                l.append([str(nt),str(pt)])
                nt += DIV
            pt += len(line)

with open(hashtb_loc["gnomAD"]+'gnomAD.chr'+c_num+'.hashtb', 'w') as file:
    for row in l:
        file.write(' '.join([str(item) for item in row]))
        file.write('\n')
