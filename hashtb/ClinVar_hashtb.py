from sys import *
import os
import math
sys.path.append(os.path.dirname(os.getcwd()))
from control_file import *


chr = [ '1', '2', '3', '4', '5', '6', '7', '8', '9','10','11','12',\
       '13','14','15','16','17','18','19','20','21','22','23', 'X',\
        'Y','MT']
DIV = hashtb_DIV["ClinVar"]
nt  = [0 for i in range(26)]

pt = 0
l = []

def myround(x, base=DIV):
    return base * math.floor(x/base)

with open(database_loc["ClinVar"],"r", encoding="latin") as f:
    for line in f:
        try:
            if line[0] == '#':
                pt += len(line)
            else:
                tmp = line.split()
                loc = chr.index(tmp[0])
                if nt[loc] == 0:
                    nt[loc] = myround(int(tmp[1]))
                if tmp[0] in chr and myround(int(tmp[1])) == nt[loc]:
                    l.append([tmp[0],str(nt[loc]),str(pt)])
                    nt[loc] += DIV
                elif myround(int(tmp[1])) > nt[loc]:
                    while myround(int(tmp[1])) > nt[loc]:
                        nt[loc] += DIV
                    l.append([tmp[0],str(nt[loc]),str(pt)])
                    nt[loc] += DIV
                pt += len(line)
        except Exception:
            pass

with open(hashtb_loc["ClinVar"], 'w') as file:
    for row in l:
        file.write(' '.join([str(item) for item in row]))
        file.write('\n')
