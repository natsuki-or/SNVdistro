from sys import *
import os
import math
sys.path.append(os.path.dirname(os.getcwd()))
from control_file import *

chr = [ "NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12", "NC_000005.10", \
        "NC_000006.12", "NC_000007.14", "NC_000008.11", "NC_000009.12", "NC_000010.11", \
        "NC_000011.10", "NC_000012.12", "NC_000013.11", "NC_000014.9", "NC_000015.10", \
        "NC_000016.10", "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11", \
        "NC_000021.9", "NC_000022.11", "NC_000023.11", "NC_000024.10", "NC_012920.1"]


DIV = hashtb_DIV["dbSNP"]
nt  = [0 for i in range(26)]

pt = 0
l = []

def myround(x, base=DIV):
    return base * math.floor(x/base)


with open(database_loc["dbSNP"],"r",) as f:
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
                    #print (tmp[0],nt[loc],pt)
                    l.append([tmp[0],str(nt[loc]),str(pt)])
                    nt[loc] += DIV
                elif myround(int(tmp[1])) > nt[loc]:
                    while myround(int(tmp[1])) > nt[loc]:
                        nt[loc] += DIV
                    #print (tmp[0],nt[loc],pt)
                    l.append([tmp[0],str(nt[loc]),str(pt)])
                    nt[loc] += DIV
                pt += len(line)
        except Exception:
            pass


with open(hashtb_loc["dbSNP"], 'w') as file:
    for row in l:
        file.write(' '.join([str(item) for item in row]))
        file.write('\n')
