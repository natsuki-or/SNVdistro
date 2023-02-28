#!/usr/bin/env python

from sys import *
from math import *

DIRECTORY = './key/trembl/'
#
# function: selectfile()
def selectfile(id,length,pt,header,idfile,number,maximum):
    """store id/ac and pointer to idfile list """
    ok = 0
    idlen = len(id)
    idx = 0
    for i in range(0,idlen):
        if header.find(id[i]) >= 0:
            idx += number[header.index(id[i])]
    idx = idx % maximum
    idfile[idx].extend([[id,pt]])
    return idfile

#
# function: init()
def init():
    header = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
    number = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157]
    length = len(header)
    maximum    = 2047
    idfile = [[] for i in range(0,maximum+1)]
    for p in range(0,maximum+1):
        filen = DIRECTORY
        filen += str(p)
        f = open(filen,"w")
        f.close()

    return header,length,idfile,number,maximum

#
# function: printall()
def printall(idlen,maximum):
    for p in range(0,maximum+1):
        filen = DIRECTORY
        filen += str(p)
        printout(filen,idfile[p])

    return
#
# function: printout()
def printout(filen,ids):
    f = open(filen,'a')
    for e in ids:
        f.write(e[0]+'\t'+str(e[1])+'\n')
    f.close()

    return

# main routine
# argv[1]    uniprot_trembl.dat
#
if len(argv) != 2:
   exit()

print ('Building hash table.')
header,length,idfile,number,maximum = init()

print ('Reading %s' % (argv[1]))
f = open(argv[1])
pt = 0
pre = 0
while True:
    l = f.readline()
    if l == '':
        break
    l = l.strip()
    if l[0:3] == 'ID ':
        pt = f.tell()-(len(l)+1)
        id = l.split()[1].upper()
        idfile = selectfile(id,length,pt,header,idfile,number,maximum)
    if l[0:3] == 'AC ':
        ac = l[2:].strip(' ').upper().split(';')
        for x in ac:
            x = x.strip()
            if len(x) != 0:
                idfile = selectfile(x,length,pt,header,idfile,number,maximum)
    if pt - pre >  100000000:
        printall(idfile,maximum)
        idfile = [[] for i in range(0,maximum+1)]
        print ('\b%c' % '.'),
        stdout.flush()
        pre = pt

f.close()
print
print ('Closing %s' % (argv[1]))
printall(idfile,maximum)
print ('Hash table completed.')
#EOF
