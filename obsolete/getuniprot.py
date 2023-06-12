#!/usr/bin/env python3.6

from sys import *
from math import *
#from ast import literal_eval

# function: selectfile()
def selectfile(id,header,idfile,idlist,number,maximum):
    idlen = len(id)
    idx = 0
    for i in range(0,idlen):
        if header.find(id[i]) >= 0:
            idx += number[header.index(id[i])]
    idx = idx % maximum
    if idx not in idfile:
        idfile.append(idx)
        idlist.append([id])
    else:
        idlist[idfile.index(idx)].append(id)
    return idfile,idlist

#
# function: retrieve()
def retrieve(db,idfile,idlist,path,header,output):
    """ Retrieve Uniprot Data """
    g = open(db,'r')
    length = len(idfile)
    for num in range(0,length):
        if len(idlist[num]) == 0:
            continue

        f = open(path+str(idfile[num]),'r')
        dic = {}
        for l in f:
            tmp = l.strip().split()
            dic[tmp[0]] = int(tmp[1])
        f.close() 

        tmp = idlist[num][:]
        for id in idlist[num]:
            if id in dic:
                tmp.remove(id)                
                pt = dic[id]
                g.seek(pt)
                for l in g:
                    output.append(l)
                    if l[0:2] == '//':
                        break
        idlist[num] = tmp
    g.close()
    
    return idlist, output

#
# function: init()
def init():
    header = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
    number = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157]
    maximum    = 2047
    return header,number,maximum

#
# main function
def getuniprot(ids):
    header, number, maximum = init()
    idfile = []
    idlist = []
    output = []
    for id in ids:
        idfile,idlist = selectfile(id,header,idfile,idlist,number,maximum)
    
    db   = '/db/swissprot/uniprot_sprot.dat'
    path = '/db/swissprot/key/swissprot/'
    idlist,output = retrieve(db,idfile,idlist,path,header,output)
    
    empty = 0
    for ids in idlist:
        if len(ids) != 0:
            empty += len(ids)
    
    if empty == 0:
        return output
    
    db   = '/db/swissprot/uniprot_trembl.dat'
    path = '/db/swissprot/key/trembl/'
    idlist,output = retrieve(db,idfile,idlist,path,header,output)
    
    for ids in idlist:
        if len(ids) != 0:
            for e in ids:
                print ('%s does not exist' % (e), file=stderr)

    return output
#EOF
