#!/usr/bin/env python3.6

from sys import *
from pypdb import *

def pdbblast(seq):
    
    q = Query(seq,
              query_type="sequence", 
              return_type="polymer_entity")

    ID = list()
    iden = list()
    evalue = list()
    for result in q.search()['result_set']:
        ID.append(result['identifier'])
        iden.append(float(result['services'][0]['nodes'][0]['match_context'][0]['sequence_identity']))
        evalue.append(float(result['services'][0]['nodes'][0]['match_context'][0]['evalue']))

    return ID,iden,evalue


if __name__ == '__main__':
    if len(argv) != 2:
        exit ()

    id,iden,evalue = pdbblast(argv[1])
    l = len(id)
    for i in range(l):
        print(id[i], iden[i], evalue[i])

