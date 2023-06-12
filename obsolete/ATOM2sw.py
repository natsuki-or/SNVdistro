#!/usr/bin/env python3.6

from sys import *
import datetime

RESDAT = "./residue.dat"

ATOM = ['_atom_site.group_PDB',
        '_atom_site.id',
        '_atom_site.type_symbol',
        '_atom_site.label_atom_id',
        '_atom_site.label_alt_id',
        '_atom_site.label_comp_id',
        '_atom_site.label_asym_id',
        '_atom_site.label_entity_id',
        '_atom_site.label_seq_id',
        '_atom_site.pdbx_PDB_ins_code',
        '_atom_site.Cartn_x',
        '_atom_site.Cartn_y',
        '_atom_site.Cartn_z',
        '_atom_site.occupancy',
        '_atom_site.B_iso_or_equiv',
        '_atom_site.pdbx_formal_charge',
        '_atom_site.auth_seq_id',
        '_atom_site.auth_comp_id',
        '_atom_site.auth_asym_id',
        '_atom_site.auth_atom_id',
        '_atom_site.pdbx_PDB_model_num']

#
#
def readmmcif(file):
    pdb = list()
    with open(file,"r") as f:
        for line in f:
            pdb.append(line.replace('\n',''))
            if "_entry.id" in line:
                pdbid = line.strip().split()[1]
 
    return pdbid, pdb

#
#
def ext_atom(pdb):
    atom = list()
    for i in range(len(pdb)):
        if ('ATOM'in pdb[i] or 'HETATM' in pdb[i]) and \
            (pdb[i].split()[ATOM.index('_atom_site.label_alt_id')] == '.' or \
             pdb[i].split()[ATOM.index('_atom_site.label_alt_id')] == 'A' or \
             pdb[i].split()[ATOM.index('_atom_site.label_alt_id')] == '1') and \
            pdb[i].split()[ATOM.index('_atom_site.pdbx_PDB_model_num')] == '1':
            atom.append(pdb[i].split())
    return(atom)

#
#
def newnum(atom):
    sernum = list()
    num = 1
    sernum.append(num)
    for i in range(1,len(atom)):
        if  atom[i][ATOM.index('_atom_site.label_asym_id')] != \
            atom[i-1][ATOM.index('_atom_site.label_asym_id')]:
            num = 1
        elif  atom[i][ATOM.index('_atom_site.label_seq_id')] != \
              atom[i-1][ATOM.index('_atom_site.label_seq_id')]:
            num += 1
        sernum.append(num)
    return sernum

#
#
def het2atom(atom,sernum):
    total = len(atom)
    st = -1
    alt = 0
    for i in range(total):
        if st == -1 and \
           atom[i][ATOM.index('_atom_site.group_PDB')] == 'HETATM':
            st = i-1
        if st != -1 and \
           atom[i][ATOM.index('_atom_site.group_PDB')] == 'ATOM':
            fn = i
            if atom[st][ATOM.index('_atom_site.label_asym_id')] == \
               atom[fn][ATOM.index('_atom_site.label_asym_id')] and \
               sernum[st] < sernum[fn] and \
               int(atom[st][ATOM.index('_atom_site.label_seq_id')]) <= \
               int(atom[fn][ATOM.index('_atom_site.label_seq_id')]):
                for j in range(st+1,fn):
#                   print("before: ",atom[j])
                    atom[j][ATOM.index('_atom_site.group_PDB')] = 'ATOM' 
#                   print("after:  ",atom[j])
                alt = 1
                print("REMARK   HETATM -> ATOM [%s]" % atom[st+1][ATOM.index('_atom_site.label_comp_id')],file=stderr)
            st = -1

    for i in range(total):
        if atom[i][ATOM.index('_atom_site.group_PDB')] == 'ATOM' and \
           sernum[i] != int(atom[i][ATOM.index('_atom_site.label_seq_id')]):
            print("REMARK   RENUMBERING",file=stderr)
            break
    return atom

#
#
def replace_res(atom,sernum):
    total = len(atom)
    for i in range(total):
        atom[i][ATOM.index('_atom_site.label_seq_id')] = sernum[i]
    return atom

#
#
def conv(atom):
    ot1 = list()
    ot3 = list()
    with open(RESDAT,"r") as f:
        for line in f:
            tmp = line.strip().split()
            ot3.append(tmp[0])
            ot1.append(tmp[1])

    chain = list()
    ch    = atom[0][ATOM.index('_atom_site.label_asym_id')]
    chain.append(atom[0][ATOM.index('_atom_site.label_asym_id')])

    seq   = list()
    sq    = ''

    total = len(atom)
    for i in range(total):
        if i != 0 and ch != atom[i][ATOM.index('_atom_site.label_asym_id')]:
            if len(sq) != 0: 
                seq.append(sq)
            else:
                chain.pop()
            sq = ''

            ch = atom[i][ATOM.index('_atom_site.label_asym_id')]
            chain.append(ch)
  
        if atom[i][ATOM.index('_atom_site.group_PDB')] == "ATOM" and \
           atom[i][ATOM.index('_atom_site.label_atom_id')] == "CA":
            sq += ot1[ot3.index(atom[i][ATOM.index('_atom_site.label_comp_id')])]
    if len(sq) != 0: 
        seq.append(sq)
    else:
        chain.pop()
    return(chain,seq)

#
#
def atom2seq(id):
    pdbid,pdb = readmmcif(id)
    atom = ext_atom(pdb)
    sernum = newnum(atom)
    atom = het2atom(atom,sernum)
    atom = replace_res(atom,sernum)
    chain,seq = conv(atom)
    return pdbid,chain,seq

#
# main
if __name__ == '__main__':

    pdbid,chain,seq = atom2seq('3tg0.cif')
#   for i in range(len(chain)):
#       print(pdbid,chain[i],len(seq[i]))
#       print(seq[i])

    for i in range(len(chain)):
        print("ID   %s%s          STANDARD;      PRT;  %4d AA." % \
                                        (pdbid,chain[i],len(seq[i])))
        print("DT   converted from PDB_ATOM %s(%s):%s" % \
                            (pdbid,chain[i],datetime.datetime.now()))
        print("CC   This is generated from PDB ATOM rows.       ")
        print("CC   The sequence is not always correct.         ")
        print("SQ   SEQUENCE  %4d AA;      MW;          CN;" % \
                                                       (len(seq[i])))
        for j in range(int((len(seq[i])-1)/60+1)):
            print("     ",end="")
            line = ""
            for k in range(60):
                if j*60+k >= len(seq[i]):
                    break
                line += seq[i][j*60+k]
                if (k+1) % 10 == 0:
                    line += ' '
            print(line)
        print("//")

#EOF
