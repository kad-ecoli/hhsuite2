#!/usr/bin/env python
# This script is partially based on 
# https://github.com/soedinglab/uniclust-pipeline/blob/master/hhdatabase/reformat_old_cs219_ffindex.py
docstring='''
hhblitsdb3to2.py UniRef30_2020_01
    Convert hhblits3 database UniRef30_2020_01 to hhblits2 format.

    Note that it is not always possible to perfectly convert hhblits3
    format database, which can have unlimited sequence length, to
    hhblits2 format database, which has an upper limit of 32763 residues.

    Therefore, it is safer to search the converted database in two steps:
    firstly by hhblits2 and then, in case hhblits2 failed, by hhblits3
'''
import sys, os
import mmap
from collections import namedtuple
import textwrap

FFindexEntry = namedtuple("FFindexEntry", "name, offset, length")

AS62 ='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz*--'

AS90 ='!"#$%&\'()+,/0123456789:;<=?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~*--'

AS219='!"#$%&\'()+,/0123456789:;<=?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff*--'


def read_index(ffindex_filename):
    entries = []
    fh = open(ffindex_filename)
    for line in fh:
        tokens = line.split("\t")
        entries.append(FFindexEntry(tokens[0], int(tokens[1]), int(tokens[2])))
    fh.close()
    return entries

def read_data(ffdata_filename):
    fh = open(ffdata_filename, "r")
    data = mmap.mmap(fh.fileno(), 0, prot=mmap.PROT_READ)
    fh.close()
    return data

def read_entry_data(entry, data):
    return data[entry.offset:entry.offset + entry.length - 1]

def main(db):
    if not os.path.isfile(db+"_a3m.ffdata"):
        sys.stderr.write("ERROR! No such file %s_a3m.ffdata"%db)
        exit(1)

    if not os.path.isfile(db+"_a3m_db"):
        cmd="ln -s %s_a3m.ffdata %s_a3m_db"%(db,db)
        sys.stdout.write(cmd+'\n')
        os.system(cmd)
    if not os.path.isfile(db+"_hhm_db"):
        cmd="ln -s %s_hhm.ffdata %s_hhm_db"%(db,db)
        sys.stdout.write(cmd+'\n')
        os.system(cmd)

    if not os.path.isfile(db+"_hhm_db.index"):
        txt=''
        fp=open(db+"_hhm.ffindex",'r')
        for line in fp.read().splitlines():
            target,start,end=line.split()
            txt+="%s.a3m\t%s\t%s\n"%(target,start,end)
        fp.close()
        fp=open(db+"_hhm_db.index",'w')
        fp.write(txt)
        fp.close()

    if not os.path.isfile(db+"_a3m_db.index"):
        txt=''
        fp=open(db+"_a3m.ffindex",'r')
        for line in fp.read().splitlines():
            target,start,end=line.split()
            txt+="%s.a3m\t%s\t%s\n"%(target,start,end)
        fp.close()
        fp=open(db+"_a3m_db.index",'w')
        fp.write(txt)
        fp.close()

    if not os.path.isfile(db+".cs219") or \
       not os.path.isfile(db+".cs219.sizes"):

        input_data  = read_data(db+"_cs219.ffdata")
        input_index = read_index(db+"_cs219.ffindex")

        fh = open(db+".cs219", "wb")

        total_length = 0
        nr_sequences = len(input_index)
        sys.stdout.write("Converting %d sequences in total\n"%nr_sequences)
        line_break = bytearray("\n", "utf-8")[0]

        processed_seq=0
        for entry in input_index:
            entry_data = read_entry_data(entry, input_data)
            for i in range(len(entry_data)):
                if entry_data[i] == line_break:
                    entry_data = entry_data[(i+1):]
                    break
            total_length += len(entry_data)
            fh.write(bytearray(">"+entry.name+"\n", "utf-8"))
            #fh.write(entry_data)
            fh.write(textwrap.fill(
                ''.join([AS219[ord(e)] for e in entry_data]),100)+'\n')
            processed_seq+=1
            if processed_seq%10000==0:
                sys.stdout.write("Converting sequence %d\n"%processed_seq)
        fh.close()

        fh = open(db+".cs219.sizes", "w")
        fh.write(str(nr_sequences) + " " + str(total_length))
        fh.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=2:
        sys.stderr.write(docstring)
        exit(0)
    
    main(sys.argv[1])
