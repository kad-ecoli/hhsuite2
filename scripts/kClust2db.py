#!/usr/bin/env python
docstring='''
kClust2db.py db.fasta mydb/mydb
    cluster sequences in FASTA file db.fasta using kClust,
    and generate hhblits style database at mydb/mydb

Options:
    -tmpdir=/tmp/$USER/kClust_`date +%N`
        use -tmpdir as temperary folder

    -id=30
        kClust sequence identity cutoff 30%. legal values are: 
        20, 30, 40, 50, 60, 70, 80, 90, 99

    -ncpu=1
        number of CPU threads
'''
import sys,os
import subprocess
import shutil
from string import Template

from HHPaths import bin_dict

id2s_dict= { 20:0.52, 30:1.12, 40:1.73, 50:2.33,
    60:2.93, 70:3.53, 80:4.14, 90:4.74, 99:5.28}

kClust_template=Template("$kClust -i $infile -d $tmpdir/kClust")
kClust_mkAln_template=Template("$kClust_mkAln -c '$clustalo --threads=$ncpu -i $$infile -o $$outfile' -d $tmpdir/kClust --no-pseudo-headers|grep -P '^Filename:'|cut -d' ' -f2")
reformat_template=Template("$reformat fas a3m $filename $tmpdir/a3m/$basename.a3m")
hhblitsdb_template=Template("$hhblitsdb -cpu $ncpu -o $outdb -ia3m $tmpdir/a3m")

def mkdir_if_not_exist(tmpdir):
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

def make_tmpdir(tmpdir):
    ''' create tmp folder '''
    if not tmpdir:
        import random
        tmpdir="/tmp/%s/kClust_%s"%(
            os.getenv("USER"),random.randint(0,10**10))
        while(os.path.isdir(tmpdir)):
            tmpdir="/tmp/%s/kClust_%s"%(
                os.getenv("USER"),random.randint(0,10**10))
    mkdir_if_not_exist(tmpdir)
    sys.stdout.write("created folder %s\n"%tmpdir)
    return tmpdir 

def kClust2db(infile,outdb,tmpdir='.',s=1.12,ncpu=1):
    ''' cluster sequences in FASTA file "infile", and generate hhblits
    style database at outdb'''
    sys.stdout.write("#### cluster input fasta ####\n")
    cmd=kClust_template.substitute(dict(
        kClust=bin_dict["kClust"],
        infile=infile,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)

    sys.stdout.write("#### alignment within each cluster ####\n")
    cmd=kClust_mkAln_template.substitute(dict(
        kClust_mkAln=bin_dict["kClust_mkAln"],
        clustalo=bin_dict["clustalo"],
        ncpu=ncpu,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    stdout,stderr=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE).communicate()

    sys.stdout.write("#### reformat fas into a3m ####\n")
    a3mdir=os.path.join(tmpdir,"a3m")
    mkdir_if_not_exist(a3mdir)
    for filename in stdout.splitlines():
        os.system(reformat_template.substitute(dict(
            reformat=bin_dict["reformat"],
            filename=filename,
            tmpdir=tmpdir,
            basename=os.path.basename(os.path.splitext(filename)[0]),
        )))

    sys.stdout.write("#### build hhblitsdb ####\n")
    mkdir_if_not_exist(os.path.dirname(outdb))
    cmd=hhblitsdb_template.substitute(dict(
        hhblitsdb=bin_dict["hhblitsdb"],
        ncpu=ncpu,
        outdb=outdb,
        tmpdir=tmpdir,
    ))
    sys.stdout.write(cmd+'\n')
    os.system(cmd)
    return

if __name__=="__main__":
    seqID=30
    ncpu=1
    tmpdir=''
    
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-id="):
            seqID=float(arg[len("-id="):])
            if seqID<1:
                seqID=100*seqID
            seqID=int(seqID)
        elif arg.startswith("-ncpu="):
            ncpu=int(arg[len("-ncpu="):])
        elif arg.startswith("-tmpdir="):
            tmpdir=os.path.abspath(arg[len("-tmpdir="):])
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! No such option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if not seqID in id2s_dict:
        sys.stderr.write("ERROR! Illegal sequence identity cutoff %d\n"%seqID)
        exit()
    s=id2s_dict[seqID]

    if len(argv)!=2:
        sys.stderr.write(docstring)
        exit()

    infile=os.path.abspath(argv[0])
    outdb=os.path.abspath(argv[1])
    tmpdir=make_tmpdir(tmpdir)

    kClust2db(infile,outdb,tmpdir,s,ncpu)
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
