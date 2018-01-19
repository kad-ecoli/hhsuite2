#!/usr/bin/env python
# This is the python equivalent of HHPath.pm
import os,sys

HHLIB=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ["HHLIB"]=HHLIB

bin_dict=dict(
    #### upstream hhsuite executables ####
    hhblits       =os.path.join(HHLIB,"bin/hhblits"),
    cstranslate   =os.path.join(HHLIB,"bin/cstranslate"),
    hhfilter      =os.path.join(HHLIB,"bin/hhfilter"),
    reformat      =os.path.join(HHLIB,"scripts/reformat.pl"),
    hhblitsdb     =os.path.join(HHLIB,"scripts/hhblitsdb.pl"),

    #### qhmmer ####
    qhmmsearch    =os.path.join(HHLIB,"bin/qhmmsearch"),
    qhmmbuild     =os.path.join(HHLIB,"bin/qhmmbuild"),
    qjackhmmer    =os.path.join(HHLIB,"bin/qjackhmmer"),
    eslsfetch     =os.path.join(HHLIB,"bin/esl-sfetch"),

    #### MSAParser ####
    fasta2aln     =os.path.join(HHLIB,"bin/fasta2aln"),
    fastaCov      =os.path.join(HHLIB,"bin/fastaCov"),
    realignMSA    =os.path.join(HHLIB,"bin/realignMSA"),
    rmRedundantSeq=os.path.join(HHLIB,"bin/rmRedundantSeq"),
    calNf         =os.path.join(HHLIB,"bin/calNf"),

    ### sequence clustering by kClust ####
    kClust        =os.path.join(HHLIB,"bin/kClust"),
    kClust_mkAln  =os.path.join(HHLIB,"bin/kClust_mkAln"),
    clustalo      =os.path.join(HHLIB,"bin/clustalo"),
)

def check_hhsuite_binaries():
    ''' check if all binaries listed in bin_dict are present '''
    for exe in bin_dict:
        if not os.path.isfile(bin_dict[exe]):
            sys.stderr.write("ERROR! Cannot locate %s at %s\n"%(
                exe,bin_dict[exe]))
    return

if __name__=="__main__":
    check_hhsuite_binaries()
