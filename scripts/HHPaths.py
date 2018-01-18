#!/usr/bin/env python
docstring='''This is the python equivalent of HHPath.pm'''
import os

HHLIB=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ["HHLIB"]=HHLIB

bin_dict=dict(
    #### upstream hhsuite executables ####
    reformat=os.path.join(HHLIB,"scripts/reformat.pl"),
    hhblitsdb=os.path.join(HHLIB,"scripts/hhblitsdb.pl"),

    ### structure clustering by kClust ####
    kClust=os.path.join(HHLIB,"bin/kClust"),
    kClust_mkAln=os.path.join(HHLIB,"bin/kClust_mkAln"),
    clustalo=os.path.join(HHLIB,"bin/clustalo"),
)
