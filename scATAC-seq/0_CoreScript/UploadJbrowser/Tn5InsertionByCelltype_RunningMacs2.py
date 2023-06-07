## This scrip is only used when you already know the cell barcode - cell type after annotation step
## input data will be bed file --> Tn5 insertion postion per barcode
## Meta data with annotation of cell types
## outfile Name
import argparse
import sys
import os
import pybedtools ##
import pandas as pd
import numpy
from multiprocessing import Pool, Manager
import multiprocessing
from functools import partial
import subprocess
import copy
import errno
import datetime
import random
import string
import glob

def get_parser():
    parser = argparse.ArgumentParser(
        description="Call Peaks for scATAC data. \
    Requires cluster annnotations, as well as BED file ipput."
    )
    parser.add_argument(
        "-BedFile",
        "--BedFile",
        help="BedFile",
        required=True,
        dest="bed",
    )
    parser.add_argument(
        "-MetaFile",
        "--MetaFile",
        help="MetaFile",
        required=True,
        dest="m",
    )
    parser.add_argument(
        "-Outfile",
        "--Outfile",
        help="Outfile",
        required=True,
        dest="Outfile",
    )
    args = vars(parser.parse_args())
    return parser

def MakeBarcode_CellTypeDic(MetaFile):
    ##### Start
    infile = open(MetaFile,"r")
    infile.readline()
    Dic ={}
    #CellTypeDic = {}
    ## To be easy, annotation column is the last col
    for sLine in infile:
        sList =sLine.strip().split("\t")
        CellType = sList[len(sList)-1]
        Barcode = sList[0]
        #Dic.setdefault(CellType,[])
        #Dic[CellType].append(Barcode)
        Dic[Barcode]=CellType
        #CellTypeDic.setdefault(CellType,"")
    infile.close()
    return Dic

def ReadTn5BedFile(BedFile,Dic):
    Tn5BedFile = open(BedFile,"r")
    #outfile = open(Outfile,"w")
    AllDic ={}
    for sLine in Tn5BedFile:
        sList =sLine.strip().split("\t")
        # Barcode = sList[3]
        if sList[3] in Dic.keys():
            AssignedCell = Dic[sList[3]]
            AllDic.setdefault(AssignedCell,[])
            AllDic[AssignedCell].append(sLine)
    Tn5BedFile.close()
    return AllDic

def WriteBedFiles(Outfile,AllDic):
    for sCelltypes in AllDic.keys():
        outfile = open(Outfile+"_"+sCelltypes+".bed","w")
        for sNewLine in AllDic[sCelltypes]:
            outfile.write(sNewLine)
        outfile.close()


if __name__ == "__main__":
    args = get_parser().parse_args()
    BedFile = args.bed
    MetaFile = args.m
    Outfile =args.Outfile


    #Dic = MakeBarcode_CellTypeDic(MetaFile) #Dic -- Barcode : Cell type.
    #AllDic = ReadTn5BedFile(BedFile,Dic)
    #WriteBedFiles(Outfile,AllDic)

    for sBedFiles in glob.glob(Outfile+"/*"):
        print(sBedFiles)
        Cmd = "macs2 callpeak -t %s -f BED --nomodel \
                    --keep-dup all --extsize 150 --shift -50 --qvalue .05 --outdir {final_output_dir_name} --bdg \
                    -n %s"%(sBedFiles,sBedFiles.replace("bed",""))
        print(Cmd)
        os.system(Cmd)
