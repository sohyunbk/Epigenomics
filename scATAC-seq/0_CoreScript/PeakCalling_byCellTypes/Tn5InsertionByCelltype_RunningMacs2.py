## This scrip is to call the peaks by cell types from scATAC-seq data
## This script can only be used when you already know the cell barcode - cell type after annotation step
## input data will be bed file
## Meta data with annotation of cell types is needed

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
        "-GenomeSize",
        "--GenomeSize",
        help="GenomeSize",
        required=True,
        dest="fai",
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

def RunMACS2(Outfile):
    for sBedFiles in glob.glob(Outfile+"*"):
        print(sBedFiles)
        Cmd = "macs2 callpeak -t %s -f BED --nomodel \
                    --keep-dup all --extsize 150 --shift -50 --qvalue .05 --bdg \
                    -n %s"%(sBedFiles,sBedFiles.replace(".bed",""))
        print(Cmd)
        os.system(Cmd)

def Normalize_bdg(Outfile,FaiFile):
    infile = open(FaiFile,"r")
    Fai = {}
    for sLine in infile:
        Fai[sLine.split("\t")[0]] = int(sLine.split("\t")[1])
    ## Make Depth Dic ###
    WD = "/".join(Outfile.split("/")[0:len(Outfile.split("/"))-1])+"/"
    FileList = []
    TotalReadDic = {}
    for FileName in glob.glob(WD+"*_treat_pileup.bdg"):
        FileList.append(FileName)
            #infile = open(FileName,"r")
    for d in glob.glob(WD+"*.bed"):
        #print(Dir.replace(".bed",""))
        TotalBed = open(d,"r")
        Length = len(TotalBed.readlines())
        #print(Length)
        TotalReadDic[d.replace(".bed","")] = Length
        TotalBed.close()
    ### Outfile the statistics of the reads numbers!
    StatOut = open(WD+"NumberofTn5_byReplicates_byCT.txt","w")
    print("DoneStats!")
    for i in TotalReadDic:
        StatOut.write(i+"\t"+str(TotalReadDic[i])+"\n")
    StatOut.close()

    for Files in FileList:
        infile = open(Files,"r")
        outfile = open(Files.replace(".bdg","_CPM.bdg"),"w")
        DicName = Files.replace("_treat_pileup.bdg","")
        print(DicName)
        for sLine in infile:
            sList = sLine.strip().split("\t")
            if int(sList[2]) < Fai[sList[0]]:
                nAbundance = float(sList[3])
                Normalized = (nAbundance/int(TotalReadDic[DicName]))*1000000
                outfile.write("\t".join(sList[0:3])+"\t"+str(Normalized)+"\n")
        infile.close()
        outfile.close()

if __name__ == "__main__":
    args = get_parser().parse_args()
    BedFile = args.bed
    MetaFile = args.m
    Outfile =args.Outfile
    FaiFile = args.fai

    Normalize_bdg(Outfile,FaiFile)
