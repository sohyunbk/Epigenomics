### This script has the whole functions for whatever input like bam or bdg or bed
### It should load "ml Anaconda3/2020.02" "source activate /home/sb14489/.conda/envs/Jbrowse"
#ml Anaconda3/2020.02
#source activate /home/sb14489/.conda/envs/Jbrowse


import argparse
import sys
import os
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
        description="Make JBrowseUpload File.\
        1: From bdg file to bw file: \
            python /home/sb14489/Epigenomics/Jbrowse/Make_JBrowseUploadFiles.py \
            -Step bdgTobw -bdgFile {Path+Name} -Fai {chrFai} -OutputName {Path+NamePreFix} \
        2: python /home/sb14489/Epigenomics/Jbrowse/Make_JBrowseUploadFiles.py \
          -Step BedToTrack -bed /scratch/sb14489/3.scATAC/4.Bif3Ref/InsertedSeq.bed --OutputName InsertedSeq"
    )
    parser.add_argument(
        "-Step",
        "--Step",
        help="Step: bdgTobw or BedToTrack",
        required=False,
        dest="Step",
    )
    parser.add_argument(
        "-OutputName",
        "--OutputName",
        help="OutputName",
        required=False,
        dest="OutputName",
    )
    parser.add_argument(
        "-bdgFile",
        "--bdgFile",
        help="bdgFile",
        required=False,
        dest="bdgFile",
    )
    parser.add_argument(
            "-Fai",
            "--Fai",
            help="Fai",
            required=False,
            dest="Fai",
    )
    parser.add_argument(
            "-bed",
            "--bed",
            help="bed",
            required=False,
            dest="bed",
    )
    parser.add_argument(
            "-sam",
            "--sam",
            help="sam",
            required=False,
            dest="sam",
    )
    parser.add_argument(
            "-readlength",
            "--readlength",
            help="readlength",
            required=False,
            dest="readlength",
    )
    args = vars(parser.parse_args())
    return parser

def From_bdgfile_to_bwfile(BdgFile,OutFileName,Fai):
    Cmd_sort = "bedSort %s %s"%(BdgFile,BdgFile+"_Sorted")
    Cmd = "bedGraphToBigWig %s %s %s"%(BdgFile+"_Sorted",Fai, OutFileName+".bw")
    os.system(Cmd_sort)
    os.system(Cmd)

def From_bedfile_to_dirforTrack(BedFile,OutFileName):
    Path=os.path.dirname(BedFile)
    #/home/sb14489/jbrowse/bin/flatfile-to-json.pl --bed ./"$SampleName"/"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks --trackLabel "$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".reproducible_narrow_peaks --out ./
    Cmd = "/home/sb14489/jbrowse/bin/flatfile-to-json.pl --bed %s --trackLabel %s --out %s"%(BedFile,OutFileName,Path)
    os.system(Cmd)

def Make_bed_fromSamfile(Samfile,readlength,OutFileName):
    outfile = open(OutFileName,"w")
    Infile = open(Samfile,"r")
    for i in Infile:
        List = i.strip().split("\t")
        nFragment = abs(int(List[8]))
        if nFragment < int(readlength):
            nlength = nFragment
        else:
            nlength = int(readlength)
        nStart=int(List[3])
        sChr=List[2]
        outfile.write(sChr+"\t"+str(nStart)+"\t"+str(nStart+nlength)+"\n")
    outfile.close()
    Infile.close()

if __name__ == "__main__":
    args = get_parser().parse_args()
    if args.Step == "bdgTobw":
        From_bdgfile_to_bwfile(args.bdgFile,args.OutputName,args.Fai)
    if args.Step == "BedToTrack":
        From_bedfile_to_dirforTrack(args.bed,args.OutputName)
    if args.Step == "SamToBed":
        Make_bed_fromSamfile(args.sam,args.readlength,args.OutputName)
