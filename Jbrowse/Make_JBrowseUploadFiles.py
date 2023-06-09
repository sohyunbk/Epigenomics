### This script has the whole functions for whatever input like bam or bdg or bed
### It should load "ml Anaconda3/2020.02" "activate /home/sb14489/.conda/envs/ucsc"

import argparse
import sys
import os
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
        description="Make JBrowseUpload File.\
        1: From bdg file to bw file: -Step bdgTobw -bdgFile {Path+Name} -Fai {chrFai} -OutputName {Path+NamePreFix} "
    )
    parser.add_argument(
        "-Step",
        "--Step",
        help="Step: bdgTobw or ",
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
    args = vars(parser.parse_args())
    return parser

def From_bdgfile_to_bwfile(BdgFile,OutFileName,Fai):
    Cmd_sort = "bedSort %s %s"%(BdgFile,BdgFile+"_Sorted")
    Cmd = "bedGraphToBigWig %s %s %s"%(BdgFile+"_Sorted",Fai, OutFileName+".bw")

if __name__ == "__main__":
    args = get_parser().parse_args()
    if args.Step == "bdgTobw":
        From_bdgfile_to_bwfile(args.bdgFile,args.OutputName,args.Fai)
