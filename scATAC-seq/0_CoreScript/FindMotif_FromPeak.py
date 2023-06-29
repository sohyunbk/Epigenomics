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
    

if __name__ == "__main__":
    args = get_parser().parse_args()
    BedFile = args.bed
    MetaFile = args.m
    Outfile =args.Outfile
    FaiFile = args.fai


    Normalize_bdg(Outfile,FaiFile)
