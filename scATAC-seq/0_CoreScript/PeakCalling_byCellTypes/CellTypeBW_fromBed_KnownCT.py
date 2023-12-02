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


## This script is only for the genome browser not for the peak calling !
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
    parser.add_argument(
        "-Thread",
        "--Thread",
        help="Thread",
        required=True,
        dest="cores",
    )
    args = vars(parser.parse_args())
    return parser

### 1) First get bed file  by cell type
def process_bed_by_celltype(meta_file, bed_file, output_file_base):
    """
    Processes a BED file and writes separate BED files for each cell type.

    :param meta_file: Path to the metadata file.
    :param bed_file: Path to the BED file.
    :param output_file_base: Base path for the output files.
    """
    # Read metadata file and create a dictionary mapping barcodes to cell types
    list = []
    with open(meta_file, "r") as infile:
        infile.readline()  # Skip header
        dic = {}
        for sLine in infile:
            sList = sLine.strip().split("\t")
            cell_type = sList[-1]
            barcode = sList[0]
            dic[barcode] = cell_type

    # Process BED file and group data by cell type
    all_dic = {}
    with open(bed_file, "r") as tn5_bed_file:
        for sLine in tn5_bed_file:
            sList = sLine.strip().split("\t")
            if sList[3] in dic:
                assigned_cell = dic[sList[3]]
                all_dic.setdefault(assigned_cell, []).append(sLine)

    # Write output BED files for each cell type
    for cell_type in all_dic:
        with open(f"{output_file_base}_{cell_type}.bed", "w") as outfile:
            list.appned(f"{output_file_base}_{cell_type}.bed")
            for sNewLine in all_dic[cell_type]:
                outfile.write(sNewLine)

    return(list)


## 2) Run Macs2
def run_macs2_threaded(bed_files, output_directory, cores):
    with Pool(int(cores)) as pool:
            pool.map(
                partial(sub_func_macs2, output_dir=output_directory),
                bed_files,
            )


def  sub_func_macs2(bed_file, output_dir):
    output_file_name = bed_file.split("/")[-1].replace(".bed", ".macs")
    final_output_dir_name = output_dir if output_dir else "."

    generate_macs2_command = f"macs2 callpeak -t {bed_file} -f BED --nomodel \
    --keep-dup all --extsize 150 --shift -50 --qvalue .05 --outdir {final_output_dir_name} --bdg \
    -n {output_file_name}"

    print(f"Running MACS2 Command: {generate_macs2_command}")
    subprocess.run(generate_macs2_command, shell=True, check=True)

    print("Done Running MACS2 Calls")


if __name__ == "__main__":
    args = get_parser().parse_args()
    BedFile = args.bed
    MetaFile = args.m
    Outfile =args.Outfile
    ### 1) First get bed file  by cell type
    bed_files = process_bed_by_celltype(args.m, args.bed, args.outfile)
    print(bed_files)
    run_macs2_threaded(bed_files, args.Outfile, args.cores)
    ## 2) Run Macs2
