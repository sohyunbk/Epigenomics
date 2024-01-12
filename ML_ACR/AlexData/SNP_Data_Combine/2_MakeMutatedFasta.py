#conda activate pytorch
from pyfaidx import Fasta
import argparse
import sys
import os, glob

ef get_parser():
    parser = argparse.ArgumentParser(
        description="Call Peaks for scATAC data. \
    Requires cluster annnotations, as well as BED file ipput."
    )
    parser.add_argument(
        "-FastaFile",
        "--FastaFile",
        help="FastaFile",
        required=True,
        dest="Fasta",
    )
    parser.add_argument(
        "-SNPFile",
        "--SNPFile",
        help="SNPFile",
        required=True,
        dest="SNPFile",
    )
    args = vars(parser.parse_args())
    return parser


args = get_parser().parse_args()

#'/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/control_SNVs_curated.txt'
#'/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/ControlSNPChange_MaizeV5.fa'
with open(args.SNPFile) as mut_table:
    with Fasta(args.Fasta, mutable=True) as fasta:
        next(mut_table)  # Skip the header line
        for line in mut_table:
            _, snpID, ref, alt, *_ = line.rstrip().split()
            rname, pos = snpID.split('_')[0], snpID.split('_')[1]
            #print(rname +"\t"+str(pos)+"\t"+alt)
            #chr9	160483840	G
            if fasta[rname][int(pos) - 1]!= ref.upper():
                print(f"Mismatch at {rname} position {pos}: expected {ref}, found {fasta[rname][int(pos) - 1]}")
            else:
                fasta[rname][int(pos) - 1] = alt
