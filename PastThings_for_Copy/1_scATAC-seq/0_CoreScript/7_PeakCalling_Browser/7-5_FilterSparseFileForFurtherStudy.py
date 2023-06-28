#!/usr/bin/env python

import argparse

#python 7-5_FilterSparseFileForFurtherStudy.py -SparseBed  -MetaFile /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt -OutputName Temp

def get_parser():
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-SparseBed', "--SparseBed", help="SparseBed", required=True, dest='bed')
    parser.add_argument('-MetaFile', "--MetaFile", help="MetaFile", required=True, dest='meta')
    parser.add_argument('-OutputName', "--outputName", help="outputName", required=True, dest='output')
    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    MetaFile = open(args.meta,"r")
    #print(MetaFile)
    Dic = {}
    for sLine in MetaFile:
        Dic.setdefault(sLine.split("\t")[0])

    MetaFile.close()
    SparseFile = open(args.bed,"r")
    outfile=open(args.output,"w")
    for sLine in SparseFile:
        if sLine.split("\t")[3] in Dic.keys():
            outfile.write(sLine)

    SparseFile.close()
    outfile.close()
