#ml Anaconda3/2020.02
#source activate r_env
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
from pybedtools import BedTool
import pybedtools

def Make_Fai_Dic(GenomeSizeFile):
    infile = open(GenomeSizeFile,"r")
    Dic = {}
    for sLine in infile:
        if sLine.startswith("chr"):
            sList = sLine.strip().split("\t")
            Dic[sList[0]] = int(sList[1])
    return Dic
def ReadBedFiles_onlyChr(InfileName):
    infile = open(InfileName,"r")
    List=[]
    for sLine in infile:
        sList = sLine.strip().split("\t")
        if sList[0] in GenomeSize_Dic.keys() and int(sList[2]) < GenomeSize_Dic[sList[0]]:
            List.append(sLine)
    infile.close()
    BedtoolRead = BedTool(List)
    return BedtoolRead

def Merge_TwoBedListSameOrder(List1,List2,Cells,OutputPath,OutFileSampleName):
    List = []
    for i in range(0,len(List1)):
        OutFileName = OutputPath+"/"+OutFileSampleName+"_"+Cells[i]+"_MergedPeak.bed"
        #print(OutFileName)
        MergedPeak = List1[i].cat(List2[i],c=4, o="collapse").moveto(OutFileName)
        #MergedPeak.moveto(OutFileName)
        #print(Cells[i])
        #print(len(MergedPeak))
        List.append(MergedPeak)
    return(List)

def Make_Fake_Summit(MergedPeaks): ## As peaks are merged by sample, it's hard to choose the summit.#So I just made the summit in the middle of peak.
    NewList = []
    for i in open(MergedPeaks.fn):
        sList = i.strip().split("\t")
        nPeakLength = int(sList[2]) - int(sList[1])
        nNewStart = int(sList[1])+((nPeakLength//2)-1)
        nNewEnd = nNewStart+1
        NewList.append(sList[0]+"\t"+str(nNewStart)+"\t"+str(nNewEnd)+"\t"+"]]".join(sList)+"\n")
    BedtoolRead = BedTool(NewList)
    return BedtoolRead

def Outfile_forthe4thColumn(BedToolsInput,OutfileName):
    outfile = open(OutfileName,"w")
    for sLine in open(BedToolsInput.fn):
        OutLine = "\t".join(sLine.strip().split("\t")[3].split("]]"))+"\n"
        outfile.write(OutLine)
    outfile.close()

def Outfile_GenicIntergenic(k, Summit, AnnBed,cellnames,OutPutDir,OutfileName):
     SummitFile = Summit[k]
     CommonPeak_Intergenic = SummitFile-AnnBed
     CommonPeak_genic = SummitFile + AnnBed
     print(cellnames[k])
     print(len(SummitFile))
     print(len(CommonPeak_Intergenic))
     print(len(CommonPeak_genic))
     IntergenicOutname= OutPutDir+"/"+OutfileName+"_"+cellnames[k]+"_MergedPeak_Intergenic.bed"
     GenicOutname= OutPutDir+"/"+OutfileName+"_"+cellnames[k]+"_MergedPeak_Genic.bed"
     Outfile_forthe4thColumn(CommonPeak_Intergenic,IntergenicOutname)
     Outfile_forthe4thColumn(CommonPeak_genic,GenicOutname)

## Functions for Method2 ###

def Make_peak_with_the_SizeYouWant(MergedPeaks,NewPeakSize):
    NewPeaks =[]
    for Peak in MergedPeaks:
        nPeakLength = int(Peak.end) - int(Peak.start)
        # 1) Peaks < Size you want --> going to extend their edgeList
        if nPeakLength < NewPeakSize:
            LengthYouneed1 = (NewPeakSize-nPeakLength)//2
            LengthYouneed2 = (NewPeakSize-nPeakLength) - LengthYouneed1
            NewStart = int(Peak.start) - LengthYouneed1
            NewEnd = int(Peak.end) + LengthYouneed2
            NewPeaks.append(Peak.chrom+"\t"+str(NewStart)+"\t"+str(NewEnd)+"\t"+Peak.name+"\n")
        elif nPeakLength == NewPeakSize:
            NewStart = int(Peak.start)
            NewEnd = int(Peak.end)
            NewPeaks.append(Peak.chrom+"\t"+str(NewStart)+"\t"+str(NewEnd)+"\t"+Peak.name+"\n")
        elif nPeakLength > NewPeakSize:
            ## 1) When there is no remain
            if nPeakLength%NewPeakSize == 0:
                nDevide = nPeakLength//NewPeakSize
                NewStart = int(Peak.start)
                for s in range(0,nDevide):
                    #print(s)
                    NewEnd = NewStart + NewPeakSize
                    NewPeaks.append(Peak.chrom+"\t"+str(NewStart)+"\t"+str(NewEnd)+"\t"+Peak.name+"\n")
                    NewStart=NewEnd
            ## 2) Whre there is remain. --> 2-1)extend the seq 2-2)Divide
            elif nPeakLength%NewPeakSize != 0:
                SeqLengthSouldbeAdded2 = (NewPeakSize - (nPeakLength%NewPeakSize))//2
                SeqLengthSouldbeAdded1 = (NewPeakSize - (nPeakLength%NewPeakSize)) - SeqLengthSouldbeAdded2
                #print(SeqLengthSouldbeAdded1)
                ExtendStart = int(Peak.start) - SeqLengthSouldbeAdded1
                #ExtendEnd = int(Peak.end) + SeqLengthSouldbeAdded2
                nDevide = (SeqLengthSouldbeAdded1+nPeakLength+SeqLengthSouldbeAdded2)//NewPeakSize
                for s in range(0,nDevide):
                    ExtendEnd= ExtendStart + NewPeakSize
                    NewPeaks.append(Peak.chrom+"\t"+str(ExtendStart)+"\t"+str(ExtendEnd)+"\t"+Peak.name+"\n")
                    ExtendStart=ExtendEnd
    BedtoolRead = BedTool(NewPeaks)
    return BedtoolRead
    #print(a.merge(d=-1))
## Inputs! ###
def Check_Length_andFilterPeaksOutsideChromosome(feature):
    if(len(feature)==NewPeakSize) and \
     feature.start < GenomeSize_Dic[feature.chrom] and \
     feature.end < GenomeSize_Dic[feature.chrom]:
        return feature

def get_parser():
    parser = argparse.ArgumentParser(
        description="Cut peaks to fix humped peaks"
    )
    parser.add_argument(
        "-method",
        "--method",
        help="Method1: Requires two InputPath. Just Merge the peaks. use: MergeCutPeaks_toFixHumpedPeaks.py -method Method2 \
        -input1 Path1 -input2 Path2 -fai Chromosome faifile -OutputFileName OutputFileName -OutputPath OutputPath \
        Method2: Needs only one InputPath. ExtendPeaks Cut with the size you set up",
        required=True,
        dest="method",
    )
    parser.add_argument(
        "-inputpath1",
        "--inputpath1",
        help="inputpath1",
        required=False,
        dest="Dir1",
    )
    parser.add_argument(
        "-inputpath2",
        "--inputpath2",
        help="inputpath2",
        required=False,
        dest="Dir2",
    )
    parser.add_argument(
        "-fai",
        "--fai",
        help="Chromosome faifile",
        required=False,
        dest="fai",
    )
    parser.add_argument(
        "-OutputFileName",
        "--OutputFileName",
        help="OutputFileName",
        required=False,
        dest="SampleName",
    )
    parser.add_argument(
        "-Ann",
        "--Ann",
        help="Ann bed file to get intergenic and genic peak!",
        required=False,
        dest="Ann",
    )
    parser.add_argument(
        "-OutputPath",
        "--OutputPath",
        help="OutputPath",
        required=False,
        dest="OutputDir",
    )
    parser.add_argument(
            "-Sparse1",
            "--Sparse1",
            help="Sparse1",
            required=False,
            dest="Sparse1",
    )
    parser.add_argument(
            "-Sparse2",
            "--Sparse2",
            help="Sparse2",
            required=False,
            dest="Sparse2",
    )
    ##### Method 2 arg
    parser.add_argument(
        "-inputpath",
        "--inputpath",
        help="inputpath",
        required=False,
        dest="Dir",
    )
    parser.add_argument(
        "-size",
        "--size",
        help="Peak size you want should be larger than 150bp at least",
        required=False,
        dest="size",
    )
    args = vars(parser.parse_args())
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    GenomeSize_Dic = Make_Fai_Dic(args.fai)
    #============= 1st option: Easiest way! Just merge the peaks - with different sizes by cell type ============#
    if args.method == "Method1":
        cellnames = [dirnames for root, dirnames, filenames in os.walk(args.Dir1)][0]
        print(cellnames)
        InputSample1 = args.Dir1.split("/")[len(args.Dir1.split("/"))-1]
        InputSample2 = args.Dir2.split("/")[len(args.Dir2.split("/"))-1]
        AllBeds_list1 = [ReadBedFiles_onlyChr(args.Dir1+"/"+cellname+"/"+InputSample1+"_"+cellname+".reproducible_summits.passing_FDR") for cellname in cellnames]
        AllBeds_list2 = [ReadBedFiles_onlyChr(args.Dir2+"/"+cellname+"/"+InputSample2+"_"+cellname+".reproducible_summits.passing_FDR") for cellname in cellnames]
        MergedPeakList = Merge_TwoBedListSameOrder(AllBeds_list1,AllBeds_list2,cellnames,args.OutputDir,args.SampleName)
        FakeSummitList = [Make_Fake_Summit(MergedPeakSet) for MergedPeakSet in MergedPeakList]
        for k in range(0,len(FakeSummitList)): Outfile_GenicIntergenic(k,FakeSummitList,args.Ann,cellnames,args.OutputDir,args.SampleName)
        ## Give up getting sparse ## Should make chr9_144857976_144858376	CB:Z:AAACGAAAGAGCGAAA-1_A619_2	2
        if args.Sparse1 == None:
            pass
        else:
            SPARSE1 = BedTool(args.Sparse1)
            SPARSE1_Chr = SPARSE1.filter(lambda b: b.chrom.startswith("chr"))
            SPARSE2 = BedTool(args.Sparse2).filter(lambda b: b.chrom.startswith("chr"))
            SPARSE2_Chr = SPARSE2.filter(lambda b: b.chrom.startswith("chr"))
            #MergedPeaks = MergedPeakList[0]
            #for sFiles in MergedPeakList[1:]: MergedPeaks = MergedPeaks.cat(sFiles)
            ## Should make chr9_144857976_144858376	CB:Z:AAACGAAAGAGCGAAA-1_A619_2	2
            #Temp = MergedPeaks.intersect(SPARSE1_Chr,wo=True) # chr1	610	1051 chr 614 615 CB:Z:AAACGAAAGAGCGAAA-1_A619_2
            #print(Temp[0])
            #print(Temp[1])
            #Temp = (SPARSE1_Chr + MergedPeaks).saveas(args.OutputDir+"/A619_toAllMergedPeak.sparse")
            #Temp = (SPARSE2_Chr + MergedPeaks).saveas(args.OutputDir+"/Bif3_toAllMergedPeak.sparse")


    #============== 2nd option: Just cut the peak with the size you want!
    #============== Peak size will be all same but can not guarantee not to cut some important motifs :(
    elif args.method == "Method2":
        os.system("mkdir "+args.OutputDir)
        cellnames = [dirnames for root, dirnames, filenames in os.walk(args.Dir)][0]
        print(cellnames)
        InputSample = args.Dir.split("/")[len(args.Dir.split("/"))-1]
        AllBeds_list = [ReadBedFiles_onlyChr(args.Dir+"/"+cellname+"/"+InputSample+"_"+cellname+".reproducible_summits.passing_FDR") for cellname in cellnames]
        MergedPeaks = AllBeds_list[0]
        for sFiles in AllBeds_list[1:]: MergedPeaks = MergedPeaks.cat(sFiles,c=4, o="collapse")
        NewPeakSize = int(args.size)
        Peaks_Cut = Make_peak_with_the_SizeYouWant(MergedPeaks,NewPeakSize)
        print(len(Peaks_Cut))
        ##Check
        for Peaks_Split in Peaks_Cut:
            if len(Peaks_Split) != NewPeakSize: print("Error")
        print("Check the length is done")
        ## Filter Some peaks that are overlapped each other because of extension
        Peaks_Cut_Merged = Peaks_Cut.merge(d=-1,c=4, o="collapse")
        OutfileName = args.OutputDir+"/"+args.SampleName+"_CutPeaksWith"+str(NewPeakSize)+"bp.bed"
        Final_Peaks_Cut_Merged= Peaks_Cut_Merged.each(Check_Length_andFilterPeaksOutsideChromosome).saveas(OutfileName)
        print("TheFinal Numberof peak:%s"%(len(Final_Peaks_Cut_Merged)))

    #============= 3rd option:
    elif args.method == "Method3":
        print("Not yet")
