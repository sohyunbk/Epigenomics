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

def read_bed_file(bed_file):
    try:
        os.path.isfile(bed_file)
        bed_file_load = pybedtools.BedTool(bed_file)
        return bed_file_load
    except:
        FileNotFoundError
        pass

#One of the lines of intersect would be : chr1	610	1051	A619_G2_M.pool.macs_peak_172,A619_CalloseRelated.pool.macs_peak_173,A619_IM_SPM_SM.pool.macs_peak_183,A619_L1atFloralMeristem.pool.macs_peak_171,A619_PhloemPrecursor.pool.macs_peak_177,A619_XylemParenchyma_PithParenchyma.pool.macs_peak_183,A619_BundleSheath_VascularSchrenchyma.pool.macs_peak_174,A619_ProcambialMeristem_ProtoXylem_MetaXylem.pool.macs_peak_176	chr1	752	753	A619_G2_M.pool.macs_peak_172	1
def Get_OverlappedSummitUniqueCell(Intersect):
    Dic = {}
    for sLine in open(Intersect.fn):
        sList = sLine.strip().split("\t")
        PeakPos ="_".join(sList[0:3])
        #PeakName = sList[3]
        SummitName_list = sList[7].split(",") #	A619_G2_M.pool.macs_peak_172, ...
        Dic.setdefault(PeakPos,[])
        Dic[PeakPos] += SummitName_list
    return Dic

def Count_SummitNumberPerOneCelltype(Dic):
    NewDic ={}
    for sPos in Dic.keys():
        for Celltypes in Dic[sPos]: #Celltype = "A619_ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma.pool.macs_peak_70923"
            CellName = Celltypes.split(".")[0]
            NewDic
            List.append(CellName)

def GetPeaks_MoreThanOneSummit_PerOneCelltype(Dic):
    List_MorethanOneSummitPerOne = []
    for sPos in Dic.keys():
        Temp = {}
        for Celltypes in Dic[sPos]: #Celltype = "A619_ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma.pool.macs_peak_70923"
            CellName = Celltypes.split(".")[0]
            Temp.setdefault(CellName,0)
            Temp[CellName] +=1
        for sKeys in Temp.keys():
            if Temp[sKeys] > 1:
                if sPos not in List:
                    List.append(sPos)
    return(List)

def Calcuate_Distance_betweenSummits(Intersect,Peaks_withMoreThanOneSummit,OutfileName):
    Dic = {}
    for sLine in open(Intersect.fn):
        sList = sLine.strip().split("\t")
        PeakPos ="_".join(sList[0:3])
        nSummitStart = int(sList[5])
        if PeakPos not in Peaks_withMoreThanOneSummit:
            Dic.setdefault(PeakPos,[])
            Dic[PeakPos].append(nSummitStart)
    outfile = open(OutfileName,"w")
    for sPos in Dic:
        nDistance = max(Dic[sPos]) - min(Dic[sPos])
        outfile.write(sPos+'\t'+str(nDistance)+'\n')
    outfile.close()

def PeakLengthALl(AllBeds_list,Outfile):
    outfile = open(Outfile,"w")
    for Beds in AllBeds_list:
        for sLine in open(Beds.fn):
            nLen = int(sLine.strip().split("\t")[2])-int(sLine.strip().split("\t")[1])
            outfile.write("\t".join(sLine.strip().split("\t")[0:3])+"\t"+str(nLen)+"\n")
    outfile.close()

Dir = "c"
cellnames = [dirnames for root, dirnames, filenames in os.walk(Dir)][0]
SampleName = "A619"
AllBeds_list = [read_bed_file(Dir+"/"+cellname+"/"+SampleName+"_"+cellname+".reproducible_narrow_peaks_Onlychr") for cellname in cellnames]

PeakLengthALl(AllBeds_list,"/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_SplitPeaks/PeakLength")

MergedPeaks = AllBeds_list[0]
for sFiles in AllBeds_list[1:]: MergedPeaks = MergedPeaks.cat(sFiles,c=4, o="collapse")

#for sFiles in AllBeds_list[1:]: CombinedBed = CombinedBed.cat(*[sFiles],postmerge=False)
#print(a.cat(*[b],postmerge=False))

AllBeds_summit_list = [read_bed_file(Dir+"/"+cellname+"/"+SampleName+"_"+cellname+".reproducible_summits_Onlychr") for cellname in cellnames]
Merged_summit = AllBeds_summit_list[0]
for sFiles in AllBeds_summit_list[1:]: Merged_summit = Merged_summit.cat(sFiles,c=4, o="collapse")

Intersect_MergedPeaks_Summits = MergedPeaks.intersect(Merged_summit,wo=True)
IntersectDic = Get_OverlappedSummitUniqueCell(Intersect_MergedPeaks_Summits)
Peaks_List_withMoreThanOneSummit = GetPeaks_MoreThanOneSummit_PerOneCelltype(IntersectDic)
OutfileName="/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_SplitPeaks/Distance.txt"
Calcuate_Distance_betweenSummits(Intersect_MergedPeaks_Summits,Peaks_List_withMoreThanOneSummit,OutfileName)
## >>> len(Intersect_MergedPeaks_Summits) == len(Merged_summit)
