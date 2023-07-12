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
import random

MappableRegionsFile = "/scratch/sb14489/8.ML_ACR/1.InputBed/B73_v5_sim.mq30.merge_within_75bp.mappable_Sorted_OnlyChr.bed"
AllACRBedFile = "/scratch/sb14489/8.ML_ACR/1.InputBed/sorted_Seedling_Peaks.bed"
OutfileName = "/scratch/sb14489/8.ML_ACR/1.InputBed/B73_v5_sim.mq30.merge_within_75bp.mappable_Sorted_OnlyChr_SameNumberNeg.bed"
MappableRegions = BedTool(MappableRegionsFile)
AllACRBed = BedTool(AllACRBedFile)

##CheckNumber: 82098
#Merged = AllACRBed.merge(d=-1,c=4, o="collapse")

NotACRRegions = MappableRegions-AllACRBed
ACRRegions = MappableRegions-NotACRRegions

#len(MappableRegions)
#len(NotACRRegions)
random.seed(5)
RandomnumberList=random.sample(range(0,len(NotACRRegions)),len(ACRRegions))
#NotACRRegions.saveas(OutfileName)
Outfile= open(OutfileName,"w")
NewList = []
for s in ACRRegions:
    Outfile.write(str(s))
for j in RandomnumberList:
    Outfile.write(str(NotACRRegions[j]))
    #NewList.append(NotACRRegions[j])
#NewList.saveas(OutfileName)

#for i in NewList:
    #Outfile.write(i)
Outfile.close()
