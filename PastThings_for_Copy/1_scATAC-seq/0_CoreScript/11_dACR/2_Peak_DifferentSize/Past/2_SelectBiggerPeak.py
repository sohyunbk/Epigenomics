Path = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/"


## 2) CheckThrough All peaks
def MakeBedPosList(Infile,Dic):
    List = []
    infile = open(Infile,"r")
    for sLine in infile:
        sList = sLine.strip().split("\t")
        if "\t".join(sList[0:3]) not in Dic.keys():
            List.append("\t".join(sList[0:3]))
    infile.close()
    return(List)
def OverlappedSmallPeaks(Dic):
    List = []
    for i in Dic.keys():
        if len(Dic[i]) > 1:
            for s in Dic[i]:
                List.append(s)
    return(List)

def RunAll(CellTypeName):
    infile = open(Path+"./A619_bif3_For_dACR/"+CellTypeName+"_A619_bif3_Intersect.txt","r")

    Dic_Overlapped_WTKey = {}  # WTPos : [Mutant Pos]
    Dic_Overlapped_MutantKey = {}  # MutantPos : [WT_Pos]

    ## 1) Make OverlapDic!
    for sLine in infile:
        sList = sLine.strip().split("\t")
        WTPos="\t".join(sList[0:3])
        MutantPos="\t".join(sList[10:13])
        #print(MutantPos)
        Dic_Overlapped_WTKey.setdefault(WTPos,[])
        Dic_Overlapped_WTKey[WTPos].append(MutantPos)
        Dic_Overlapped_MutantKey.setdefault(MutantPos,[])
        Dic_Overlapped_MutantKey[MutantPos].append(WTPos)

    infile.close()
    #print(len(Dic_Overlapped_WTKey.keys()))
    #print(len(Dic_Overlapped_MutantKey.keys()))

    WTPeak_Extra =  MakeBedPosList(Path + "./A619/"+CellTypeName+"/A619_"+CellTypeName+".reproducible_narrow_peaks_Onlychr",Dic_Overlapped_WTKey)
    MutantPeak_Extra =  MakeBedPosList(Path + "./bif3/"+CellTypeName+"/bif3_"+CellTypeName+".reproducible_narrow_peaks_Onlychr",Dic_Overlapped_MutantKey)
    MutantPeak_small =OverlappedSmallPeaks(Dic_Overlapped_WTKey)
    WTPeak_small =OverlappedSmallPeaks(Dic_Overlapped_MutantKey)

    NewPeakList = []

    for sKey in Dic_Overlapped_WTKey.keys():
        if len(Dic_Overlapped_WTKey[sKey]) == 1:
            WT_Pos = sKey
            Mutant_Pos = Dic_Overlapped_WTKey[sKey][0]
            nWTLen = int(WT_Pos.split("\t")[2]) - int(WT_Pos.split("\t")[1])
            nMutantLen = int(Mutant_Pos.split("\t")[2]) - int(Mutant_Pos.split("\t")[1])
            if nWTLen > nMutantLen:
                NewPeakList.append(WT_Pos)
            else:
                NewPeakList.append(Mutant_Pos)

    AllPeaks = NewPeakList + WTPeak_Extra + MutantPeak_Extra + MutantPeak_small + WTPeak_small

    Temp = []
    outfile = open(Path+"./A619_bif3_For_dACR/"+CellTypeName+"_A619_bif3_LongerPeak.bed","w")
    for P in AllPeaks:
        if P not in Temp:
            outfile.write(P+"\n")
        Temp.append(P)
    print(len(Temp))
    outfile.close()

if __name__ == "__main__":
    Celltypes = ["L1","BundleSheath_VascularSchrenchyma","L1atFloralMeristem","CalloseRelated","PhloemPrecursor","FloralMeristem_SuppressedBract","ProcambialMeristem_ProtoXylem_MetaXylem","G2_M", "ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma","IM-OC","SPM-base_SM-base","IM_SPM_SM","XylemParenchyma_PithParenchyma"]
    for Cells in Celltypes:
        print(Cells)
        RunAll(Cells)
