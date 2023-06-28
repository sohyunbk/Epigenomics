def MakeDic(InfileName):
    infile = open(InfileName,"r")
    Dic = {}
    for sLine in infile:
        Dic.setdefault("\t".join(sLine.strip().split("\t")[0:3]),"")
    infile.close()
    return Dic
def MakeDic_byChr(InfileName):
    infile = open(InfileName,"r")
    Dic = {}
    for sLine in infile:
        Dic[sLine.strip().split("\t")[0]] = sLine.strip().split("\t")[1]
    infile.close()
    return Dic
def OutDic(Dic,OutputName):
    outfile = open(OutputName,"w")
    for sKey in Dic:
        outfile.write(sKey+"\n")
    outfile.close()

## 1) BlackList & Filter the end of the genome
def RunAll(CellType):
    WD = "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619_bif3_For_dACR/"
    CommonPeak =  MakeDic(WD+CellType+"_A619_bif3_LongerPeak_Sorted.bed")
    CommonPeak_BlackList =  MakeDic(WD+CellType+"_Overlap_ComA619Bif3_BlackList.bed")
    #print(CommonPeak_BlackList)
    GenomeSize =  MakeDic_byChr("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai")

    CommonPeakWithoutBlack = {}
    for Peaks in CommonPeak:
        if Peaks not in CommonPeak_BlackList.keys():
            Chr = Peaks.split("\t")[0]
            nEndPos = int(GenomeSize[Chr])
            if int(Peaks.split("\t")[2]) >nEndPos:
                print(Peaks)
            else:
                CommonPeakWithoutBlack.setdefault(Peaks,"")
    #print(CommonPeakWithoutBlack)
    print("total")
    print(len(CommonPeak.keys()))
    print("Remove BL")
    print(len(CommonPeakWithoutBlack.keys()))

    #### 2) Genic

    Dic_InsideGene={}
    Dic_StartOveralp_Morethan300bpOverlap = {}
    Dic_StartOveralp_Lessthan300bpOverlap = {}
    Dic_EndOverlap = {}

    infile = open(WD+CellType+"_Overlap_ComA619Bif3_Genes.bed","r")
    RepeatList = []
    for sLine in infile:
        sList = sLine.strip().split("\t")
        sPos = "\t".join(sLine.strip().split("\t")[0:3])
        sStrand = sList[8]
        FullOverlapLength = int(sList[2])-int(sList[1])
        if (int(sList[len(sList)-1])) == FullOverlapLength and (sPos in CommonPeakWithoutBlack) and (sPos not in RepeatList):
            #print("Check")
            Dic_InsideGene.setdefault(sPos,"")
        #elif int(sList[len(sList)-1]) > 400 and int(sList[len(sList)-1])  < 501 and sPos in CommonPeakWithoutBlack and sPos not in RepeatList:
        elif int(sList[len(sList)-1])  < FullOverlapLength and sPos in CommonPeakWithoutBlack and sPos not in RepeatList:
            if sStrand == "+" and int(sList[1]) < int(sList[5]): ## Start overlap
                if int(sList[len(sList)-1]) >= 300:
                    Dic_StartOveralp_Morethan300bpOverlap.setdefault(sPos,"")
                elif int(sList[len(sList)-1]) < 300:
                    Dic_StartOveralp_Lessthan300bpOverlap.setdefault(sPos,"")
            elif   sStrand == "+" and int(sList[1]) > int(sList[5]):  ## End part overlap
                Dic_EndOverlap.setdefault(sPos,"")
            elif   sStrand == "-" and int(sList[1]) > int(sList[5]):  ## Start overlap
                if int(sList[len(sList)-1]) >= 300:
                    Dic_StartOveralp_Morethan300bpOverlap.setdefault(sPos,"")
                elif int(sList[len(sList)-1]) < 300:
                    Dic_StartOveralp_Lessthan300bpOverlap.setdefault(sPos,"")
            elif   sStrand == "-" and int(sList[1]) < int(sList[5]):  ## End part overlap
                Dic_EndOverlap.setdefault(sPos,"")
        elif sPos in CommonPeakWithoutBlack and sPos not in RepeatList:
            Dic_StartOveralp_Lessthan400bpOverlap.setdefault(sPos,"")

        RepeatList.append(sPos)
    infile.close()

    #print("\t".join(sLine.split("\t")[0:4]))
    print("InsideGenes")
    print(len(Dic_InsideGene.keys()))
    print("Overlap_Morethan300 StartSite")
    print(len(Dic_StartOveralp_Morethan300bpOverlap.keys()))
    print("Overlap_Lessthan300 TSS")
    print(len(Dic_StartOveralp_Lessthan300bpOverlap.keys()))
    print("Overlap End")
    print(len(Dic_EndOverlap.keys()))

    Dic_Intergenic = {}
    for Peaks in CommonPeakWithoutBlack:
        if (Peaks not in Dic_InsideGene) and (Peaks not in Dic_StartOveralp_Morethan300bpOverlap) and (Peaks not in  Dic_StartOveralp_Lessthan300bpOverlap) and (Peaks not in  Dic_EndOverlap):
            Dic_Intergenic.setdefault(Peaks,"")
    print("Intergenic")
    print(len(Dic_Intergenic.keys()))

    OutputPath = WD
    OutDic(CommonPeakWithoutBlack,OutputPath+CellType+"_ComA619Bif3_BLRemove.bed")
    OutDic(Dic_InsideGene,OutputPath+CellType+"_ComA619Bif3_BLRemove_Genic.bed")
    OutDic(Dic_StartOveralp_Morethan300bpOverlap,OutputPath+CellType+"_ComA619Bif3_BLRemove_TSSStart_more300over.bed")
    OutDic(Dic_StartOveralp_Lessthan300bpOverlap,OutputPath+CellType+"_ComA619Bif3_BLRemove_TSSStart_less300over.bed")
    OutDic(Dic_EndOverlap,OutputPath+CellType+"_ComA619Bif3_BLRemove_TSSEnd.bed")
    OutDic(Dic_Intergenic,OutputPath+CellType+"_ComA619Bif3_BLRemove_Intergenic.bed")

if __name__ == "__main__":
    Celltypes = ["L1","BundleSheath_VascularSchrenchyma","L1atFloralMeristem","CalloseRelated","PhloemPrecursor","FloralMeristem_SuppressedBract","ProcambialMeristem_ProtoXylem_MetaXylem","G2_M", "ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma","IM-OC","SPM-base_SM-base","IM_SPM_SM","XylemParenchyma_PithParenchyma"]
    for Cells in Celltypes:
        print(Cells)
        RunAll(Cells)
