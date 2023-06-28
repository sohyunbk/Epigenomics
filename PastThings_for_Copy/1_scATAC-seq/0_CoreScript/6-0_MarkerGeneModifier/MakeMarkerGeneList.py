################### Make Markergene


# 1) Read the geneSymbol - gene # IDEA:
def MakeGeneID_SymbolDic(GeneInfo): ## Zm001eb000:er1
    infile = open(GeneInfo,"r") ## Zm001eb000:er1
    infile.readline()
    Dic = {}
    for sLine in infile:
        sList = sLine.replace(" ","").strip().split("\t")
        if len(sList) > 9 and len(sList[9]) >1 :
            #print(sList)
            Dic[sList[0]] = sList[9]
            #print(sList[0]+"   "+ sList[9])
    infile.close()
    return(Dic)

def MakeSymbol_GeneIDDic(GeneInfo):
    infile = open(GeneInfo,"r")
    infile.readline()
    Dic = {}
    for sLine in infile:
        sList = sLine.replace(" ","").strip().split("\t")
        if len(sList) > 9 and len(sList[9]) >1 :
            #print(sList)
            Dic[sList[9].upper()] = sList[0]
            #print(sList[0]+"   "+ sList[9])
    infile.close()
    return(Dic)

## 2) Read bed file
# chr     start   end     geneID  name    type
def BedCoordinateDic(BEDFile):
    infile = open(BEDFile,"r")
    Dic = {}
    for sLine in infile:
        sList = sLine.strip().split("\t")
        Dic[sList[3]] = "\t".join(sList[0:3])
    infile.close()
    return(Dic)
def ReadInput(Input):
    infile = open(Input,"r")
    LIST = []
    for sLine in infile:
        LIST.append(sLine.strip())
    return(LIST)

GeneInfo = "/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.txt"
BEDFile = "/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene.bed"
BedDic = BedCoordinateDic(BEDFile)
Path = "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/"
Switch = "GeneName"
Input = "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/CalloseGenelist.txt"
Inputlist = ReadInput(Input)
OutFileName = "CalloseGenelist.bed"
outfile = open(Path+"/"+OutFileName,"w")
outfile.write("chr\tstart\tend\tgeneID\tname\ttype\n")
if Switch == "GeneName":
    for i in Inputlist:
        SymbolDic = MakeSymbol_GeneIDDic(GeneInfo)
        if i.upper() in SymbolDic:
            GeneID = SymbolDic[i.upper()]
            if GeneID in BedDic:
                CoordinateInfo = BedDic[GeneID]
                outfile.write(CoordinateInfo+"\t"+GeneID+"\t"+i+"\tNA\n")

        else:
            print(i)
        #print(GeneID)

outfile.close()
#print(Dic)
