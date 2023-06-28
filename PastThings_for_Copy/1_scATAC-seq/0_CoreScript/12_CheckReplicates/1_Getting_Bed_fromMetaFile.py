import sys

BedFile = sys.argv[1]
MetaFile = sys.argv[2]
Outfile = sys.argv[3]
print(sys.argv)

infile = open(MetaFile,"r")
infile.readline()
Dic ={}
CellTypeDic = {}
## To be easy, annotation column is the last col
for sLine in infile:
    sList =sLine.strip().split("\t")
    CellType = sList[len(sList)-1]
    Barcode = sList[0]
    #Dic.setdefault(CellType,[])
    #Dic[CellType].append(Barcode)
    Dic[Barcode]=CellType
    CellTypeDic.setdefault(CellType,"")
infile.close()

Tn5BedFile = open(BedFile,"r")
#outfile = open(Outfile,"w")
AllDic ={}
for sLine in Tn5BedFile:
    sList =sLine.strip().split("\t")
    # Barcode = sList[3]
    if sList[3] in Dic.keys():
        AssignedCell = Dic[sList[3]]
        AllDic.setdefault(AssignedCell,[])
        AllDic[AssignedCell].append(sLine)
Tn5BedFile.close()

#outfile.close()
for sCelltypes in AllDic.keys():
    outfile = open(Outfile+"_"+sCelltypes+".bed","w")
    for sNewLine in AllDic[sCelltypes]:
        outfile.write(sNewLine)
    outfile.close()
