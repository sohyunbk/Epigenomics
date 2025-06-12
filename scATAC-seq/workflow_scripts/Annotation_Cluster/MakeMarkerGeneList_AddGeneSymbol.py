################### Make Markergene


# 1) Read the geneSymbol - gene # IDEA:
def MakeSymbol_GeneIDDic(GeneInfo):
    infile = open(GeneInfo,"r")
    infile.readline()
    Dic = {}
    for sLine in infile:
        sList = sLine.replace(" ","").strip().split("\t")
        if len(sList) > 9 and len(sList[9]) >1 :
            #print(sList)
            Dic[sList[0]] = sList[9].upper()
            #print(sList[0]+"   "+ sList[9])
    infile.close()
    return(Dic)

## 2) Read bed file
# chr     start   end     geneID  name    type
def GeneIDSymbolMatch(BEDFile,SymbolDic,OutfilePath,OutfileName):
    infile = open(BEDFile,"r")
    outfile = open(OutfilePath+"/"+OutfileName,"w")
    outfile.write(infile.readline())
    for sLine in infile:
        sList = sLine.strip().split("\t")
        #print(sList[3])
        if sList[3] in SymbolDic:
            FinalLine= "\t".join(sList[0:4])+"\t"+SymbolDic[sList[3]]+"\t"+sList[5]+"\n"
        else:
            FinalLine= sLine
        outfile.write(FinalLine)
    infile.close()
    outfile.close()


def get_parser():
    parser = argparse.ArgumentParser(description='Add Gene Symbol.')
    parser.add_argument('-GeneInfo', "--GeneInfo", help="Zm00001eb.1.fulldata.txt", required=True, dest='GeneInfo')
    parser.add_argument('-MarkerGeneFile', "--MarkerGeneFile", help="MarkerGeneFile", required=True, dest='Marker')
    parser.add_argument('-OutPath', "--OutPath", help="OutPath", required=True, dest='OutPath')
    parser.add_argument('-OutFileName', "--OutFileName", help="OutFileName", required=True, dest='OutFileName')
    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    GeneInfo = "/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.txt"
    BEDFile = "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/231113_Top5DenovoGenesinA619_NoRedundant.txt"
    OutfilePath = "/scratch/sb14489/3.scATAC/0.Data/MarkerGene/"
    OutfileName = "231113_Top5DenovoGenesinA619_NoRedundant_withGeneSymbol.txt"
    SymbolDic = MakeSymbol_GeneIDDic(GeneInfo)
    GeneIDSymbolMatch(BEDFile,SymbolDic,OutfilePath,OutfileName)
#print(Dic)
