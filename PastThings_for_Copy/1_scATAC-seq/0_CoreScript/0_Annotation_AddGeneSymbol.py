import glob, os, sys

infile = open("/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.txt","r")
BED = open("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed","r")
infile.readline()
Dic = {}
for sLine in infile:
    sList = sLine.replace(" ","").strip().split("\t")
    if len(sList) > 9 and len(sList[9]) >1 :
        #print(sList)
        Dic[sList[0]] = sList[9]
        #print(sList[0]+"   "+ sList[9])
infile.close()
NewBED = open("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr_AddGeneSymbol.bed","w")

for sBed in BED:
    sBed_list = sBed.strip().split("\t")

    #sGeneSymbol = Dic[sBed_list[3]]
    #NewLine = "\t".join(sBed_list[0:len(sBed_list)-1])
    if sBed_list[3] in Dic:
        sGeneSymbol = Dic[sBed_list[3]]
    else:
        sGeneSymbol = "DK"
    NewLine = "\t".join(sBed_list[0:len(sBed_list)-1])+"\t"+sGeneSymbol+"\n"
    NewBED.write(NewLine)
BED.close()
NewBED.close()
