import glob, os, sys

infile = open("/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.txt","r")
Path = "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/AnnV3/"
infile.readline()
Dic = {}
for sLine in infile:
    sList = sLine.replace(" ","").strip().split("\t")
    if len(sList) > 9 and len(sList[9]) >1 :
        #print(sList)
        Dic[sList[0]] = sList[9]
        #print(sList[0]+"   "+ sList[9])
infile.close()

for sFiles in glob.glob(Path+"*_de_novo.visualize.bed"):
    print(sFiles)
    infile = open(sFiles,"r")
    outfile = open(sFiles.replace(".bed",".GeneSymbol.bed"),"w")
    i = 0
    infile.readline()
    outfile.write("chr\tstart\tend\tgeneID\tname\ttype\n")
    for sLine in infile:
        sList = sLine.strip().split("\t")
        i+=1
        #print(sList)
        if sList[3] in Dic.keys():
            NewLine = "\t".join(sList[0:4])+"\t"+Dic[sList[3]]+"\tDenovo\n"
        else:
            NewLine = "\t".join(sList[0:4])+"\t"+sList[3]+"\tDenovo\n"
        #if i < 100:
        outfile.write(NewLine)
    infile.close()
    outfile.close()
