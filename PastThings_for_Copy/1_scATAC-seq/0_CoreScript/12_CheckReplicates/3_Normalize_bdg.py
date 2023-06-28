## Normalize with CPM: (Amount / library depth) * 1000,000
import os, glob
## Make Fai Dic
infile = open("/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai","r")
Fai = {}
for sLine in infile:
    Fai[sLine.split("\t")[0]] = int(sLine.split("\t")[1])
## Make Depth Dic ###
WD = "/scratch/sb14489/3.scATAC/2.Maize_ear/12.Macs_GenomeBrowserByReplicates/"
FileList = []
TotalReadDic = {}
for Dir in os.listdir(WD):
    d = os.path.join(WD, Dir)
    if os.path.isdir(d) and "Re" in Dir:
        #print(d)
        FileName = d+"/"+Dir+"_treat_pileup.bdg"
        FileList.append(FileName)
        #infile = open(FileName,"r")
    else :
        if "Re" in Dir:
            #print(Dir.replace(".bed",""))
            TotalBed = open(d,"r")
            Length = len(TotalBed.readlines())
            #print(Length)
            TotalReadDic[Dir.replace(".bed","")] = Length
            TotalBed.close()

### Outfile the statistics of the reads numbers!
StatOut = open(WD+"NumberofTn5_byReplicates_byCT.txt","w")
print("DoneStats!")
for i in TotalReadDic:
    StatOut.write(i+"\t"+str(TotalReadDic[i])+"\n")
StatOut.close()

for Files in FileList:
    infile = open(Files,"r")
    outfile = open(Files.replace(".bdg","_CPM.bdg"),"w")
    DicName = Files.split("/")[len(Files.split("/"))-1].replace("_treat_pileup.bdg","")
    print(DicName)
    for sLine in infile:

        sList = sLine.strip().split("\t")
        if int(sList[2]) < Fai[sList[0]]:
            nAbundance = float(sList[3])
            Normalized = (nAbundance/int(TotalReadDic[DicName]))*1000000
            outfile.write("\t".join(sList[0:3])+"\t"+str(Normalized)+"\n")
    infile.close()
    outfile.close()
