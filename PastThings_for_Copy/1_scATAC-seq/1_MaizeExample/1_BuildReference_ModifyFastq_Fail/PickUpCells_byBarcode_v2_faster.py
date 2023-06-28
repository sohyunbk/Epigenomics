import sys, glob
import string

SampleName = sys.argv[1] #A619
##########################################################################################
## R2 file Make Dic
List = []
Switch="OFF"
for sFiles in glob.glob("/scratch/sb14489/3.scATAC/1.MaizeExample/1.Rawdata/OriginalData/Ex_"+SampleName+"_S*_R2_*"):
    print(sFiles)
    infile = open(sFiles,"r")
    j = 0
    for sLine in infile:
        j+=1
        if j%4 ==1:
            Name=sLine.strip().split(" ")[0]
            #print(Name)
            List.append(Name)
    infile.close()
print("Done Read R2")

for sFiles in glob.glob("/scratch/sb14489/3.scATAC/1.MaizeExample/1.Rawdata/OriginalData/"+SampleName+"_S*"):
    if "_R2_" not in sFiles:
    #print(sFiles)
        infile = open(sFiles,"r")
        outfileName = sFiles.replace(SampleName,"Ex_"+SampleName)
        print(outfileName)

        j = 0
        Dic = {}
        for sLine in infile:
            j+=1
            if j%4 == 1:
                Name=sLine.strip().split(" ")[0]
                Dic.setdefault(Name,"")
                Dic[Name]+=sLine
            else:
                Dic[Name]+=sLine
        infile.close()
        ############### Makde Dic #############
        outfile = open(outfileName,"w")
        print("Done MakeDic")
        for key in List:
            outfile.write(Dic[key])
        outfile.close()
#########################################################################################
