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
        outfile = open(outfileName,"w")
        j = 0
        k =0
        Switch = "Off"
        for sLine in infile:
            j+=1
            #k+=1
            print(j)
            if j%4 == 1:
                Name=sLine.strip().split(" ")[0]
                if Name in List:
                    #print(Name)
                    Switch = "On"
                    outfile.write(sLine)
                    #print(sLine)
                else:
                    Switch = "Off"
            elif Switch == "On":
            #if j%4 !=1 and Switch = "On":
                outfile.write(sLine)
                #print(sLine)
        infile.close()
        outfile.close()
'''
            elif j%4 ==2:
                if Switch =="On":
                    outfile.write(sLine)
                    print(sLine)
            elif j%4 ==3:
                if Switch == "On":
                    outfile.write(sLine)
                    print(sLine)
            elif j%4 ==0:
                if Switch == "On":
                    outfile.write(sLine)
                    print(sLine)
                j =0
                Switch == "Off"'''

#########################################################################################
