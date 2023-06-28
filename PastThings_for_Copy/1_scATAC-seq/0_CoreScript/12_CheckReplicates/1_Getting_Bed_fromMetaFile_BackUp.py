import sys

BedFile = sys.argv[1]
MetaFile = sys.argv[2]
Outfile = sys.argv[3]
print(sys.argv)

infile = open(MetaFile,"r")
infile.readline()
List=[]
for sLine in infile:
    sList =sLine.strip().split("\t")
    List.append(sList[0])

infile.close()

Tn5BedFile = open(BedFile,"r")
outfile = open(Outfile,"w")
for sLine in Tn5BedFile:
    sList =sLine.strip().split("\t")
    if sList[3] in List:
        outfile.write(sLine)

Tn5BedFile.close()
outfile.close()
