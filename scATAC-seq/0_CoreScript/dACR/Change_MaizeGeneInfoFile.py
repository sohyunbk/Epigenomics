infile = open("/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.txt","r")
outfile =  open("/scratch/sb14489/0.Reference/Maize_B73/Zm00001eb.1.fulldata.Curated.txt","w")
FirstLine = infile.readline().replace(' ','').split("\t")
ColNum = len(FirstLine)
outfile.write("\t".join(FirstLine))
for sLine in infile:
    sList = sLine.replace(" ","").strip().split("\t")
    if len(sList) > ColNum:
        print("Error")
    else:
        for j in range(0,ColNum):
            if j < len(sList):
                if len(sList[j])<1:
                    outfile.write("NA")
                else:
                    outfile.write(sList[j])
                if j == (ColNum-1):
                    outfile.write("\n")
                else:
                    outfile.write("\t")
            else:
                outfile.write("NA")
                if j == (ColNum-1):
                    outfile.write("\n")
                else:
                    outfile.write("\t")
outfile.close()
