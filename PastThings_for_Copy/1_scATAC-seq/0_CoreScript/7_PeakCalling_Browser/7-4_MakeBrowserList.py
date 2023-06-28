import os, glob

BringFilesPath="/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/BwFiles"
## GitHubInfo
GithubDir="atac"
NewDir="scATACMutantsV3" ##

for sFiles in glob.glob(BringFilesPath+"/*.bw"):
    FileNamewithE = os.path.basename(sFiles)
    FileName = FileNamewithE.split(".")[0]
    SampleName = FileName.split("_")[0]
    CellName = "_".join(FileName.split("_")[1:])
    cmd = "%s,%s,%s,%s/%s,,,%s/%s,,,,,,,"%(GithubDir,FileName,FileName,NewDir,CellName,NewDir,FileNamewithE)
    print(cmd)
