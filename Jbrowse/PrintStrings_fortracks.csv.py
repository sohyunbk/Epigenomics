# Don't need to attach any conda envs
import glob, os
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=""
    )
    parser.add_argument(
        "-DirForFiles",
        "--DirForFiles",
        help="DirForFiles",
        required=True,
        dest="DirForFiles",
    )
    parser.add_argument(
            "-bw",
            "--bw",
            help="bw",
            required=False,
            dest="bw",
    )
    parser.add_argument(
            "-TiedName",
            "--TiedName",
            help="TiedName",
            required=False,
            dest="TiedName",
    )
    parser.add_argument(
            "-GitHubDir",
            "--GitHubDir",
            help="GitHubDir: atac",
            required=False,
            dest="GitHubDir",
    )
    parser.add_argument(
                "-bed",
                "--bed",
                help="bed",
                required=False,
                dest="bed",
    )
    parser.add_argument(
                "-GitHubSubDir",
                "--GitHubSubDir",
                help="GitHubSubDir",
                required=False,
                dest="GitHubSubDir",
    )
    args = vars(parser.parse_args())
    return parser

def BWFiles(Dir,BW,GithubDir,TiedName,SubDir):
    for sFiles in glob.glob(Dir+"/*.bw"):
        FileNamewithE = os.path.basename(sFiles)
        FileName = FileNamewithE.split(".")[0]
        #SampleName = FileName.split("_")[0]
        #CellName = "_".join(FileName.split("_")[1:]).split(".")[0]
        #atac,3_bif3_Re2_G2_M,3_bif3_Re2_G2_M,scATACMutantsV3/G2_M,,,scATACMutantsV3/3_bif3_Re2_G2_M.bw,,,,,,,
        cmd = "%s,%s,%s,%s,,,%s/%s,,,,,,,"%(GithubDir,FileName,FileName,TiedName,SubDir,FileNamewithE)
        print(cmd)

def bedFiles(Dir,Bed,GithubDir,TiedName):
    #peaks,HB122_WUS2_B73v5_Q30_default_finalBl.GEM,WUS2_peak,DAP_Sohyun,,,,,,,,,,
    for sFiles in glob.glob(Dir+"/*"):
        if sFiles.endswith(".bed"):
            FileNamewithE = os.path.basename(sFiles)
            FileName = FileNamewithE.split(".")[0]
        else:
            FileName = os.path.basename(sFiles)
        cmd = "%s,%s,%s,%s,,,,,,,,,,"%(GithubDir,FileName,FileName,TiedName)
        print(cmd)

if __name__ == "__main__":
    args = get_parser().parse_args()
    if args.bw == "yes":
        BWFiles(args.DirForFiles,args.bw,args.GitHubDir,args.TiedName,args.GitHubSubDir)

        #python //home/sb14489/Epigenomics/Jbrowse/PrintStrings_fortracks.csv.py -DirForFiles /scratch/sb14489/3.scATAC/4.Bif3Ref/5.Jbrowse_MACS2 -bw yes -TiedName scATAC -GitHubDir atacc
    elif args.bed == "yes":
        bedFiles(args.DirForFiles,args.bed,args.GitHubDir,args.TiedName)
        #python /home/sb14489/Epigenomics/Jbrowse/PrintStrings_fortracks.csv.py -DirForFiles /scratch/sb14489/3.scATAC/4.Bif3Ref/ -bed yes -TiedName InsertedSeq119+445bp -GitHubDir peak
