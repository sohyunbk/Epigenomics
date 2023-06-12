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
    args = vars(parser.parse_args())
    return parser

def BWFiles(Dir,BW,GithubDir,TiedName):
    for sFiles in glob.glob(Dir+"/*.bw"):
        FileNamewithE = os.path.basename(sFiles)
        FileName = ".".join(FileNamewithE.split(".")[0])
        #SampleName = FileName.split("_")[0]
        #CellName = "_".join(FileName.split("_")[1:]).split(".")[0]
        #atac,3_bif3_Re2_G2_M,3_bif3_Re2_G2_M,scATACMutantsV3/G2_M,,,scATACMutantsV3/3_bif3_Re2_G2_M.bw,,,,,,,
        cmd = "%s,%s,%s,%s,,,%s,,,,,,,"%(GithubDir,FileName,FileName,TiedName,FileNamewithE)
        print(cmd)

if __name__ == "__main__":
    args = get_parser().parse_args()
    if args.bw != None:
        BWFiles(args.DirForFiles,args.bw,args.GitHubDir,args.TiedName)
