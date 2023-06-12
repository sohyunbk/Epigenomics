# Don't need to attach any conda envs
import glob
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
        cmd = "%s,%s,%s,%s,,,,,,,,,,"%(GithubDir,FileNamewithE,FileNamewithE,TiedName)
        print(cmd)

if __name__ == "__main__":
    args = get_parser().parse_args()
    if args.bw != None:
        BWFiles(args.DirForFiles,args.bw,args.GitHubDir,args.TiedName)
