from selene_sdk.utils import load_path
from selene_sdk.utils import parse_configs_and_run
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('-wmlFile', "--wmlFile", help="It should be Model_wmlFile.", required=True, dest='wml')
    parser.add_argument('-learningRate', "--learningRate", help="learningRate.", required=True, dest='lr')
    parser.add_argument('-bedfile', "--bedfile", help="bedfile.", required=True, dest='bed')
    parser.add_argument('-feaeturefile', "--featurefile", help="featurefile.", required=True, dest='feature')
    parser.add_argument('-OutwmlfileName', "--OutwmlfileName", help="OutwmlfileName", required=True, dest='Outwml')
    parser.add_argument('-NewOutputDir', "--NewOutputDir", help="NewOutputDir", required=True, dest='OutDir')
    args = vars(parser.parse_args())
    return parser


def Modify_ymlFiles():
    infile = open(args.wml,"r")
    outfile = open(args.Outwml,"w")
    nFeatures = len(open(args.featurefile,"r").readlines())
    for sLine in infile:
        if "n_targets" in sLine:
            outfile.write(sLine.replace("n_targets:","n_targets: "+str(nFeatures))
        if "input_path" in sLine:
            outfile.write(sLine.replace("input_path:","input_path: "+str(args.feature))
        if "target_path:" in sLine:
            outfile.write(sLine.replace("target_path:","target_path: "+str(args.bed))
        if "output_dir:" in sLine:
            outfile.write(sLine.replace("output_dir:","output_dir: "+str(args.OutDir)))
        else:
            outfile.write(sLine)
    infile.close()
    outfile.close()

## Main Run!
args = get_parser().parse_args()
Modify_ymlFiles()

#configs = load_path(args.Outwml)
# deeperdeepsea was lr=0.01
#parse_configs_and_run(configs, lr=args.lr)
