from selene_sdk.utils import load_path
from selene_sdk.utils import parse_configs_and_run
import argparse

'''
/home/sb14489/miniconda3/envs/pytorch/bin/python Run_Selene_MakeYMLFile.py \
 --wmlFile /home/sb14489/Epigenomics/ML_ACR/Selene_AllModels_BashRun_ymlFiles/Standard_DanQ_WithoutCuda_SeqLength1000bp.yml \
--learningRate 0.005 \
--bedfile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/Seedling_18Celltypes.500.RestrictACR6CT_Sorted.bed.gz \
--featurefile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/Seedling_18Celltypes.500.RestrictACR6CT_distinctfeatures.txt \
--OutwmlfileName /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/Test.wml \
--NewOutputDir /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/Test
'''

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
    nFeatures = len(open(args.feature,"r").readlines())
    for sLine in infile:
        if "n_targets" in sLine:
            outfile.write(sLine.replace("n_targets:","n_targets: "+str(nFeatures)))
        elif "input_path_forfeatures" in sLine:
            outfile.write(sLine.replace("input_path_forfeatures:","input_path: "+str(args.feature)))
        elif "target_path:" in sLine:
            outfile.write(sLine.replace("target_path:","target_path: "+str(args.bed)))
        elif "output_dir:" in sLine:
            outfile.write(sLine.replace("output_dir:","output_dir: "+str(args.OutDir)))
        else:
            outfile.write(sLine)
    infile.close()
    outfile.close()

## Main Run!
args = get_parser().parse_args()
Modify_ymlFiles()

configs = load_path(args.Outwml)
# deeperdeepsea was lr=0.01
parse_configs_and_run(configs, lr=args.lr)
