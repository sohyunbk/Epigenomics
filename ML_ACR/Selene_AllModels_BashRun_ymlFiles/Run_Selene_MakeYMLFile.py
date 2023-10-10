from selene_sdk.utils import load_path
from selene_sdk.utils import parse_configs_and_run
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('-wmlFile', "--wmlFile", help="It should be Model_wmlFile.", required=True, dest='wml')
    parser.add_argument('-learningRate', "--learningRate", help="learningRate.", required=True, dest='lr')
    parser.add_argument('-learningRate', "--learningRate", help="learningRate.", required=True, dest='lr')
    parser.add_argument('-learningRate', "--learningRate", help="learningRate.", required=True, dest='lr')

    args = vars(parser.parse_args())
    return parser


def Modify_ymlFiles():


## Main Run!
args = get_parser().parse_args()
configs = load_path(args.wml)
# deeperdeepsea was lr=0.01
parse_configs_and_run(configs, lr=args.lr)
