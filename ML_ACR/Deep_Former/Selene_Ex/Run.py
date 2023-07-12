from selene_sdk.utils import load_path
from selene_sdk.utils import parse_configs_and_run
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('-wmlFile', "--wmlFile", help="wmlFile.", required=True, dest='wml')
    args = vars(parser.parse_args())
    return parser

args = get_parser().parse_args()


configs = load_path("/home/sb14489/ACR_ML_caQTLs/Deep_Former/Selene_Ex/"+args.wml)

parse_configs_and_run(configs, lr=0.01)
