""" Command-line interface

Sub-commands:
    1. model - train the BMR model.
    2. infer - test for driver elements.
"""


import logging
import argparse
import sys
from driverpower import __version__
from driverpower.BMR import run_bmr
from driverpower.inference import make_inference

logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(asctime)s | %(levelname)s: %(message)s',
                    datefmt='%m/%d/%Y %H:%M:%S')
# create logger
logger = logging.getLogger('DP')


def get_args():
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.MetavarTypeHelpFormatter):
        pass
    parser = argparse.ArgumentParser(prog='driverpower',
                                     description='DriverPower v{}: Combined burden and functional impact '
                                                 'tests for coding and non-coding cancer driver discovery'.format(__version__))
    # global argument
    parser.add_argument('-v', '--version', dest='version', action="store_true",
                        help='Print the version of DriverPower')
    subparsers = parser.add_subparsers(title='The DriverPower sub-commands include',
                                       dest='subcommand')
    #
    # Background model
    #
    parser_bmr = subparsers.add_parser('model',
                                       help='Build the background mutation model',
                                       formatter_class=CustomFormatter)
    # Input data
    dat_bmr = parser_bmr.add_argument_group(title="Input data")
    dat_bmr.add_argument('--feature', dest='X_path', required=True, type=str,
                         help='Path to the training feature table')
    dat_bmr.add_argument('--response', dest='y_path', required=True, type=str,
                         help='Path to the training response table')
    dat_bmr.add_argument('--featImp', dest='fi_path', required=False, type=str,
                         help='Path to the feature importance table', default=None)
    # Parameters
    par_bmr = parser_bmr.add_argument_group(title="Parameters")
    par_bmr.add_argument('--method', dest='model_name', required=True, type=str,
                         help='Algorithms to use', choices=['GLM', 'GBM'])
    par_bmr.add_argument('--featImpCut', dest='fi_cut', required=False, type=float,
                         help='Cutoff of feature importance score', default=0.5)
    par_bmr.add_argument('--name', dest='project_name', required=False, type=str,
                         help='Identifier for output files', default='DriverPower')
    par_bmr.add_argument('--outDir', dest='out_dir', type=str,
                         help='Directory of output files', default='./output/')
    #
    # Inference
    #
    parser_infer = subparsers.add_parser('infer',
                                         help='Infer driver elements',
                                         formatter_class=CustomFormatter)
    # Input data
    dat_infer = parser_infer.add_argument_group(title="Input data")
    # Required data
    dat_infer.add_argument('--feature', dest='X_path', required=True, type=str,
                           help='Path to the test feature table')
    dat_infer.add_argument('--response', dest='y_path', required=True, type=str,
                           help='Path to the test response table')
    dat_infer.add_argument('--model', dest='model_path', required=True, type=str,
                           help='Path to the trained model')
    dat_infer.add_argument('--modelInfo', dest='model_info_path', required=True, type=str,
                           help='Path to the model information')
    # Optional data
    dat_infer.add_argument('--featImp', dest='fi_path', required=False, type=str,
                           help='Path to the feature importance table', default=None)
    dat_infer.add_argument('--funcScore', dest='fs_path', required=False, type=str,
                           help='Path to the functional score table', default=None)
    # Parameters
    par_infer = parser_infer.add_argument_group(title="Parameters")
    par_infer.add_argument('--method', dest='test_method', required=False, type=str,
                           help='Test method to use', choices=['auto', 'binomial', 'negative_binomial'], default='auto')
    par_infer.add_argument('--featImpCut', dest='fi_cut', required=False, type=float,
                           help='Cutoff of feature importance score', default=0.5)
    par_infer.add_argument('--scale', dest='scale', required=False, type=float,
                           help='Scaling factor for theta in negative binomial distribution', default=1)
    par_infer.add_argument('--funcScoreCut', dest='fs_cut', required=False, type=str,
                           help='Score name:cutoff pairs for all scores e.g., "CADD:0.01;DANN:0.03;EIGEN:0.3"',
                           default=None)
    par_infer.add_argument('--geoMean', dest='use_gmean', required=False, type=bool,
                           help='Use geometric mean in test', default=True)
    par_infer.add_argument('--name', dest='project_name', required=False, type=str,
                           help='Identifier for output files', default='DriverPower')
    par_infer.add_argument('--outDir', dest='out_dir', type=str,
                           help='Directory of output files', default='./output/')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    logger.info('DriverPower {}'.format(__version__))
    if args.subcommand == 'model':
        run_bmr(model_name=args.model_name,
                X_path=args.X_path,
                y_path=args.y_path,
                fi_cut=args.fi_cut,
                fi_path=args.fi_path,
                project_name=args.fi_name,
                out_dir=args.out_dir)
    elif args.subcommand == 'infer':
        make_inference(args)


if __name__ == '__main__':
    main()