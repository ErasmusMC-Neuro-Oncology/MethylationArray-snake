#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PythonScript.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# <Objective>
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: 
# Usage: 
"""
python3 scripts/PythonScript.py \
        -i {input} \
        -o {output}
"""
#
# TODO:
# 1) 
#
# History:
#  05-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
import argparse

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 PythonScript.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Create tiled WSI for Prov_GigaPath  ')
    parser.add_argument('-i', help='path to input file',
                        dest='input',
                        type=str)
    parser.add_argument('-o', help='path to output file',
                        dest='output',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

"""
args.input = '~/mnt/neuro-genomic-1-ro/'
args.output = 'output/'
"""

#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
