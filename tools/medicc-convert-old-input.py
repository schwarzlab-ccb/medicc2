#### Transform fasta format of copy number profiles to the one for new MEDICC
import sys
import glob
import pandas as pd
import os
import re
import random
import argparse 
import csv 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import medicc

parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str, help="The name of the input descriptor file (desc.txt) of the old MEDICC input.")
parser.add_argument("output_file", type=str, help="The name of the output TSV file where to store the results.")
args = parser.parse_args()

infile = args.input_file
outfile = args.output_file

result = medicc.io.read_fasta_as_dataframe(infile)

result.to_csv(outfile, sep='\t', index=False)
