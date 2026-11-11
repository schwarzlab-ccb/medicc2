import argparse
import logging
import os

import fstlib
import medicc

logger = logging.getLogger('tools')

objects_dir = os.path.join(os.path.dirname(os.path.realpath(medicc.__file__)), "objects")

parser = argparse.ArgumentParser()
parser.add_argument("input_file", 
                    help = "a path to the input file")
parser.add_argument("output_file", 
                    help ="a path to the output file")
parser.add_argument("--input-type", "-t", 
                    type = str, 
                    dest = "input_type", 
                    default = "t", 
                    choices = ["f", "t"], 
                    required = False, 
                    help = "Choose the type of input: f for FASTA, t for TSV (default: TSV)")
parser.add_argument("--input-allele-columns", "-a",
                    type=str,
                    dest='input_allele_columns',
                    default='cn_a, cn_b',
                    required=False,
                    help = """Name of the CN columns (comma separated) if using TSV input format (default: 'cn_a, cn_b').
                    This also adjusts the number of alleles considered (min. 1, max. 2).""")
parser.add_argument("--normal-name", "-n",
                    default = "diploid", 
                    type = str,
                    dest = "normal_name", 
                    help = """ID of the sample to be treated as the normal sample. 
                    Trees are rooted at this sample for ancestral reconstruction (default: \"diploid\").
                    If the sample ID is not found, an artificial normal sample of the same name is created 
                    with CN states = 1 for each allele.""",
                    required = False)
parser.add_argument("--prefix", '-p', type=str, dest='prefix', default=None, 
                    help='Output prefix to be used (default: input filename).', required=False)
parser.add_argument("--wgd", '-w', action='store_true', default=False, 
                    help='Enable whole-genome doubling events (default: false).', required=False)
parser.add_argument("--fst", type=str, dest='fst', default=None,
                    help='Expert option: path to an alternative FST.')
parser.add_argument("--fst-chr-separator", type=str, dest='fst_chr_separator', default='X',
                    help = 'Expert option: character used to separate chromosomes in the FST (default: \"X\").')
parser.add_argument("--maxcn", type=int, dest='maxcn', default=8,
                    help='Expert option: maximum CN supported by the supplied FST.')

args = parser.parse_args()

input_file = args.input_file
output_file = args.output_file
normal_name = args.normal_name 
allele_columns = [x.strip() for x in args.input_allele_columns.split(',')]

## Determine prefix for output files 
if args.prefix is None:
    output_prefix = os.path.basename(os.path.splitext(input_file)[0])
else:
    output_prefix = args.prefix

## Read in symbol table and transducers (phasing uses an fstlib FST)
logger.info("Reading FST.")
if args.fst is None:
    fst_path = os.path.join(objects_dir, 'wgd_asymm.fst' if args.wgd else 'no_wgd_asymm.fst')
else:
    fst_path = args.fst
fst = medicc.io.read_fst(user_fst=fst_path)

## Load data
logger.info("Reading and parsing input data.")
input_df = medicc.io.read_and_parse_input_data(
    args.input_file,
    normal_name,
    args.input_type.strip(),
    allele_columns=allele_columns)

## Run phasing method
logger.info("Phasing.")
output_df = medicc.phase(input_df, fst, normal_name, args.fst_chr_separator)

## Write ouput table
logger.info("Writing output table.")
output_df.to_csv(output_file, sep='\t')
