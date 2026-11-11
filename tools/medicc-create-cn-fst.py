import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import argparse
import logging

import medicc

# prepare logger and parse command line arguments
logger = logging.getLogger('tools')

parser = argparse.ArgumentParser()
parser.add_argument("output_folder", type=str, help="output folder")
parser.add_argument("--max-cn", "-m", type=int, required=False, default=8, help="Maximum copy-number per allele (default:8).")
parser.add_argument("--sep", "-s", type=str, required=False, default='X', help="Chromosome separator (default:\"X\")")
parser.add_argument("--wgd", action='store_true', required=False, default=False)
parser.add_argument("--wgd_x2", action='store_true', required=False, default=False)
parser.add_argument("--total_cn", action='store_true', required=False, default=False)
parser.add_argument("--max-num-wgds", type=int, required=False, default=3, help="Maximum number of WGD events (Default: 3)")
parser.add_argument("--prefix", "-p", action='store', required=False, default='fst')
parser.add_argument("--write-symbol-table", action='store_true', required=False, default=False)
parser.add_argument("--legacy-loh", action='store_true', required=False, default=False)
args = parser.parse_args()

separator = args.sep

logger.info('Creating symbol table.')
symbol_table = medicc.create_symbol_table(args.max_cn, separator)
logger.info('Symbol table: %s', str(list(symbol_table)))

logger.info('Creating FSTs.')
fst = medicc.create_copynumber_fst(symbol_table=symbol_table, sep=separator, 
                                   enable_wgd=args.wgd, max_num_wgds=args.max_num_wgds,
                                   wgd_x2=args.wgd_x2, total_cn=args.total_cn,
                                   exact_nowgd=not args.legacy_loh)
logger.info('FST: %d states.', fst.num_states())

logger.info('Writing.')
if args.write_symbol_table:
    symbol_table.write_text(os.path.join(args.output_folder, "%s_symbols.txt" % args.prefix))
fst.write(os.path.join(args.output_folder, "%s_asymm.fst" % args.prefix))
logger.info('Done.')

