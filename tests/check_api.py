"""General MEDICC inputs, normal creation, supplied topology, and validation errors."""
import csv
from io import StringIO
from pathlib import Path
import sys
import tempfile
ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT))
import medicc
from medicc.full_fst import reconstruct
from Bio import Phylo

def main():
    with tempfile.TemporaryDirectory(prefix='medicc-api-check-') as tmp:
        path = Path(tmp)/'input.tsv'
        with path.open('w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['sample_id','chrom','start','end','cn_a','cn_b'])
            for name, cn in [('alpha',1),('beta',2),('gamma',3)]:
                writer.writerow([name,'chr1',0,99,cn,1])
                writer.writerow([name,'chr2',0,99,1,cn])
        frame = medicc.io.read_and_parse_input_data(path, normal_name='normal')
        assert 'normal' in frame.index.get_level_values('sample_id')
        tree = Phylo.read(StringIO('(alpha,(beta,gamma));'), 'newick')
        for nni in [False, True]:
            result = reconstruct(frame, normal_name='normal', input_tree=tree, n_cores=2, nni=nni)
            assert result.pairwise_distances is None
            assert result.tree.root.name == 'normal'
            assert len(result.profiles) == 6
            for sample in ['alpha','beta','gamma','normal']:
                assert result.copy_numbers.loc[sample].astype(str).equals(frame.loc[sample].astype(str))
        wrong = frame.copy(); wrong.iloc[0,0] = '-1'
        try: reconstruct(wrong, normal_name='normal', input_tree=tree)
        except ValueError: pass
        else: raise AssertionError('Negative CN accepted')
        try: reconstruct(frame.iloc[1:], normal_name='normal', input_tree=tree)
        except ValueError: pass
        else: raise AssertionError('Inconsistent grid accepted')
    print('PASS: arbitrary sample labels, automatic normal, supplied tree, NNI, observed profiles, invalid inputs')
if __name__ == '__main__': main()
