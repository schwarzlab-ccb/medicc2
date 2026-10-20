"""Real FastME regression for long taxon IDs, restoration and short-name parity."""
import csv
from pathlib import Path
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from Bio import Phylo
from medicc.full_fst.core import infer_tree

FASTME = os.environ.get('FASTME', '/projects/ag-schwarzr/chenxi_conda/medicc3_dev_env/bin/fastme')


def matrix(labels):
    return pd.DataFrame([[0,2,6,6],[2,0,6,6],[6,6,0,2],[6,6,2,0]],index=labels,columns=labels,dtype=float)


class FastmeLabels(unittest.TestCase):
    def test_long_names_are_bijective_and_restored(self):
        labels=['diploid','x'*64+'a','x'*64+'b','T00000001']
        with tempfile.TemporaryDirectory(prefix='fastme-labels-') as tmp:
            out=Path(tmp);tree,command=infer_tree(matrix(labels),'diploid',out,FASTME)
            self.assertEqual(tree.root.name,'diploid')
            self.assertEqual({n.name for n in tree.get_terminals()},set(labels)-{'diploid'})
            restored=Phylo.read(out/'fastme_raw.nwk','newick')
            self.assertEqual({n.name for n in restored.get_terminals()},set(labels))
            short=Phylo.read(out/'fastme_short_ids.nwk','newick')
            self.assertEqual({n.name for n in short.get_terminals()},{f'T{i:08d}' for i in range(4)})
            with (out/'fastme_taxa.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([r['sample_id'] for r in rows],labels)
            self.assertEqual(len(set(r['fastme_id'] for r in rows)),4)
            for line in (out/'distances.phy').read_text().splitlines()[1:]:self.assertLessEqual(len(line.split()[0].encode()),64)
            for i,a in enumerate(labels):
                for j,b in enumerate(labels):self.assertAlmostEqual(restored.distance(a,b),float(matrix(labels).iloc[i,j]))

    def test_63_bytes_keeps_original_interface(self):
        labels=['diploid','x'*63,'beta','gamma']
        with tempfile.TemporaryDirectory(prefix='fastme-labels-') as tmp:
            out=Path(tmp);infer_tree(matrix(labels),'diploid',out,FASTME)
            self.assertFalse((out/'fastme_taxa.tsv').exists())
            self.assertFalse((out/'fastme_short_ids.nwk').exists())
            self.assertEqual({n.name for n in Phylo.read(out/'fastme_raw.nwk','newick').get_terminals()},set(labels))

    def test_exactly64_bytes_is_also_safely_mapped(self):
        labels=['diploid','x'*64,'beta','gamma']
        with tempfile.TemporaryDirectory(prefix='fastme-labels-') as tmp:
            out=Path(tmp);infer_tree(matrix(labels),'diploid',out,FASTME)
            self.assertTrue((out/'fastme_taxa.tsv').exists())
            self.assertEqual({n.name for n in Phylo.read(out/'fastme_raw.nwk','newick').get_terminals()},set(labels))

    def test_duplicate_original_names_are_rejected(self):
        labels=['diploid','x'*65,'x'*65,'gamma']
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'unique'):infer_tree(matrix(labels),'diploid',tmp,FASTME)

    def test_column_order_mismatch_is_rejected(self):
        d=matrix(['diploid','alpha','beta','gamma']);d=d[d.columns[::-1]]
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'match in order'):infer_tree(d,'diploid',tmp,FASTME)

    def test_unexpected_or_duplicate_fastme_tips_are_rejected(self):
        def fake(command,**kwargs):
            Path(command[command.index('-o')+1]).write_text('(T00000000,T00000001,(T00000001,T00000003));')
            class Result:
                stdout=stderr=''
                def check_returncode(self):pass
            return Result()
        with tempfile.TemporaryDirectory() as tmp,patch('medicc.full_fst.core.subprocess.run',side_effect=fake):
            with self.assertRaisesRegex(ValueError,'leaf labels differ'):
                infer_tree(matrix(['diploid','a'*65,'beta','gamma']),'diploid',tmp,FASTME)


if __name__=='__main__':unittest.main(verbosity=2)
