#!/usr/bin/env python3
"""Exact event costs, replay, boundaries, and unchanged reconstruction outputs."""
import itertools
import math
from io import StringIO
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT))
import pandas as pd
import medicc
from Bio import Phylo
from medicc.full_fst import reconstruct
from medicc.full_fst.cli import main as cli
from medicc.full_fst.events import EventCaller, replay, call_tree_events
from oracle import get_parts, Primitive, relation, minplus


def exhaustive():
    parts, symbols = get_parts()
    count = 0
    with EventCaller() as caller:
        for separator in ['', 'X']:
            words = [separator.join(pair) for pair in itertools.product('012345678', repeat=2)]
            p = {k: relation(Primitive(parts[k]), words, symbols) for k in ['L_LOH', 'G', 'W1']}
            a = minplus(p['L_LOH'], p['G'])
            t = a
            for _ in range(2):
                t = minplus(minplus(t, p['W1']), a)
            for i, x in enumerate(words):
                for j, y in enumerate(words):
                    if math.isfinite(t[i, j]):
                        h = caller.branch(x, y, expected_score=t[i, j])
                        assert h.cost == len(h.events) == t[i, j], (x, y, h)
                        assert replay(x, h.events) == y, (x, y, h)
                        assert h.wgds == sum(e.type == 'wgd' for e in h.events) <= 2
                    else:
                        try:
                            caller.branch(x, y)
                        except ValueError as exc:
                            assert 'No feasible' in str(exc)
                        else:
                            raise AssertionError((x, y, 'accepted infeasible history'))
                    count += 1
                caller.scorer.clear()
            print('Checked', count, 'literal-oracle event histories', flush=True)
        for x, y, cost, wgds in [
            ('21', '10', 1, 0), ('303', '000', 3, 0),
            ('102', '203', 1, 0), ('1X1', '4X4', 2, 2),
            ('555X1', '888X2', 1, 1), ('12', '12', 0, 0),
            ('111X111', '222X222', 1, 1),
        ]:
            h = caller.branch(x, y, expected_score=cost)
            assert h.wgds == wgds and replay(x, h.events) == y
        for high in '5678':
            assert caller.branch(high+'X1X1', '8X2X2').wgds == 1
        for x, y in [('1X1', '111'), ('9', '8'), ('', ''), ('1', '11')]:
            try:
                caller.branch(x, y)
            except ValueError:
                pass
            else:
                raise AssertionError('Accepted malformed genomes')
        try:
            caller.branch('1', '2', expected_score=7)
        except ValueError:
            pass
        else:
            raise AssertionError('Accepted mismatching branch score')


def workflow():
    with tempfile.TemporaryDirectory(prefix='medicc-events-check-') as temp:
        tmp = Path(temp)
        path, treefile = tmp/'input.tsv', tmp/'tree.new'
        rows = [(name, chrom, 10, 19, a, b)
                for name, a, b in [('alpha', 1, 1), ('beta', 2, 2), ('gamma', 4, 4)]
                for chrom in ['chr1', 'chr2']]
        pd.DataFrame(rows, columns=['sample_id','chrom','start','end','cn_a','cn_b']).to_csv(path, sep='\t', index=False)
        treefile.write_text('(alpha,(beta,gamma));\n')
        frame = medicc.io.read_and_parse_input_data(path, normal_name='normal')
        for nni in [False, True]:
            args = dict(normal_name='normal', input_tree=treefile, n_cores=2, nni=nni)
            original = reconstruct(frame, **args)
            result = reconstruct(frame, events=True, **args)
            assert original.profiles == result.profiles
            pd.testing.assert_frame_equal(original.copy_numbers, result.copy_numbers)
            assert original.tree.format('newick') == result.tree.format('newick')
            assert original.details['tree_length'] == result.details['tree_length'] == len(result.events)
            if nni:
                for key in ['trace', 'best_score', 'step_records', 'best_ancestors']:
                    assert original.nni_result[key] == result.nni_result[key], key
            else:
                assert original.nni_result is result.nni_result is None
            assert list(result.events.reset_index().columns[:7]) == [
                'sample_id','chrom','start','end','allele','type','cn_child']
        try:
            reconstruct(frame, normal_name='normal', events=True, topology_only=True)
        except ValueError:
            pass
        else:
            raise AssertionError('Accepted events without ancestors')
        for name, flag in [('plain', []), ('events', ['--events'])]:
            cli([str(path), str(tmp/name), '--tree', str(treefile), '--normal-name', 'normal',
                 '--prefix', 'check', '--no-plot-tree', *flag])
        old, new = tmp/'plain', tmp/'events'
        assert {p.name for p in new.iterdir()} - {p.name for p in old.iterdir()} == {'check_copynumber_events_df.tsv'}
        for p in old.iterdir():
            if p.suffix not in ('.json', '.gz'):
                assert p.read_bytes() == (new/p.name).read_bytes(), p.name
        table = pd.read_csv(new/'check_copynumber_events_df.tsv', sep='\t')
        assert len(table) == result.details['tree_length']
        # An identity tree still has a usable empty table with the same schema.
        tree = Phylo.read(StringIO('(child:0)parent;'), 'newick')
        empty = call_tree_events(tree, {'parent':'1X1', 'child':'1X1'}, [('chr1',10,19)], ['cn_a','cn_b'])
        assert empty.empty and 'cn_child' in empty.columns
        # Independent alleles share one genome-wide WGD, never two local events.
        tree.root.clades[0].branch_length = 1
        table = call_tree_events(tree, {'parent':'5X1', 'child':'8X2'}, [('chr1',10,19)], ['cn_a','cn_b'])
        assert len(table) == 1 and table.iloc[0]['type'] == 'wgd'
        assert table.index[0] == ('child', 'chr0', 10, 19)
        assert table.iloc[0]['allele'] == 'both'
    print('PASS: API/NNI/CLI unchanged outputs, legacy event filename/schema, coordinates and empty histories', flush=True)


if __name__ == '__main__':
    exhaustive()
    workflow()
