#!/usr/bin/env python3
"""BED coverage, WGD multiplicity, and CLI annotation compatibility."""
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT))
import pandas as pd
from medicc.event_reconstruction import overlap_events
from medicc.full_fst.cli import main as cli


def event_table(rows):
    columns = ['sample_id','chrom','start','end','allele','type','cn_child']
    return pd.DataFrame(rows, columns=columns).set_index(columns[:4])


def checks(tmp):
    arms, regions = tmp/'arms.bed', tmp/'regions.bed'
    arms.write_text('chr1\t0\t50\tchr1p\nchr1\t50\t100\tchr1q\n')
    regions.write_text('chr1\t0\t100\tGENE\nchr2\t0\t100\tOTHER\n')
    events = event_table([
        ('a','chr1',0,90,'cn_a','gain',2),  # Exactly 90% does NOT qualify.
        ('b','chr1',0,91,'cn_a','loss',0),
        ('c','chr1',0,100,'cn_b','gain',2),
        ('w1','chr0',0,100,'both','wgd',0),
        ('w1','chr0',0,100,'both','wgd',0),
        ('w2','chr0',0,100,'both','wgd',0),
    ])
    result = overlap_events(events_df=events.copy(), chromosome_bed=str(arms), regions_bed=str(regions))
    genes = result.loc[result['name']=='GENE']
    assert set(genes.index) == {'b','c'}
    assert genes.loc['b','FractionOverlaps'] == 0.91
    assert 'chr1p' in set(result.loc[['a'],'name'])
    assert set(result.loc[['c'],'name']) == {'chr1','GENE'}  # Whole chromosome replaces arms.
    wgds = result.loc[result['event']=='WGD']
    assert list(wgds.index) == ['w1','w1','w2']
    # Coverage is the union of same-type events, not a sum of overlap lengths.
    union = event_table([('u','chr1',0,60,'cn_a','gain',2),
                         ('u','chr1',40,100,'cn_a','gain',2)])
    result = overlap_events(events_df=union, chromosome_bed='none', regions_bed=str(regions))
    assert len(result) == 1 and result.iloc[0]['FractionOverlaps'] == 1


def edge_checks(tmp):
    regions = tmp/'regions.bed'
    rows = [('u','chr1',0,60,'cn_a','gain',2), ('u','chr1',40,100,'cn_b','gain',2)]
    result = overlap_events(events_df=event_table(rows), chromosome_bed='none', regions_bed=str(regions))
    assert result.empty
    for rows, expected in [([],0), ([('w','chr0',0,100,'both','wgd',0)]*2,2),
                           ([('x','chr9',0,100,'cn_a','gain',2)],0)]:
        result = overlap_events(events_df=event_table(rows), chromosome_bed='none', regions_bed=str(regions))
        assert len(result) == expected
        assert list(result.columns) == ['Chromosome','Event-Start','Event-End','name','FractionOverlaps','event','final_name']


def workflow(tmp):
    input_file, tree = tmp/'input.tsv', tmp/'tree.new'
    rows = [(sample,'chr1',start,start+50,cn,1)
            for sample,cn in [('alpha',1),('beta',2),('gamma',3)] for start in [0,50]]
    pd.DataFrame(rows, columns=['sample_id','chrom','start','end','cn_a','cn_b']).to_csv(input_file, sep='\t', index=False)
    tree.write_text('(alpha,(beta,gamma));\n')
    common = ['--tree',str(tree),'--normal-name','normal','--prefix','check','--no-plot-tree']
    variants = {
        'plain': [],
        'custom': ['--chromosomes-bed',str(tmp/'arms.bed'),'--regions-bed',str(tmp/'regions.bed')],
        'default': ['--chromosomes-bed','default','--regions-bed','default'],
        'hg19': ['--chromosomes-bed',str(ROOT/'medicc/objects/hg19_chromosome_arms.bed')],
    }
    for name, options in variants.items():
        cli([str(input_file),str(tmp/name),*common,'--events',*options])
    plain, annotated = tmp/'plain', tmp/'custom'
    assert {p.name for p in annotated.iterdir()} - {p.name for p in plain.iterdir()} == {'check_events_overlap.tsv'}
    for p in plain.iterdir():
        if p.suffix not in ('.json','.gz'):
            assert p.read_bytes() == (annotated/p.name).read_bytes(), p.name
    output = pd.read_csv(annotated/'check_events_overlap.tsv', sep='\t')
    assert 'GENE +' in set(output['final_name'])
    assert (tmp/'default/check_events_overlap.tsv').is_file()
    assert (tmp/'hg19/check_events_overlap.tsv').is_file()
    # Invalid combinations and missing BED files fail before creating outputs.
    for name, options in [
        ('noevents',['--regions-bed','default']),
        ('missing',['--events','--regions-bed',str(tmp/'missing.bed')]),
        ('topology',['--events','--topology-only','--regions-bed','default']),
    ]:
        try:
            cli([str(input_file),str(tmp/name),*common,*options])
        except SystemExit as exc:
            assert exc.code == 2
        else:
            raise AssertionError('Accepted invalid BED arguments')
        assert not (tmp/name).exists()


if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='medicc-bed-check-') as temp:
        tmp = Path(temp)
        checks(tmp)
        edge_checks(tmp)
        workflow(tmp)
    print('PASS: BED coverage/alleles, repeated WGDs, empty overlaps, bundled files and unchanged CLI outputs')
