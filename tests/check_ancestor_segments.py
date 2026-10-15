"""Compare full weighted languages and assert segmented/fallback path coverage."""
import argparse,itertools,json,sys
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('tool');ap.add_argument('out');ap.add_argument('--expected');ap.add_argument('--forced',action='store_true');args=ap.parse_args()
root=Path(args.tool);sys.path[:0]=[str(root),str(root/'tests')]
from check_ancestor_bounds import weighted
from medicc.full_fst.model import Model,ExecutionBudget
cases={
 'factorized_two':({a+'X'+b:i+2*j for i,a in enumerate('24') for j,b in enumerate('48')},2),
 'correlated_two':({'2X4':0,'4X8':1},0),
 'factorized_three':({a+'X'+b+'X'+c:i+2*j+k for i,a in enumerate('24') for j,b in enumerate('48') for k,c in enumerate('28')},3),
 'only_first_cut':({a+'X'+b+'X'+c:i+j for i,a in enumerate('24') for j,(b,c) in enumerate([('2','4'),('4','8')])},2),
 'only_second_cut':({a+'X'+b+'X'+c:i+j for i,(a,b) in enumerate([('2','4'),('4','8')]) for j,c in enumerate('28')},2),
 'mixed_separator_layout':({'2X42':0,'24X2':1},0),
 'consecutive_separators':({'2XX4':0,'4XX4':2},3),
 'trailing_separator':({'2X':0,'4X':2},0),
 'leading_separator':({'X2X4':0,'X4X4':2},3),
}
model=Model();records={}
try:
 for name,(language,cuts) in cases.items():
  left,lv=weighted(model,language,3);right,rv=weighted(model,{w:c+1 for w,c in language.items()},7)
  outputs=[]
  for workers in [1,2]:
   model.execution_budget=ExecutionBudget(workers);out=model.join(left,right);outputs.append(out.language());out.close()
   branches=model.stats()['join_trace']['branches']
   if args.forced:
    counts=[b['segments'] for b in branches if b.get('segment_restarts')]
    if cuts: assert counts and all(n==cuts for n in counts),(name,cuts,branches)
    else: assert not counts,(name,branches)
  assert outputs[0]==outputs[1],name
  records[name]=outputs[0];left.close();right.close();model.execution_budget=None
finally:model.close()
if args.expected: assert records==json.loads(Path(args.expected).read_text())
Path(args.out).write_text(json.dumps(records,sort_keys=True)+'\n')
print(json.dumps(dict(status='ok',cases=len(cases),weighted_all_tie_joins=2*len(cases),trace_coverage=args.forced)))
