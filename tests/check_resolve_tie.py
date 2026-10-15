"""Preserve the historical representative when compact WGD histories tie."""
from pathlib import Path
import ctypes as C
import json,sys
ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from medicc.full_fst.model import Model
from medicc.full_fst.message import Message,INF


def main():
    # Both 2264 and 3186 cost 8. Minimized AD factors have a different state
    # ordering, so downward resolution deliberately retains original A D A.
    language={'2264':0,'5572':1,'7662':2,'3186':1,'5781':2}
    nodes=[dict(final=INF,arcs={})]
    for word,cost in language.items():
        state=0
        for digit in word:
            symbol=int(digit)
            if symbol not in nodes[state]['arcs']:
                nodes[state]['arcs'][symbol]=len(nodes)
                nodes.append(dict(final=INF,arcs={}))
            state=nodes[state]['arcs'][symbol]
        nodes[state]['final']=cost
    rows=[]
    for node in nodes:
        row=[node['final']]+[-1,0]*10
        for symbol,target in node['arcs'].items():row[1+2*symbol]=target
        rows.extend(row)
    data=(C.c_int*len(rows))(*rows)
    model=Model()
    message=Message(model,model.lib.am_import(4,len(nodes),0,0,data))
    try:
        assert message.language()==language
        result=model.resolve('5448',message)
        assert result==('2264',8),result
        cost=C.c_int();stats=(C.c_double*3)()
        word=model.lib.am_direct_resolve(model.handle,b'5448',message.handle,INF,C.byref(cost),stats)
        assert (word.decode(),cost.value)==('2264',8),(word,cost.value)
    finally:
        message.close();model.close()
    print(json.dumps(dict(status='ok',parent='5448',representative='2264',cost=8)))


if __name__=='__main__':main()
