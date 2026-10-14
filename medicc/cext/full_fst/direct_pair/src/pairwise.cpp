// Exact fixed-string distances using the validated full L_LOH primitive cascade.
// No eager pullbacks, weighted subset construction, or incumbent-engine dependency.
#include "direct_trace.hpp"

namespace pair_direct {
static int changes(const vector<int>& a, const vector<int>& b) {
    int previous_loss=0,previous_gain=0,value=0;
    for(size_t i=0;i<a.size();++i) {
        if(a[i]<0) {previous_loss=previous_gain=0;continue;}
        int loss=max(0,a[i]-b[i]),gain=max(0,b[i]-a[i]);
        value+=max(0,loss-previous_loss)+max(0,gain-previous_gain);
        previous_loss=loss;previous_gain=gain;
    }
    return value;
}
// A concrete history: choose a binwise ancestor, apply k exact global doublings,
// then interval unit losses. All intermediate copy numbers remain in 0..8.
static int witness(const string& x,const string& y,bool directed) {
    vector<int>a(x.size()),b(y.size());
    for(size_t i=0;i<x.size();++i) {
        a[i]=x[i]=='X'?-1:x[i]-'0';b[i]=y[i]=='X'?-1:y[i]-'0';
        if(directed && a[i]==0 && b[i]>0)return INF;
    }
    if(directed)return changes(a,b); // zero-WGD losses then gains, no resurrection
    int best=INF;
    for(int ka=0;ka<=2;++ka)for(int kb=0;kb<=2;++kb) {
        int fa=1<<ka,fb=1<<kb;vector<int>pa(a.size()),pb(b.size());bool valid=true;
        for(size_t i=0;i<a.size();++i) {
            if(a[i]<0){pa[i]=pb[i]=-1;continue;}
            int z=max((a[i]+fa-1)/fa,(b[i]+fb-1)/fb);
            if(z*max(fa,fb)>8){valid=false;break;}
            pa[i]=z*fa;pb[i]=z*fb;
        }
        if(valid)best=min(best,ka+kb+changes(pa,a)+changes(pb,b));
    }
    return best;
}
}
extern "C" int am_pair_score_only(){return 1;}
extern "C" int am_pair_search_uncapped(){return 1;}
extern "C" double am_pair_score_whole(void *handle,const char *xx,const char *yy,int directed,double *stat) {
    try {
        string x=xx,y=yy;
        if(x.size()!=y.size())return INFINITY;
        for(size_t i=0;i<x.size();++i) {
            for(char c:{x[i],y[i]})if(c!='X'&&(c<'0'||c>'8'))throw runtime_error("invalid symbol");
            if((x[i]=='X')!=(y[i]=='X'))return INFINITY;
        }
        int best=pair_direct::witness(x,y,directed);
        if(best==0||best>=INF)return best>=INF?INFINITY:double(best);
        auto *m=(AM*)handle;
        DFA a=linear(x),b=linear(y);
        if(directed) {
            auto branches=direct_nfa::branches(m,b,0,false);
            vector<pair<int,int>> order;
            for(int k=0;k<int(branches.size());++k) {
                auto&side=*branches[k];order.push_back({side.wg+side.h(side.start()),k});
            }
            sort(order.begin(),order.end());
            for(auto [lower,k]:order) {
                if(lower>=best)continue;
                stat[5]++;
                auto r=direct_trace::search(*branches[k],x,best-1,stat);
                if(r.cost<best)best=r.cost;
            }
        } else {
            auto left=direct_nfa::branches(m,a,0),right=direct_nfa::branches(m,b,0);
            auto lm=lazy_cascade::layout(a),rm=lazy_cascade::layout(b);
            vector<int>mask(a.len);
            for(int i=0;i<a.len;++i)mask[i]=(x[i]=='X')?512:(510|((x[i]=='0'&&y[i]=='0')?1:0));
            vector<array<int,3>>order;
            vector<vector<int>>futures(9);
            for(int ka=0;ka<3;++ka)for(int kb=0;kb<3;++kb) {
                auto&l=*left[ka];auto&r=*right[kb];
                auto&future=futures[3*ka+kb];future=chromosome_bound::future(lm,rm,3*ka+kb);
                int lower=ka+kb+max(l.h(l.start())+r.h(r.start()),l.chrom_h(l.start())+r.chrom_h(r.start())+future[0]);
                order.push_back({lower,ka,kb});
            }
            sort(order.begin(),order.end());
            for(auto [lower,ka,kb]:order) {
                if(lower>=best)continue;
                stat[5]++;
                direct_nfa::Stats stats{stat};int cost=INF;
                direct_nfa::search(*left[ka],*right[kb],mask,best,0,stats,futures[3*ka+kb],nullptr,&cost);
                best=min(best,cost);
            }
        }
        return double(best);
    }catch(exception&e){error=e.what();return -1;}
}

namespace factorized {
// X resets every primitive to its start with zero weight. For fixed ka,kb,
// all chromosomes are independent. Sum their event costs, add ka+kb ONCE,
// and only then minimize the nine global WGD choices.
struct Chunk {
    string x,y;
    array<int,9> lower{};
    array<int,9>* exact;
    bool flipped;
    DFA a,b;
    vector<unique_ptr<direct_nfa::Side>> left,right;
    vector<int> mask,future;
    int cache_branch(int branch) const {return flipped?3*(branch%3)+branch/3:branch;}
    Chunk(AM*m,string xx,string yy):x(move(xx)),y(move(yy)) {
        flipped=y<x;
        string key=flipped?y:x;key.push_back('\0');key+=flipped?x:y;
        array<int,9> unset;unset.fill(-1);
        exact=&m->scalar_cache.emplace(move(key),unset).first->second;
        for(int k=0;k<9;++k) {
            int previous=-1,total=0;
            for(size_t i=0;i<x.size();++i) {
                int u=x[i]-'0',v=y[i]-'0',state=9*u+v;
                if(i==0)total=chromosome_bound::END[k][state];
                else total+=chromosome_bound::EDGE[k][(previous/9*9+u)*81+previous%9*9+v];
                previous=state;
            }
            total+=chromosome_bound::END[k][previous];
            lower[k]=(total+1)/2;
            int cached=(*exact)[cache_branch(k)];
            if(cached>=0)lower[k]=cached;
        }
    }
    void prepare(AM*m) {
        if(!left.empty())return;
        a=linear(x);b=linear(y);
        left=direct_nfa::branches(m,a,0);right=direct_nfa::branches(m,b,0);
        mask.resize(x.size());future.assign(x.size()+1,0);
        for(size_t i=0;i<x.size();++i)mask[i]=510|((x[i]=='0'&&y[i]=='0')?1:0);
    }
    // Return the exact local cost if it is strictly below room; otherwise INF.
    // A cutoff failure is never cached as an exact score.
    int score(AM*m,int branch,int room,double*stat) {
        int&cached=(*exact)[cache_branch(branch)];
        if(cached>=0){stat[3]++;return cached;}
        stat[4]++;
        prepare(m);
        int ka=branch/3,kb=branch%3,cost=INF;
        direct_nfa::Stats stats{stat};
        direct_nfa::search(*left[ka],*right[kb],mask,room+ka+kb,0,stats,future,nullptr,&cost);
        if(cost<INF){cached=cost-ka-kb;lower[branch]=cached;return cached;}
        return INF;
    }
};
static int score(AM*m,const string&x,const string&y,double*stat) {
    int best=pair_direct::witness(x,y,false);
    if(best==0)return best;
    vector<unique_ptr<Chunk>> chunks;
    for(size_t i=0;i<x.size();) {
        if(x[i]=='X'){++i;continue;}
        size_t end=x.find('X',i);if(end==string::npos)end=x.size();
        chunks.push_back(make_unique<Chunk>(m,x.substr(i,end-i),y.substr(i,end-i)));
        i=end;
    }
    vector<pair<int,int>> branches;
    for(int k=0;k<9;++k) {
        int bound=k/3+k%3;
        for(auto&c:chunks)bound+=c->lower[k];
        branches.push_back({bound,k});
    }
    sort(branches.begin(),branches.end());
    for(auto [initial,k]:branches) {
        if(initial>=best)continue;
        stat[5]++;
        int bound=k/3+k%3;
        vector<Chunk*> order;
        for(auto&c:chunks){bound+=c->lower[k];order.push_back(c.get());}
        sort(order.begin(),order.end(),[k](Chunk*a,Chunk*b){return a->lower[k]>b->lower[k];});
        for(auto*c:order) {
            if(bound>=best)break;
            int before=c->lower[k],room=best-(bound-before);
            int value=c->score(m,k,room,stat);
            if(value>=INF){bound=INF;break;}
            bound+=value-before;
        }
        best=min(best,bound);
    }
    return best;
}
}
extern "C" void am_pair_clear(void*h){((AM*)h)->scalar_cache.clear();}
extern "C" size_t am_pair_cache_size(void*h){return ((AM*)h)->scalar_cache.size();}
extern "C" double am_pair_score(void*handle,const char*xx,const char*yy,int directed,double*stat) {
    if(directed)return am_pair_score_whole(handle,xx,yy,directed,stat);
    try {
        string x=xx,y=yy;
        if(x.size()!=y.size())return INFINITY;
        for(size_t i=0;i<x.size();++i) {
            for(char c:{x[i],y[i]})if(c!='X'&&(c<'0'||c>'8'))throw runtime_error("invalid symbol");
            if((x[i]=='X')!=(y[i]=='X'))return INFINITY;
        }
        return factorized::score((AM*)handle,x,y,stat);
    }catch(exception&e){error=e.what();return -1;}
}
