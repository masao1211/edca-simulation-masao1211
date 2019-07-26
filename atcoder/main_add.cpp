// eratos, is_prime, def_prime_of_order, def_order_of_prime
VI poo,oop(1);int poo_max=0,oop_max=0,oop_ord=1;VI eratos(C int&n){VI p(n+1,0);fou(i,0,2)if(n>=i)p[i]=-1;fou(i,2,n)if(p[i]==0)for(int j=i*2;j<=n;j+=i)if(p[j]==0)p[j]=i;return p;}bool is_prime(int n){int m=sqrt(n-1)+1;fou(i,2,m+1)if(n%i==0)return false;return n>1;}void def_poo(C int&n){if(poo_max>=n)return;VI e=eratos(n);fou(i,poo_max+1,n+1)if(!e[i])poo.emb(i);poo_max=n;}void def_oop(C int&n){if(oop_max>=n)return;oop.resize(n+1);VI e=eratos(n);fou(i,oop_max,n+1)oop[i]=!e[i]?oop_ord++:0;oop_max=n;}

// FFT: VI convolution(VI,VI)
C int FMOD=924844033;using Fint=Mod<integral_constant<decay<decltype(FMOD)>::type,FMOD>>;C Fint FGEN=5;void fft(VI&a){int n=si(a);fou(i,0,n){int j=0,x=i,y=n-1;while(y>0){j=(j<<1)+(x&1);x>>=1;y>>=1;}if(i<j)swap(a[i],a[j]);}for(int len=1;len<n;len*=2){Fint root=power(FGEN,(FMOD-1)/(2*len));for(int i=0;i<n;i+=2*len){Fint w=1;fou(j,0,len){Fint u=a[i+j],v=a[i+j+len]*w;a[i+j]=(int)(u+v);a[i+j+len]=(int)(u-v);w*=root;}}}}VI convolution(VI a,VI b){int abn=max(si(a),si(b)),need=si(a)+si(b)-1,nn=1;while(nn<2*abn)nn<<=1;a.resize(nn);b.resize(nn);fft(a);fft(b);a=a*b%FMOD;reverse(a.begin()+1,a.end());fft(a);a=a*(int)power(Fint(nn),FMOD-2)%FMOD;a.resize(need);return a;}


// BIT: add(v,x), add(v,l,r), sum(x), sum(l,r), get(x), set(x)
TTT class BIT{public:VT t1,t2;int n;void resize(C int&_n){n=_n;t1.resize(n);t2.resize(n);}int size(){return n;}constexpr BIT(){}BIT(C int&_n){resize(_n);}void add(C T&v,C int&x){for(int i=x;i<n;i|=(i+1))t1[i]+=v;}void add(C T&v,C int&l,C int&r){for(int i=l;i<n;i|=(i+1)){t1[i]-=v*(l-1);t2[i]+=v;}for(int i=r;i<n;i|=(i+1)){t1[i]+=v*(r-1);t2[i]-=v;}}T sum(C int&x){T v=0;for(int i=x;i>=0;i=(i&(i+1))-1){v+=t1[i]+t2[i]*x;}return v;}T sum(C int&l,C int&r){return sum(r-1)-sum(l-1);}T get(C int&x){return sum(x,x+1);}void set(C T&v,C int&x){add(-get(x)+v,x);}};

// BIT2: add(v,x,y), add(v,lx,rx,ly,ry), sum(x,y), sum(lx,rx,ly,ry), get(x,y), set(x,y)
TTT class BIT2{public:VVT t1,t2x,t2y,t2xy;int n,m;void resize(VVT&t){t.resize(n);tra(e,t)e.resize(m);}void resize(C int&_n,C int&_m){n=_n;m=_m;resize(t1);resize(t2x);resize(t2y);resize(t2xy);}int size(){return n*m;}constexpr BIT2(){}BIT2(C int&_n,C int&_m){resize(_n,_m);}void add(T v,int x,int y){for(int i=x;i<n;i|=(i+1))for(int j=y;j<m;j|=(j+1)){t1[i][j]+=v;}}void _add(C T&v,C int&lrx,C int&lry,C int&pm){for(int i=lrx;i<n;i|=(i+1))for(int j=lry;j<n;j|=(j+1)){t2xy[i][j]+=pm*v;t2x[i][j]-=pm*v*(lry-1);t2y[i][j]-=pm*v*(lrx-1);t1[i][j]+=pm*v*(lrx-1)*(lry-1);}}void add(T v,C int&lx,C int&rx,C int&ly,C int&ry){_add(v,lx,ly,1);_add(v,rx,ly,-1);_add(v,lx,ry,-1);_add(v,rx,ry,1);}T sum(C int&x,C int&y){T v=0;for(int i=x;i>=0;i=(i&(i+1))-1)for(int j=y;j>=0;j=(j&(j+1))-1){v+=t1[i][j]+t2x[i][j]*x+t2y[i][j]*y+t2xy[i][j]*x*y;}return v;}T sum(C int&lx,C int&rx,C int&ly,C int&ry){return sum(rx-1,ry-1)-sum(lx-1,ry-1)-sum(rx-1,ly-1)+sum(lx-1,ly-1);}T get(C int&x,C int&y){return sum(x,x+1,y,y+1);}void set(C T&v,C int&x,C int&y){add(-get(x,y)+v,x,y);}};

// unionfind: unite(x,y), separate(x,y), same(x,y)=1(united),0(ununited,unseparated),-1(separated), diff(x,y)
TTT class unionfind{public:VI par;VT wei;int n;void resize(C int&_n){par.resize(_n);fou(i,n,_n)par[i]=i+1;n=_n;wei.resize(n);}int size(){return n;}constexpr unionfind(){}unionfind(C int&_n){resize(_n);}int get(C int&x){if(par[x]==(x+1))return x+1;else{int parent=sign(par[x])*get(abs(par[x])-1);wei[x]+=wei[abs(par[x])-1];return par[x]=parent;}}void get(int&x,int&y){x=get(x);y=get(y);}int weight(int x){get(x);return wei[x];}bool unite(int x,int y,T w=0){if(w<0){swap(x,y);w*=-1;}w+=weight(x)-weight(y);get(x,y);if(abs(x)==abs(y))return w==weight(x)-weight(y);par[abs(y)-1]=x*sign(y);wei[abs(y)-1]=w;return true;}bool separate(int x,int y,T w=0){if(w<0){swap(x,y);w*=-1;}w+=weight(x)-weight(y);get(x,y);if(abs(x)==abs(y))return w==weight(x)-weight(y);par[abs(y)-1]=-x*sign(y);wei[abs(y)-1]=w;return true;}
int same(int x,int y){get(x,y);return (x==y)-(x==-y);}T diff(int x,int y){return weight(y)-weight(x);}};


TTT int binarysearch(C VT&vec,C T&n,C char&c){
  int left=-1;int right=si(vec);bool outleft=(c=='d'||c=='U');int getmin=(c=='u'||c=='d')?1:-1;
  while(right-left>1){
    int mid=left+(right-left)/2;
    bool OK=vec[mid]*getmin>=n*getmin;
    if(outleft==OK) left=mid;else right=mid;
  }
  return outleft?left:right;
}

// graph class
TTT class graph{public:struct edge{int from,to;T cost;};int n;vector<edge>edges;VVI g;function<bool(int)>ignore;graph(int _n):n(_n){g.resize(n);ignore=nullptr;}virtual int add(int from,int to,T cost)=0;void clear(){edges.clear();g=VVI(n);};virtual void set_ignore_edge_rule(C function<bool(int)> &f){ignore=f;}virtual void reset_ignore_edge_rule(){ignore = nullptr;}};
// undigraph: add(x,y,w), lowling() -> VI bridge, crunode;
TTT class undigraph:public GT{public:VI bridge,crunode;using GT::edges;using GT::g;using GT::n;undigraph(int n):GT(n){}int add(int from,int to,T cost=1){int id=si(edges);g[from].emb(id);g[to].emb(id);edges.pub({from,to,cost});return id;}int lowlink_dfs(C int&from,int&k,C int&par,VI&ord,VI&low){ord[from]=low[from]=k++;bool is_crunode=false;int cnt=0;tra(id,g[from]){auto e=edges[id];int to=e.from^e.to^from;if(ord[to]==-1){cnt++;k=lowlink_dfs(to,k,from,ord,low);chmin(low[from],low[to]);is_crunode|=par!=-1&&low[to]>=ord[from];if(ord[from]<low[to])bridge.emb(id);}elif(to!=par)chmin(low[from],ord[to]);}is_crunode|=par==-1&&cnt>1;if(is_crunode)crunode.emb(from);return k;}void lowlink(){VI ord(n,-1),low(n,-1);int k=0;fou(i,0,n)if(ord[i]==-1)k=lowlink_dfs(i,k,-1,ord,low);}};
// digraph: add(x,y,w), scc() -> VVI strong_nodes, strong_edgeid; VI weak_edgeid, strong_id;, topsort()
TTT class digraph:public GT{public:VVI strong_nodes,strong_edgeid;VI weak_edgeid,strong_id;using GT::edges;using GT::g;using GT::n;using GT::ignore;digraph(int _n):GT(_n){}int add(int from,int to,T cost=1){int id=si(edges);g[from].emb(id);edges.pub({from,to,cost});return id;}DGT revgraph()const{DGT rev(n);tra(e,edges)rev.add(e.to,e.from,e.cost);if(ignore!=nullptr){rev.set_ignore_edge_rule([&](int id){return ignore(id);});}return rev;}void scc_dfs(C int&from,VI&ord,VI&strong_id){if(strong_id[from]==-1)return;strong_id[from]=-1;tra(id,g[from]){auto&e=edges[id];int to=e.from^e.to^from;scc_dfs(to,ord,strong_id);}ord.emb(from);}void scc_rdfs(C int&from,C int&cnt,C DGT&revg){if(strong_id[from]!=-1)return;strong_id[from]=cnt;strong_nodes[cnt].emb(from);tra(id,revg.g[from]){auto&e=revg.edges[id];int to=e.from^e.to^from;if(strong_id[to]==-1){strong_edgeid[cnt].emb(id);scc_rdfs(to,cnt,revg);}elif(strong_id[to]==cnt){strong_edgeid[cnt].emb(id);}elif(strong_id[to]<cnt)weak_edgeid.emb(id);}}void scc(){VI ord(n);strong_id.assign(n,0);DGT revg=revgraph();fou(from,0,n)scc_dfs(from,ord,strong_id);reverse(all(ord));int cnt=0;tra(from,ord)if(strong_id[from]==-1){strong_edgeid.emb();strong_nodes.emb();scc_rdfs(from,cnt,revg);cnt++;}}
VI topsort(){VI to_node(n,0);fou(id,0,si(edges)){if(ignore!=nullptr&&ignore(id))continue;to_node[edges[id].to]++;}VI sorted;fou(i,0,n)if(to_node[i]==0)sorted.emb(i);fou(i,0,si(sorted)){tra(id,g[sorted[i]]){if(ignore!=nullptr&&ignore(id))continue;if(--to_node[edges[id].to]==0)sorted.emb(edges[id].to);}}if(si(sorted)!=n)return VI();return sorted;}};
// forest <- undigraph: search_all() -> VI root, parent, treesize, depth, dist;
TTT class forest:public UGT{public:using UGT::edges;using UGT::g;using UGT::n;VI root,parent,treesize,depth,dist;forest(int _n):UGT(_n){init();}void init(){root=VI(n,-1);parent=VI(n,-1);treesize=VI(n,1);depth=VI(n);dist=VT(n);}void dfs(int from){tra(id,g[from]){auto&e=edges[id];if(from==e.from)swap(e.from,e.to);int to=e.from;if(root[to]!=-1){continue;}parent[to]=from;dist[to]=dist[from]+e.cost;depth[to]=depth[from]+1;root[to]=root[from];dfs(to);treesize[from]+=treesize[to];}}void set_root(int v){if(root[v]==-1){root[v]=v;dfs(v);}}void search_all(){init();fou(v,0,n)set_root(v);}T diameter(){init();set_root(0);int end=argmax(dist);init();set_root(end);return dist[argmax(dist)];}};
// shortest path
TTT VT dijkstra(C GT&g,C int&start){VT d(g.n,MAXT);priority_queue<pair<T,int>,vector<pair<T,int>>,greater<pair<T,int>>>s;d[start]=0;s.emplace(d[start],start);while(!s.empty()){T expected=s.top().fi;int from=s.top().se;s.pop();if(d[from]<expected)continue;tra(id,g.g[from]){auto e=g.edges[id];int to=e.from^e.to^from;if(d[from]+e.cost<d[to]){d[to]=d[from]+e.cost;s.emplace(d[to],to);}}}return d;}
TTT VT bellman(C GT&g,C int&start){VT d(g.n,MAXT);d[start]=0;bool nl=false;VI upd_b(1,start);int i=0;while(!upd_b.empty()&&i<g.n){VI upd;tra(from,upd_b){tra(id,g.g[from]){auto e=g.edges[id];int to=e.from^e.to^from;if(d[to]>d[from]+e.cost){d[to]=d[from]+e.cost;upd.emb(to);if(i==g.n-1){d[to]=MINT;nl=true;}}}}upd_b=upd;i++;}if(nl){i=0;while(!upd_b.empty()&&i<g.n-1){VI upd;tra(from,upd_b){tra(id,g.g[from]){auto e=g.edges[id];int to=e.from^e.to^from;if(d[to]!=MINT){d[to]=MINT;upd.emb(to);}}}upd_b=upd;i++;}}return d;}
// matrix graph
TTT VVT matrix_of_graph(C GT&g){VVT mat(g.n,VT(g.n,MAXT));fou(from,0,g.n){tra(id,g.g[from]){auto e=g.edges[id];int to=e.from^e.to^from;mat[from][to]=e.cost;}}return mat;}
TTT VVT warshall(VVT mat){fou(k,0,si(mat))fou(i,0,si(mat))fou(j,0,si(mat))mat[i][j]=(mat[i][k]==MAXT||mat[k][j]==MAXT)?mat[i][j]:min(mat[i][j],mat[i][k]+mat[k][j]);return mat;}
// minimum spanning tree
TTT T prim(GT&g,C int&start=0){vector<pair<T,PII>> cft;T total = 0;VI used(g.n,0);priority_queue<pair<T,PII>,vector<pair<T,PII>>,greater<pair<T,PII>>>s;s.emplace(mp(0,-1,start));while(!s.empty()){T cost=s.top().fi;int fromfrom=s.top().se.fi,from=s.top().se.se;s.pop();if(used[from])continue;used[from]=true;total+=cost;if(fromfrom>-1)cft.emb(mp(cost,fromfrom,from));tra(id,g.g[from]){auto e=g.edges[id];int to=e.from^e.to^from;s.emplace(mp(e.cost,from,to));}}g.clear();tra(a,cft)g.add(a.sec,a.thi,a.fir);return total;}

TTT class flow_graph {public:
  static constexpr T eps = (T) 1e-9;struct edge{int from;int to;T c;T f;};
  VVI g;vector<edge> edges;VI ptr,d,q;int n,start,goal;T flow;
  flow_graph(int _n,int start,int goal):n(_n),start(start),goal(goal),flow(0){g.resize(_n);ptr.resize(_n);d.resize(_n);q.resize(_n);}
  void clear_flow(){tra(e,edges)e.f=0;flow=0;}
  void add(C int&from,C int&to,C T&for_cap,C T&back_cap){g[from].emb(si(edges));edges.pub({from,to,for_cap,0});g[to].emb(si(edges));edges.pub({to,from,back_cap,0});}
  bool expath(){fill(all(d),-1);q[0]=goal;d[goal]=0;int beg=0,end=1;while(beg<end) {int i=q[beg++];tra(id,g[i]){edge&e=edges[id];if(edges[id^1].c-edges[id^1].f>eps&&d[e.to]==-1) {d[e.to]=d[i]+1;if(e.to==start)return true;q[end++]=e.to;}}}return false;}
  T dfs(int v,T w){if(v==goal)return w;while(ptr[v]>=0){int id=g[v][ptr[v]];edge&e=edges[id];if(e.c-e.f>eps&&d[e.to]==d[v]-1){T t=dfs(e.to,min(e.c-e.f,w));if(t>eps){edges[id].f+=t;edges[id^1].f-=t;return t;}}ptr[v]--;}return 0;}
  T max_flow(){while(expath()){fou(i,0,n)ptr[i]=si(g[i])-1;T big_add=0;while(true){T add=dfs(start,MAXT);if(add<=eps)break;big_add+=add;}if(big_add<=eps)break;flow+=big_add;}return flow;}
  vector<bool>min_cut(){max_flow();vector<bool>ret(n);fou(i,0,n)ret[i]=(d[i]!=-1);return ret;}
};









class time_measure{chrono::system_clock::time_point start;public:time_measure():start(chrono::system_clock::now()){}void time(){out(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-start).count(),"ms");}};
// ifstream ifs("ifile.txt");ifs>>a;
// ofstream ofs("ofile.txt");ofs<<a;




TTT VI euler_path(C UGT&ug){
  VI euler,used(ug.n);VVI g_next(ug.n,VI(ug.n));
  int beg=-1,end=-1;
  fou(i,0,ug.n){
    out(ug.g[i]);
    if(si(ug.g[i])%2==1){
      if(beg==-1) beg=i;
      elif(end==-1) end=i;
      else return VI(0);
    }
  }
  // 閉路を作れたら追加，それまでバッファ
  euler.emb(beg);
  euler.emb(end);
  return euler;
}

// undigraph<int> ug(6);
//   fou(i,0,5)ug.add(i,i+1);
//   ug.add(5,0);
//   ug.add(3,5);
//   ug.add(5,2);
//   ug.add(2,0);
//   // ug.add(1,4);
//   out(euler_path(ug));
