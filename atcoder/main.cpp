#include <bits/stdc++.h>
#define int int64_t
#define double long double
#define C const
#define INIT ios_base::sync_with_stdio(0);cin.tie(0);cout<<fixed<<setprecision(10)
using namespace std;C int MAX=numeric_limits<int>::max();C int MIN=numeric_limits<int>::min();C double EPS=1e-7;C int MOD=1e9+7;
// template
#define TTT template<typename T>
#define TT12 template<typename T1,typename T2>
#define TTTT template<typename T,typename TT>
#define TTHT template<typename Hd,typename ...Tl>
#define TT12T template<typename T1,typename T2, typename T>
#define VT vector<T>
#define VVT vector<VT>
#define PTT pair<T1,T2>
#define MT Mod<T>
#define VMT vector<MT>
#define VVMT vector<VMT>
#define GT graph<T>
#define DGT digraph<T>
#define UGT undigraph<T>
#define MAXT numeric_limits<T>::max()
#define MINT numeric_limits<T>::min()
#define emb emplace_back
#define pub push_back
#define pob pop_back
#define fou(i,a,n) for(int i=a;i<n;i++)
#define fod(i,a,n) for(int i=n-1;i>=a;i--)
#define tra(e,v) for(auto&e:v)
#define elif(c) else if(c)
#define si(v) ((int)(v).size())
#define all(v) v.begin(),v.end()
#define fi first
#define se second
#define fir fi
#define sec se.fi
#define thi se.se
#define firs fi
#define seco sec
#define thir thi.fi
#define four thi.se
#define secon sec
#define third thir
#define fourt four.fi
#define fifth four.se
// mod
#define md1(op) TTTT MT operator op(C MT&l,C TT&r){return l op MT(r);}TTTT MT operator op (C TT&l,C MT&r){return MT(l) op r;}
#define mdcp(op) TTT bool operator op(C MT&l,C MT&r){return (int)l op (int)r;}md1(op)
#define mdcl(op) TTT MT operator op(C MT&l,C MT&r){return MT((int)l op (int)r);}md1(op)
template<typename M>class Mod{public:int value;constexpr Mod():value(){}TTT Mod(C T&a){value=normalize(a);}TTT static int normalize(C T&a){return (a%mod()+mod())%mod();}constexpr explicit operator int()C{return value;}constexpr static int mod(){return M::value;}Mod inv()C{int a=value;int b=mod(),u=0,v=1;while(a>1){u-=b/a*v;b%=a;swap(a,b);swap(u,v);}return Mod(v);}Mod operator++(){return Mod(++value);}Mod operator--(){return Mod(--value);}Mod operator++(signed){return Mod(value++);}Mod operator--(signed){return Mod(value--); }Mod operator-()C{return Mod(-value);}};TTT ostream&operator<<(ostream&os,C MT&m){return os<<(int)m;}TTT istream&operator>>(istream&is,MT&m){is>>m.value;m.value=MT::normalize((int)m);return is;}TTT string to_string(C MT&m){return to_string((int)m);}mdcp(==)mdcp(!=)mdcp(>=)mdcp(<=)mdcp(>)mdcp(<)mdcl(+)mdcl(-)mdcl(*)mdcl(<<)mdcl(>>)mdcl(&)mdcl(|)mdcl(^)TTT MT operator/(C MT&l,C MT&r){return l*r.inv();}md1(/)using Mint=Mod<integral_constant<decay<decltype(MOD)>::type,MOD>>;
// typedef VVVPIDCS, make_pair/vector calculation
#define defall(xdef) xdef(int,I)xdef(double,D)xdef(char,C)xdef(string,S)xdef(bool,B)xdef(Mint,M)
#define vdef(t,T) typedef vector<t>V##T;typedef vector<V##T>VV##T;typedef vector<VV##T>VVV##T;
#define pdef(T) T;vdef(T,T)
#define pdef0(t1,T1,t2,T2,t3,T3,t4,T4,t5,T5) typedef pair<t1,t2>pdef(P##T1##T2)typedef pair<t1,pair<t2,t3>>pdef(P##T1##T2##T3)typedef pair<t1,pair<t2,pair<t3,t4>>>pdef(P##T1##T2##T3##T4)typedef pair<t1,pair<t2,pair<t3,pair<t4,t5>>>>pdef(P##T1##T2##T3##T4##T5);
#define pdef1(...) pdef0(int,I,__VA_ARGS__)pdef0(double,D,__VA_ARGS__)pdef0(char,C,__VA_ARGS__)pdef0(string,S,__VA_ARGS__)pdef0(bool,B,__VA_ARGS__)pdef0(Mint,M,__VA_ARGS__)
#define pdef2(...) pdef1(int,I,__VA_ARGS__)pdef1(double,D,__VA_ARGS__)pdef1(char,C,__VA_ARGS__)pdef1(string,S,__VA_ARGS__)pdef1(bool,B,__VA_ARGS__)pdef1(Mint,M,__VA_ARGS__)
#define pdef3(...) pdef2(int,I,__VA_ARGS__)pdef2(double,D,__VA_ARGS__)pdef2(char,C,__VA_ARGS__)pdef2(string,S,__VA_ARGS__)pdef2(bool,B,__VA_ARGS__)pdef2(Mint,M,__VA_ARGS__)
#define pdef4(...) pdef3(int,I,__VA_ARGS__)pdef3(double,D,__VA_ARGS__)pdef3(char,C,__VA_ARGS__)pdef3(string,S,__VA_ARGS__)pdef3(bool,B,__VA_ARGS__)pdef3(Mint,M,__VA_ARGS__)
#define opall(calc) calc(+)calc(-)calc(*)calc(/)calc(%)calc(<<)calc(>>)calc(&)calc(|)calc(^)
#define eq_pair_vec(op) TT12 T1 operator op##=(T1&l,C T2&r){return l=l op r;}TT12T PTT operator op(C PTT&p,C T&c){return mp(p.fi op c,p.se op c);}TT12T PTT operator op(C T&c,C PTT&p){return mp(c op p.fi,c op p.se);}TT12 PTT operator op(C PTT l,C PTT&r){return mp(l.fi op r.fi,l.se op r.se);}TTTT VT operator op(VT vec,C TT&c){tra(v,vec)v=v op c;return vec;}TTTT VT operator op(C TT&c,VT vec){tra(v,vec)v=c op v;return vec;}TTT VT operator op(VT l,C VT&r){fou(i,0,min(si(l),si(r)))l[i]=l[i] op r[i];return l;}
defall(vdef)defall(pdef4)opall(eq_pair_vec)TT12 istream&operator>>(istream&is,PTT&p){return is>>p.fi>>p.se;}TTT istream&operator>>(istream&is,VT&v){tra(e,v)is>>e;return is;}TT12 ostream&operator<<(ostream&os,C PTT&p){return os<<p.fi<<" "<<p.se;}TTT ostream&operator<<(ostream&os,C VT&v){if(v.empty())return os;fou(i,0,si(v)-1)os<<v[i]<<" ";return os<<v[si(v)-1];}TT12 PTT mp(C T1&a,C T2&b){return make_pair(a,b);}TTHT auto mp(C Hd&hd,C Tl&...tl){return mp(hd,mp(tl...));}
// iostream, argmax, argmin, chmax, chmin
#define in(T,...) T __VA_ARGS__;_in(__VA_ARGS__)
void _in(){}TTHT void _in(Hd&hd,Tl&&...tl){cin>>hd;_in(forward<Tl>(tl)...);}TTT void out(C T&a){cout<<a<<"\n";}TTHT void out(C Hd&hd,C Tl&...tl){cout<<hd<<" ";out(tl...);}TTT void vout(C VT&v){tra(e,v)out(e);}TTT int argmax(C T&v){int m=0;fou(i,1,si(v))if(v[i]>v[m])m=i;return m;}TTT int argmin(C T&v){int m=0;fou(i,1,si(v))if(v[i]<v[m])m=i;return m;}TT12 void chmax(T1&a,C T2&b){a=max(a,(T1)b);}TT12 void chmin(T1&a,C T2&b){a=min(a,(T1)b);}
// math basic -> sign, gcd, lcm, fact, parm, comb, homo, power, is_pow
#define IT typename conditional<is_same<T,signed>::value,int,T>::type
TTT T sign(T a){return (a>0)-(a<0);}int _gcd(C int&a,C int&b){return b?_gcd(b,a%b):a;}int gcd(C int&a,C int&b){return _gcd(max(a,b),min(a,b));}int lcm(C int&a,C int&b){return a/gcd(a,b)*b;}TTT IT fact(C T&n){IT r=1;fou(i,1,(int)n+1)r*=i;return r;}TTTT IT parm(C T&m,C TT&n){IT r=1;fou(i,0,(int)n)r*=m-i;return r;}TTT int comb(C int&m,T n){int r=1;chmin(n,m-n);fou(i,0,(int)n){r*=m-i;r/=i+1;}return r;}TTTT MT comb(C MT&m,TT n){MT p=1,q=1,k=min(MT(n),m-n);fou(i,0,(int)n){p*=m-i;q*=i+1;}return p/q;}TTTT IT homo(C T&m,C TT&n){return comb(IT(m+n-1),min(IT(m-1),IT(n)));}TTTT IT power(T a,TT b){IT r=1;while(b>0){if((b&1)==1)r=r*a;a=a*a;b>>=1;}return r;}bool is_pow(C int&n,C int&m){int r;if(m==2)r=sqrt(n);if(m==3)r=cbrt(n);if(m==4)r=sqrt(sqrt(n));else r=pow(n,(double)1/m);return n==power(r,m);}VM nfact(1,1),ifact(1,1);int fact_max=0;void def_fact(C int&n){if(fact_max>=n)return;nfact.resize((int)n+1);fou(i,fact_max+1,(int)n+1)nfact[i]=i*nfact[i-1];ifact.resize((int)n+1);ifact[n]=nfact[n].inv();fod(i,fact_max+1,(int)n)ifact[i]=(i+1)*ifact[i+1];fact_max=n;}TTT Mint facts(C T&n){if(n>fact_max)def_fact((int)n);return nfact[(int)n];}TT12 Mint parms(C T1&m,C T2&n){if(m<n)return Mint(0);if(m>fact_max)def_fact((int)m);return nfact[(int)m]*ifact[(int)(m-n)];}TT12 Mint combs(C T1&m,C T2&n){if(m<n)return Mint(0);if(m>fact_max)def_fact((int)m);return nfact[(int)m]*ifact[(int)(m-n)]*ifact[(int)n];}TT12 Mint homos(C T1&m,C T2&n){return combs(m+n-1,min((int)m-1,(int)n));}TTTT VMT powers(C MT&a,C TT&b){VMT r((int)b+1,1);fou(i,1,(int)b+1)r[i]=a*r[i-1];return r;}
// cumsum, partsum, decumsum // sum, partof, dot, power // sortup, sortdown, bucket_sort // erase_unique, rotation
#define sortcomp(v,compXY) function<bool(decltype(v[0]),decltype(v[0]))>([&](auto&X,auto&Y){return compXY;})
TTT void cumsum(VT&v){fou(i,1,si(v))v[i]+=v[i-1];}TTT void cumsum(VVT&m){fou(i,1,si(m))m[i]+=m[i-1];fou(i,0,si(m))cumsum(m[i]);}TTT T partsum(VT&v,C int&l,C int&r){return v[r-1]-(l==0?0:v[l-1]);}TTT T partsum(VVT&m,C int&l1,C int&r1,C int&l2,C int&r2){return m[r1-1][r2-1]-(l1==0?0:m[l1-1][r2-1])-(l2==0?0:m[r1-1][l2-1])+(l1==0||l2==0?0:m[l1-1][l2-1]);}TTT void decumsum(VT&v){fod(i,1,si(v))v[i]-=v[i-1];}TTT void decumsum(VVT&m){fod(i,1,si(m))m[i]-=m[i-1];fod(i,0,si(m))decumsum(m[i]);}TTT T sum(C VT&v,C int&l=0,int r=-1){if(r==-1)r=si(v);T s=0;fou(i,l,r)s+=v[i];return s;}TTT T partof(T v){return v;}TTHT Hd partof(Hd v,C int&l,C int&r,C Tl&...tl){copy(v.begin()+l,v.begin()+r,v.begin());v.erase(v.begin()+r-l,v.end());fou(i,0,r-l)v[i]=partof(v[i],tl...);return v;}TTT VVT dot(C VVT&m1,C VVT&m2){VVT m(si(m1),VI(si(m2[0])));fou(i,0,si(m1))fou(j,0,si(m2[0]))fou(k,0,min(si(m1[0]),si(m2)))m[i][j]+=m1[i][k]*m2[k][j];return m;}TTT VVT power(VVT m,int n){VVT r(si(m),VT(si(m)));fou(i,0,si(m))r[i][i]=1;while(n>0){if(n&1)r=dot(r,m);m=dot(m,m);n>>=1;}return r;}TTT void sortup(VT&v,function<bool(T&,T&)>f){sort(all(v),f);}TTT void sortdown(VT&v,function<bool(T&,T&)>f){sort(v.rbegin(),v.rend(),f);}TTT void sortup(T&v){sort(all(v));}TTT void sortdown(T&v){sort(v.rbegin(),v.rend());}void bucket_sort(VI&v,C int&m,C int&M){VI b(M-m+1);tra(e,v)b[e-m]++;cumsum(b);fou(i,0,si(b))fou(j,i==0?0:b[i-1],b[i])v[j]=m+i;}TTT void erase_unique(T&v){v.erase(unique(all(v)),v.end());}TTT void rotation(T&v,C int&n){if(n>0)rotate(v.begin(),v.begin()+n,v.end());if(n<0)rotate(v.rbegin(),v.rbegin()-n,v.rend());}


////////// additional library //////////






// 10^5 -> NlogN, 3000 -> N^2, 200 -> N^3, 50 -> N^4, 20 -> 2^N 






signed main(){INIT;
  // unionfind<int> uf;
  // uf.resize(10);
  // out(uf.unite(0,1,-2));
  // out(uf.unite(1,2,-3));
  // out(uf.unite(0,2,-5));
  // out(uf.unite(0,2,-6));
  // fou(i,0,10) uf.get(i);
  // out(uf.par);
  // out(uf.wei);
  // out(uf.diff(1,2));
  VI a{1,2,3,4,5};
  rotate(a.begin(),a.begin()+1,a.end());
  out(a);

}