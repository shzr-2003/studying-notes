// luogu-judger-enable-o2
# include <cstdio>
# include <iostream>
# include <cstring>
# define R register int
# define LL long long
# define nl (n<<1)
# define nr (n<<1|1)

using namespace std;

const int maxn=100005;
const LL inf=5e10;
int n,m,x,y,a,b,firs[maxn],h;
int id[maxn],Top[maxn],low[maxn],siz[maxn],son[maxn],cnt,fid[maxn],f[maxn];
char typ[10];
LL dp[maxn][2],t1,t2,v[maxn],fans;
struct edge { int too,nex; }g[maxn<<1];
struct mat 
{
    LL a[2][2];
    void write()
    {
        for (R i=0;i<=1;++i)
        {
            for (R j=0;j<=1;++j) printf("%lld ",a[i][j]);
            printf("\n");
        }
        printf("\n");
    }
}t[maxn<<2],s[maxn],ans;

int read();
void add (int x,int y); 
void dfs1 (int x);
void dfs2 (int x,int Tp);
void Dp (int x);
void build (int n,int l,int r);
mat ask_ans (int x);
mat ask (int n,int l,int r,int ll,int rr);
void change (int n,int l,int r,int pos);
void mag (int x,LL y); //magic
mat operator * (mat a,mat b);

int main()
{
    scanf("%d%d%s",&n,&m,typ);
    for (R i=1;i<=n;++i) v[i]=read();
    for (R i=1;i<n;++i)
    {
        x=read(),y=read();
        add(x,y),add(y,x);
    }
    dfs1(1); dfs2(1,1); Dp(1); 
    build(1,1,n);
    for (R i=1;i<=m;++i)
    {
        fans=0;
        a=read(),x=read(),b=read(),y=read();
        t1=v[a],t2=v[b];
        if(x==0&&y==0&&(f[a]==b||f[b]==a)) { printf("-1\n"); continue; }
        if(x==0) { mag(a,inf); } 
        if(x==1) { mag(a,-inf); fans+=inf+t1; }
        if(y==0) { mag(b,inf); }
        if(y==1) { mag(b,-inf); fans+=inf+t2; }
        ans=ask_ans(1);
        printf("%lld\n",fans+min(min(ans.a[0][1],ans.a[1][0]),ans.a[1][1]));
        mag(a,t1),mag(b,t2);
    }
    return 0;
}

int read()
{
    R x=0,f=1;
    char c=getchar();
    while (!isdigit(c)) { if(c=='-') f=-f; c=getchar(); }
    while (isdigit(c)) x=(x<<3)+(x<<1)+(c^48),c=getchar();
    return x;
}

void add (int x,int y)
{
    g[++h].nex=firs[x];
    firs[x]=h;
    g[h].too=y;
}

void dfs1 (int x)
{
    siz[x]=1;
    int maxx=-1,j;
    for (R i=firs[x];i;i=g[i].nex)
    {
        j=g[i].too;
        if(j==f[x]) continue;
        f[j]=x;
        dfs1(j); siz[x]+=siz[j];
        if(siz[j]>=maxx) maxx=siz[j],son[x]=j;
    }
}

void dfs2 (int x,int Tp)
{
    id[x]=++cnt; Top[x]=Tp; low[Tp]=cnt; fid[cnt]=x;
    if(!son[x]) return;
    dfs2(son[x],Tp);
    int j;
    for (R i=firs[x];i;i=g[i].nex)
    {
        j=g[i].too;
        if(son[x]==j||f[x]==j) continue;
        dfs2(j,j);
    }
}

void Dp (int x)
{
    int j;
    for (R i=firs[x];i;i=g[i].nex)
    {
        j=g[i].too;
        if(j==f[x]) continue;
        Dp(j);
        dp[x][0]+=dp[j][1];
        dp[x][1]+=min(dp[j][0],dp[j][1]);
    }
    dp[x][1]+=v[x];
}

void build (int n,int l,int r)
{
    if(l==r)
    {
        int x=fid[l];
        t[n].a[0][0]=inf;
        t[n].a[0][1]=dp[x][0]-dp[ son[x] ][1];
        t[n].a[1][0]=t[n].a[1][1]=dp[x][1]-min(dp[ son[x] ][0],dp[ son[x] ][1]);
        s[l]=t[n];
        return;
    }
    int mid=(l+r)>>1;
    build(nl,l,mid);
    build(nr,mid+1,r);
    t[n]=t[nl]*t[nr];
}

mat ask_ans (int x) { return ask(1,1,n,id[ Top[x] ],low[ Top[x] ]); }

mat ask (int n,int l,int r,int ll,int rr)
{
    if(ll<=l&&r<=rr) return t[n];
    int mid=(l+r)>>1;
    if(rr<=mid) return ask(nl,l,mid,ll,rr);
    if(ll>mid) return ask(nr,mid+1,r,ll,rr);
    return ask(nl,l,mid,ll,rr)*ask(nr,mid+1,r,ll,rr);
}

mat operator * (mat a,mat b)
{
    mat c;
    c.a[0][0]=min(a.a[0][0]+b.a[0][0],a.a[0][1]+b.a[1][0]);
    c.a[0][1]=min(a.a[0][0]+b.a[0][1],a.a[0][1]+b.a[1][1]);
    c.a[1][0]=min(a.a[1][0]+b.a[0][0],a.a[1][1]+b.a[1][0]);
    c.a[1][1]=min(a.a[1][0]+b.a[0][1],a.a[1][1]+b.a[1][1]);
    return c;
}

void chan (int n,int l,int r,int pos)
{
    if(l==r) { t[n]=s[l]; return; }
    int mid=(l+r)>>1;
    if(pos<=mid) chan(nl,l,mid,pos);
    if(pos>mid) chan(nr,mid+1,r,pos);
    t[n]=t[nl]*t[nr];
}

void mag (int x,LL y)
{
    s[ id[x] ].a[1][0]+=y-v[x]; 
    s[ id[x] ].a[1][1]+=y-v[x];
    v[x]=y;
    mat t1,t2;
    while(x)
    {
        t1=ask_ans(Top[x]);
        chan(1,1,n,id[x]);
        t2=ask_ans(Top[x]);
        x=f[ Top[x] ];
        if(x)
        {
            s[ id[x] ].a[0][1]+=t2.a[1][1]-t1.a[1][1];
            s[ id[x] ].a[1][0]+=min(t2.a[0][1],t2.a[1][1])-min(t1.a[0][1],t1.a[1][1]);
            s[ id[x] ].a[1][1]=s[ id[x] ].a[1][0];
        }
    }
}
