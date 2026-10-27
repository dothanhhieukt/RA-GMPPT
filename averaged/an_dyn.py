import pandas as pd, numpy as np
from scipy.stats import wilcoxon
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
D=pd.read_csv('camp_dyn_fixed.csv')
order=['72-C1>C4-S60','72-C2>C4-S60','72-C4>C2-S60','72-C4>C1-S60']
g=D.groupby(['scen','name']).agg(n=('Ebat','size'),eff=('eff','mean'),effmin=('eff','min'),E=('Ebat','mean'),Ipk=('Ipk','mean'),Ipkmax=('Ipk','max'),tOver=('tOver','mean'),tset=('tset','median'),Qover=('Qover','mean'))
print(g.round(3).to_string())
rows=[]
for sc in order:
    x=D[D.scen==sc]
    for a,b in [('RA-CS','CS'),('RA-CS-F','RA-CS')]:
        P=x[x.name.isin([a,b])].pivot_table(index='seed',columns='name',values='Ebat')
        d=100*(P[a]/P[b]-1); p=wilcoxon(d).pvalue if (d!=0).any() else 1
        rows.append((sc,a,b,d.median(),(d>0).mean(),p))
R=pd.DataFrame(rows,columns=['scen','a','b','dE','pos','p'])
for (a,b),grp in R.groupby(['a','b']):
    p=grp.p.values; o=np.argsort(p); n=len(p); adj=np.empty(n); run=0
    for k,i in enumerate(o): run=max(run,min(1,(n-k)*p[i])); adj[i]=run
    R.loc[grp.index,'pH']=adj
print(R.round(4))
