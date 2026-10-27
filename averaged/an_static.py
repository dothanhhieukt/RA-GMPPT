import pandas as pd, numpy as np
from scipy.stats import wilcoxon
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
D=pd.read_csv('camp_static.csv')
print(D.groupby('name').size())
order=['48-C4-S60','48-C7-S20','48-C7-S60','48-C8-S60','72-C4-S60','72-C6-S60','72-C7-S60','72-C7-S95','72-C8-S20','72-C8-S60']
rng=np.random.default_rng(0)
def boot_ci(d,B=5000):
    d=np.asarray(d); m=[np.median(rng.choice(d,len(d))) for _ in range(B)]; return np.percentile(m,[2.5,97.5])
rows=[]
for sc in order:
    g=D[D.scen==sc]; r={'scen':sc}
    for a in ['CS','RA-CS','CS-AR','RA-CS-noRA','RA-CS-R1','P&O','DPSO','DPSO-AR','RA-DPSO','SCAN','SCAN-AR','RA-SCAN']:
        x=g[g.name==a]
        r[a+'_eff']=x.eff.mean(); r[a+'_min']=x.eff.min(); r[a+'_succ']=x.success.mean(); r[a+'_tconv']=x.tconv.median()*1e3
        r[a+'_E']=x.Ebat.mean(); r[a+'_tunr']=x.tUnreach.mean()*1e3; r[a+'_Ipk']=x.Ipk.mean(); r[a+'_etaE']=x.etaE.mean()
    P=g[g.name.isin(['CS','RA-CS'])].pivot_table(index='seed',columns='name',values=['Ebat','eff','tUnreach','Ipk','etaE'])
    dE=100*(P[('Ebat','RA-CS')]/P[('Ebat','CS')]-1)
    r['dE_med']=dE.median(); lo,hi=boot_ci(dE); r['dE_lo']=lo; r['dE_hi']=hi; r['dE_pos']=(dE>0).mean()
    r['p_E']=wilcoxon(dE).pvalue
    for m in ['eff','tUnreach','Ipk']:
        d=P[(m,'RA-CS')]-P[(m,'CS')]; r['d_'+m]=d.median()
        r['p_'+m]=wilcoxon(d).pvalue if (d!=0).any() else 1.0
    rows.append(r)
A=pd.DataFrame(rows)
for m in ['E','eff','tUnreach','Ipk']:
    p=A['p_'+m].values; o=np.argsort(p); n=len(p); adj=np.empty(n); run=0
    for k,i in enumerate(o): run=max(run,min(1,(n-k)*p[i])); adj[i]=run
    A['pH_'+m]=adj
A.to_csv('an_static.csv',index=False)
print(A[['scen','CS_eff','CS_min','CS_succ','RA-CS_eff','RA-CS_min','RA-CS_succ','CS-AR_eff','CS-AR_min','dE_med','dE_lo','dE_hi','pH_E','d_eff','pH_eff','pH_tUnreach','d_Ipk','pH_Ipk']].round(3).to_string())
print(A[['scen','CS_tconv','RA-CS_tconv','CS_tunr','RA-CS_tunr','CS_etaE','RA-CS_etaE']].round(1).to_string())
print(A[['scen','P&O_eff','DPSO_eff','DPSO-AR_eff','RA-DPSO_eff','SCAN_eff','SCAN-AR_eff','RA-SCAN_eff','DPSO_E','DPSO-AR_E','RA-DPSO_E','SCAN_E','SCAN-AR_E','RA-SCAN_E','P&O_E']].round(3).to_string())
print('scenario-level Wilcoxon on dE_med:', wilcoxon(A.dE_med).pvalue, 'eff:', wilcoxon(A['RA-CS_eff']-A['CS_eff']).pvalue)
S6=['48-C7-S60','48-C8-S60','72-C4-S60','72-C7-S60','72-C7-S95','72-C8-S60']
for a in ['CS','RA-CS','RA-CS-R1','RA-CS-noRA','CS-AR','DPSO','DPSO-AR','RA-DPSO','SCAN','SCAN-AR','RA-SCAN']:
    x=D[(D.name==a)&D.scen.isin(S6)]; print('%-12s 6sc eff %.3f min %.3f tunr %.0f E %.1f | all10 eff %.3f min %.3f succ %.2f E %.1f'%(a,x.groupby('scen').eff.mean().mean(),x.eff.min(),x.groupby('scen').tUnreach.mean().mean()*1e3,x.groupby('scen').Ebat.mean().mean(),
          D[D.name==a].groupby('scen').eff.mean().mean(),D[D.name==a].eff.min(),D[D.name==a].success.mean(),D[D.name==a].groupby('scen').Ebat.mean().mean()))
