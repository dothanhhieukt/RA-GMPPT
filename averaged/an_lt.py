import numpy as np, pandas as pd, json, glob, os, functools
from plant import PVString
D=os.path.join(os.environ.get('RA_DATA','../data'),'one_hour','variable','')
@functools.lru_cache(maxsize=None)
def curve(G,T):
    pv=PVString(list(G),T,n=3000); V=np.linspace(1,pv.V.max(),3000); return V, V*pv.i(V)
def reach(df):
    Er=0; Eg=0
    for r in df.itertuples():
        G=tuple(round(x/10)*10 for x in (r.G1,r.G2,r.G3,r.G4)); V,P=curve(G,round(r.T*2)/2)
        Ib=r.Pb/max(r.Vb,1); vb0=(r.Vb+0.025)/0.95+0.5+0.02*Ib
        m=V>=vb0; pb=float(np.interp(vb0,V,P)); Pr=max(P[m].max() if m.any() else 0,pb)
        cap=10*r.Vb/0.97
        Er+=min(Pr,cap); Eg+=min(P.max(),cap)
    return Er,Eg
S=[json.loads(l) for l in open(D+'summary.jsonl')]
rows=[]
for s in S:
    f=D+f"a_{s['alg']}_{s['abl']}_{s['seed']}.csv"
    df=pd.read_csv(f); Er,Eg=reach(df)
    s.update(Ereach=Er,Egmpp=Eg,etaR=s['Epv']/Er); rows.append(s)
R=pd.DataFrame(rows)
nm={(5,0):'DPSO',(5,8):'DPSO*',(6,0):'RA-DPSO',(6,8):'RA-DPSO*',(6,64):'RA-DPSO-R1',(7,0):'CS',(8,0):'RA-CS',(8,64):'RA-CS-R1',(6,2):'RA-DPSO-F',(6,10):'RA-DPSO-F*',(8,2):'RA-CS-F'}
R['name']=[nm[(a,b)] for a,b in zip(R.alg,R.abl)]
R.to_csv('lt_results.csv',index=False)
pd.set_option('display.width',200)
print(R.groupby('name')[['E','Epv','etaR','tUnr','tBnd','nScan','nR2','nOver','nOver50','Ipk','socEnd']].agg(['mean','min','max']).round(4).T.to_string())
