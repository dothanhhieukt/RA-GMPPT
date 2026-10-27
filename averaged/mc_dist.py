import numpy as np, pandas as pd, sys, multiprocessing as mp
from plant import PVString, SOCtab, OCVtab
def run(args):
    idx,seed,pfull,lo=args
    rng=np.random.default_rng(seed)
    G=np.where(rng.random(4)<pfull,1000.0,rng.uniform(lo,950,4)); G=np.sort(G)[::-1]
    T=rng.uniform(25,65); soc=rng.uniform(0.1,0.95); pack=int(rng.choice([48,72]))
    pv=PVString(G,T,n=3000); Ns=13 if pack==48 else 20; R=Ns*0.045/6
    ocv=Ns*np.interp(soc,SOCtab,OCVtab); V=np.linspace(1,pv.V.max(),6000); P=V*pv.i(V)
    def vbnd(Ppv):
        Ib=0.99*Ppv/ocv; Vb=ocv+R*Ib; return (Vb+0.025)/0.95+0.5+0.02*Ib, Vb
    vb0,_=vbnd(300.0)
    for _ in range(3):
        Pb=float(np.interp(vb0,V,P)); vb0,Vbat=vbnd(Pb)
    m=V>=vb0; Preach=max(P[m].max(),Pb); Vg=V[P.argmax()]
    sit='S1' if Vg>=vb0 else ('S3' if Pb>=P[m].max()-0.5 else 'S2')
    return dict(pfull=pfull,lo=lo,pack=pack,soc=soc,sit=sit,cc=Preach*0.96>1.01*Vbat*10)
if __name__=='__main__':
    J=[(i,20000+i,pf,lo) for pf,lo in [(0.0,100),(0.6,150),(0.35,150)] for i in range(1500)]
    with mp.Pool(2) as p: R=pd.DataFrame(p.map(run,J,chunksize=50))
    R.to_csv('mc/mc_dist.csv',index=False)
    X=R[~R.cc]
    print(X.groupby(['pfull','lo','pack']).sit.value_counts(normalize=True).unstack().round(3))
    print(X[(X.pack==72)&(X.soc>0.7)].groupby(['pfull','lo']).sit.value_counts(normalize=True).unstack().round(3))
