import numpy as np, pandas as pd, sys, multiprocessing as mp
from plant import PVString, SOCtab, OCVtab
from ctrl import DPSO, CS, SCAN
VB=200.0
def sample(rng):
    G=np.where(rng.random(4)<0.35,1000.0,rng.uniform(150,950,4)); G=np.sort(G)[::-1]
    T=rng.uniform(25,65); soc=rng.uniform(0.1,0.95); pack=int(rng.choice([48,72]))
    return G,T,soc,pack
def evaluate(args):
    idx,seed=args
    rng=np.random.default_rng(seed)
    G,T,soc,pack=sample(rng)
    pv=PVString(G,T,n=3000)
    Ns=13 if pack==48 else 20; R=Ns*0.045/6
    ocv=Ns*np.interp(soc,SOCtab,OCVtab)
    V=np.linspace(1,pv.V.max(),6000); P=V*pv.i(V)
    def vbnd(Ppv):
        Ib=0.99*Ppv/ocv; Vb=ocv+R*Ib; return (Vb+0.025)/0.95+0.5+0.02*Ib, Vb
    # physical boundary at the boundary power (fixed point)
    vb0,_=vbnd(300.0)
    for _ in range(3):
        Pb=float(np.interp(vb0,V,P)); vb0,Vbat=vbnd(Pb)
    m=V>=vb0
    Preach=max(P[m].max(),Pb); Vreach=V[m][P[m].argmax()] if P[m].max()>=Pb else vb0
    ig=P.argmax(); Pg=P[ig]; Vg=V[ig]
    if Vg>=vb0: sit='S1'
    elif Pb>=P[m].max()-0.5: sit='S3'
    else: sit='S2'
    cc = Preach*0.96 > 1.01*Vbat*10
    out=dict(idx=idx,pack=pack,soc=soc,T=T,G=','.join('%d'%g for g in G),Pgmpp=Pg,Vg=Vg,Preach=Preach,Vbnd=vb0,Pbnd=Pb,sit=sit,cc=cc,
             nVmp_Vbat=(Vg/Vbat))
    if cc: return out
    # static closed loop of each core: one call = one evaluation period (TS = 2.6 ms)
    TS=2.6e-3; xmaxC=1-(Vbat/0.95-0.5)/VB
    def pw(Vref):
        Va=max(Vref,vb0); return Va, float(pv.i(Va))
    nz=np.random.default_rng(seed+99)
    for nm,cls,ra,ar in [('DPSO',DPSO,0,0),('RA-DPSO',DPSO,1,0),('DPSO-AR',DPSO,0,1),('CS',CS,0,0),('RA-CS',CS,1,0),('CS-AR',CS,0,1),
                         ('SCAN',SCAN,0,0),('RA-SCAN',SCAN,1,0),('SCAN-AR',SCAN,0,1)]:
        core=cls(TS) if cls is not CS else CS(TS,seed%1000+1)
        opt={'AR':bool(ar),'dVR':2.0,'xR2':xmaxC}
        x=0.4; acc=0; n=0; Pl=0
        for k in range(300):
            Va,Ia=pw(VB*(1-x))
            Vm=Va+0.05*nz.standard_normal(); Im=Ia+0.005*nz.standard_normal()
            x,ph,rs,xh=core.step(Vm,Im,bool(ra),xmaxC,False,True,True,(xmaxC if ra else 0.4),opt)
            Va,Ia=pw(VB*(1-x)); Pl=Va*Ia
            if k>=8: acc+=Pl; n+=1
        out[nm+'_eff']=Pl/Preach; out[nm+'_etaE']=acc/n/Preach
    return out
if __name__=='__main__':
    N=int(sys.argv[1]); 
    with mp.Pool(2) as p: R=p.map(evaluate,[(i,10000+i) for i in range(N)],chunksize=20)
    D=pd.DataFrame(R); D.to_csv('mc/mc_static.csv',index=False)
    print(D.sit.value_counts()); print(D[~D.cc].groupby('sit')[[c for c in D.columns if c.endswith('_eff')]].mean().round(3))
