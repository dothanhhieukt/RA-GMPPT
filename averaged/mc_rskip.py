import numpy as np, pandas as pd, sys, multiprocessing as mp
import mc_static as M
from plant import PVString, SOCtab, OCVtab
from ctrl import DPSO
VB=200.0
def ev(args):
    idx,seed=args
    rng=np.random.default_rng(seed); G,T,soc,pack=M.sample(rng)
    pv=PVString(G,T,n=3000); Ns=13 if pack==48 else 20; R=Ns*0.045/6
    ocv=Ns*np.interp(soc,SOCtab,OCVtab); V=np.linspace(1,pv.V.max(),6000); P=V*pv.i(V)
    def vbnd(Ppv):
        Ib=0.99*Ppv/ocv; Vb=ocv+R*Ib; return (Vb+0.025)/0.95+0.5+0.02*Ib, Vb
    vb0,_=vbnd(300.0)
    for _ in range(3):
        Pb=float(np.interp(vb0,V,P)); vb0,Vbat=vbnd(Pb)
    m=V>=vb0; Preach=max(P[m].max(),Pb)
    cc = Preach*0.96 > 1.01*Vbat*10
    out=dict(idx=idx)
    if cc: return out
    TS=2.6e-3; xmaxC=1-(Vbat/0.95-0.5)/VB
    def pw(Vref):
        Va=max(Vref,vb0); return Va, float(pv.i(Va))
    nz=np.random.default_rng(seed+99)
    for nm,ra,ar in [('DPSO',0,0),('RA-DPSO',1,0),('DPSO-AR',0,1)]:
        nz=np.random.default_rng(seed+99) if nm=='DPSO' else nz
        core=DPSO(TS); opt={'AR':bool(ar),'dVR':2.0,'xR2':xmaxC,'RSkip':True}
        x=0.4
        for k in range(300):
            Va,Ia=pw(VB*(1-x)); Vm=Va+0.05*nz.standard_normal(); Im=Ia+0.005*nz.standard_normal()
            x,ph,rs,xh=core.step(Vm,Im,bool(ra),xmaxC,False,True,True,(xmaxC if ra else 0.4),opt)
            Va,Ia=pw(VB*(1-x)); Pl=Va*Ia
        out[nm+'*_eff']=Pl/Preach
    return out
if __name__=='__main__':
    N=int(sys.argv[1])
    with mp.Pool(1) as p: R=p.map(ev,[(i,10000+i) for i in range(N)],chunksize=50)
    pd.DataFrame(R).to_csv('mc/mc_rskip.csv',index=False)
