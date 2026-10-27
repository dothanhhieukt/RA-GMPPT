import numpy as np, scipy.io as sio, sim, sys
import os
TR=os.path.join(os.environ.get('RA_DATA','../data'),'results','traces')   # switched-model waveforms (results archive)
def met(d):
    t=d['t']; m=t>=t[-1]-0.05
    Pss=np.median(d['P'][m]); E=np.trapezoid(d['Vb']*d['Ib'],t)
    # 1ms avg Ipk
    dt=t[1]-t[0]; w=max(1,int(round(1e-3/dt))); Ib=np.convolve(d['Ib'],np.ones(w)/w,'same')
    return Pss, E, Ib.max(), np.median(d['V'][m])
jobs=[('static_DPSO_72-C4-S60_r01',5,0,72,'C4','C4',0.6,0.8),('static_RA-DPSO_72-C4-S60_r01',6,0,72,'C4','C4',0.6,0.8),
      ('static_DPSO_48-C1-S60_r01',5,0,48,'C1','C1',0.6,0.8),('static_RA-DPSO_48-C1-S60_r01',6,0,48,'C1','C1',0.6,0.8),
      ('static_CS_72-C4-S60_r01',7,0,72,'C4','C4',0.6,0.8),('static_RA-CS_72-C4-S60_r01',8,0,72,'C4','C4',0.6,0.8),
      ('dynamic_DPSO_72-C4_C1-S60_r01',5,0,72,'C4','C1',0.6,1.1),('dynamic_RA-DPSO_72-C4_C1-S60_r01',6,0,72,'C4','C1',0.6,1.1)]
for nm,alg,abl,pack,cA,cB,soc,T in jobs:
    S=sio.loadmat(os.path.join(TR,nm.split('_')[0],nm+'.mat'),squeeze_me=True)
    tsw=0.5 if cA!=cB else np.inf
    A=sim.run(alg,abl,1,pack,cA,cB,soc,tsw=tsw,Tend=T)
    ms=met(S); ma=met(A)
    print('%-36s SIM Pss %6.1f E %6.1f Ipk %5.2f V %6.1f | AVG Pss %6.1f E %6.1f Ipk %5.2f V %6.1f'%((nm,)+ms+ma))
    np.savez('val_'+nm+'.npz',**A)
