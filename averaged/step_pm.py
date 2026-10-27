import numpy as np
exec(open('bode_tol.py').read().split('nom=(')[0])
import plant
def step(o,p,KP=12,N=200):
    Phi,Gam,Cav,Dav=disc(o,*p); V=o['V']; Dff=o['Vb']/V
    x=np.zeros(3); xi=0.0; u=0.0; y=[]
    for k in range(N):
        yk=float((Cav@x)[0]+Dav[0,0]*u)       # period-averaged deviation of v_pv
        e=(yk-1.0)/V                     # unit step of the reference (1 V)
        xi+=KI*T*e
        unew=Dff*(KP*e+xi)
        x=Phi@x+Gam[:,0]*u; u=unew; y.append(yk)
    y=np.array(y); t=np.arange(N)*T
    ov=(y.max()-1)*100; out=np.where(np.abs(y-1)>0.05)[0]; ts=t[out[-1]+1] if len(out) else 0
    return ov, ts*1e3
nom=(plant.L,plant.CIN,plant.RD,plant.CD); new=(plant.L,plant.CIN,1.0,470e-6)
for nm,o in [OPS[8],OPS[7],OPS[4]]:
    print(nm,'nominal overshoot %.0f%%, settle(5%%) %.2f ms'%step(o,nom),' | Rd=1 Cd=470u: %.0f%%, %.2f ms'%step(o,new))
