import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, json
from scipy.linalg import expm
from plant import PVString, L, CIN, RD, CD, RON, battery, ocv
T=50e-6; KI=1000
def op(pack, G, T_, V, soc=0.6):
    pv=PVString(G,T_); I=float(pv.i(V)); g=float(pv.slope(V,0.2)); r=-1/g if g<0 else 1e6
    b=battery(pack,soc); P=V*I
    Ib=P*0.99/ (ocv(b)+b['R']*P/ocv(b)); Vb=ocv(b)+b['R']*Ib
    D0=min(0.95,(Vb+0.025)/(V-0.5-RON*Ib)); return dict(V=V,I=I,r=r,Ib=Ib,Vb=Vb,D0=D0,Rb=b['R'],P=P)
def disc(o):
    V,r,Ib,D0,Rb=o['V'],o['r'],o['Ib'],o['D0'],o['Rb']
    A=np.array([[-1/(r*CIN)-1/(RD*CIN), 1/(RD*CIN), -D0/CIN],[1/(RD*CD), -1/(RD*CD), 0],[D0/L, 0, -(Rb+RON*D0)/L]])
    B=np.array([[-Ib/CIN],[0],[(V-0.5-RON*Ib)/L]]); C=np.array([[1.0,0,0]])
    n=3; M=np.zeros((2*n+2,2*n+2)); M[:n,:n]=A; M[:n,n:n+1]=B
    # augmented: x' = A x + B u ; w' = x (integral of state) ; u'=0
    M=np.zeros((n+n+1,n+n+1)); M[:n,:n]=A; M[:n,2*n:]=B; M[n:2*n,:n]=np.eye(n)
    E=expm(M*T); Phi=E[:n,:n]; Gam=E[:n,2*n:]; Ix=E[n:2*n,:n]; Iu=E[n:2*n,2*n:]
    Cav=C@Ix/T; Dav=C@Iu/T
    return Phi,Gam,Cav,Dav
def loop(o,f,KP=12):
    Phi,Gam,Cav,Dav=disc(o); z=np.exp(2j*np.pi*f*T)
    P=np.array([(Cav@np.linalg.solve(zk*np.eye(3)-Phi,Gam)+Dav)[0,0] for zk in z])
    Dff=o['Vb']/o['V']; Cz=Dff/o['V']*(KP+KI*T/(z-1))
    return -P*Cz/z
def margins(f,Lg):
    mag=np.abs(Lg); ph=np.degrees(np.unwrap(np.angle(Lg)))
    i=np.where((mag[:-1]>=1)&(mag[1:]<1))[0]; pm=fc=np.nan
    if len(i): k=i[-1]; pm=180+ph[k]; fc=f[k]
    j=np.where((ph[:-1]>-180)&(ph[1:]<=-180))[0]; gm=f180=np.nan
    if len(j): k=j[0]; gm=-20*np.log10(mag[k]); f180=f[k]
    return fc,pm,f180,gm
f=np.logspace(0,np.log10(9999),4000)
C1=[1000]*4; C4=[1000,850,500,300]
OPS=[('72 V, Case 1, boundary (83 V)',op(72,C1,25,83.0)),('72 V, Case 1, left flank (96 V)',op(72,C1,25,96)),
     ('72 V, Case 1, MPP (119 V)',op(72,C1,25,118.8)),('72 V, Case 1, right flank (131 V)',op(72,C1,25,131.2)),
     ('48 V, Case 1, boundary (54 V)',op(48,C1,25,54.3)),('48 V, Case 1, left flank (62 V)',op(48,C1,25,62.4)),
     ('48 V, Case 1, right flank (139 V)',op(48,C1,25,139.4)),
     ('72 V, Case 4, peak (96 V)',op(72,C4,25,95.8)),('72 V, Case 4, valley/boundary (82 V)',op(72,C4,25,82.2)),
     ('48 V, Case 4, GMPP (59 V)',op(48,C4,25,59.3))]
rows=[]
for nm,o in OPS:
    fc,pm,f180,gm=margins(f,loop(o,f)); rows.append(dict(op=nm,r=o['r'],D0=o['D0'],fc=fc,PM=pm,f180=f180,GM=gm))
    print('%-38s r=%8.2f D0=%.3f fc=%6.0f Hz PM=%5.1f f180=%6.0f GM=%5.1f dB'%(nm,o['r'],o['D0'],fc,pm,f180,gm))
json.dump(rows,open('bode.json','w'),indent=1)
# GM check vs time-domain: 48 V C1 left flank, KP 24 stable / 32 unstable
o=OPS[5][1]
for KP in [24,28,32]: print('KP',KP,'GM %.2f dB'%margins(f,loop(o,f,KP))[3])
