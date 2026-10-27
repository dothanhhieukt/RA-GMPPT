import numpy as np, itertools, json
from scipy.linalg import expm
import plant
from plant import PVString, RON, battery, ocv
T=50e-6; KI=1000
exec(open('bode.py').read().split("def disc(o):")[0].split("T=50e-6; KI=1000")[1])   # op()
def disc(o,L,CIN,RD,CD):
    V,r,Ib,D0,Rb=o['V'],o['r'],o['Ib'],o['D0'],o['Rb']
    A=np.array([[-1/(r*CIN)-1/(RD*CIN), 1/(RD*CIN), -D0/CIN],[1/(RD*CD), -1/(RD*CD), 0],[D0/L, 0, -(Rb+RON*D0)/L]])
    B=np.array([[-Ib/CIN],[0],[(V-0.5-RON*Ib)/L]]); C=np.array([[1.0,0,0]]); n=3
    M=np.zeros((2*n+1,2*n+1)); M[:n,:n]=A; M[:n,2*n:]=B; M[n:2*n,:n]=np.eye(n)
    E=expm(M*T); Phi=E[:n,:n]; Gam=E[:n,2*n:]; Ix=E[n:2*n,:n]; Iu=E[n:2*n,2*n:]
    return Phi,Gam,C@Ix/T,C@Iu/T
def loop(o,f,p,KP=12):
    Phi,Gam,Cav,Dav=disc(o,*p); z=np.exp(2j*np.pi*f*T)
    P=np.array([(Cav@np.linalg.solve(zk*np.eye(3)-Phi,Gam)+Dav)[0,0] for zk in z])
    return -P*(o['Vb']/o['V'])/o['V']*(KP+KI*T/(z-1))/z
def margins(f,Lg):
    mag=np.abs(Lg); ph=np.degrees(np.unwrap(np.angle(Lg)))
    i=np.where((mag[:-1]>=1)&(mag[1:]<1))[0]; pm=np.nan
    if len(i): pm=180+ph[i[-1]]
    j=np.where((ph[:-1]>-180)&(ph[1:]<=-180))[0]; gm=np.nan
    if len(j): gm=-20*np.log10(mag[j[0]])
    return pm,gm
f=np.logspace(0,np.log10(9999),1500)
C1=[1000]*4; C4=[1000,850,500,300]
OPS=[('72 C1 83',op(72,C1,25,83.0)),('72 C1 96',op(72,C1,25,96)),('72 C1 119',op(72,C1,25,118.8)),('72 C1 131',op(72,C1,25,131.2)),
     ('48 C1 54',op(48,C1,25,54.3)),('48 C1 62',op(48,C1,25,62.4)),('48 C1 139',op(48,C1,25,139.4)),
     ('72 C4 96',op(72,C4,25,95.8)),('72 C4 82',op(72,C4,25,82.2)),('48 C4 59',op(48,C4,25,59.3))]
nom=(plant.L,plant.CIN,plant.RD,plant.CD); print('nominal',nom)
res=[]
for nm,o in OPS:
    pms=[];gms=[]
    for a,b,c,d in itertools.product([0.8,1,1.2],repeat=4):
        p=(nom[0]*a,nom[1]*b,nom[2]*c,nom[3]*d); pm,gm=margins(f,loop(o,f,p)); pms.append(pm); gms.append(gm)
    pn,gn=margins(f,loop(o,f,nom))
    res.append((nm,pn,np.nanmin(pms),gn,np.nanmin(gms)))
    print('%-10s PM nom %5.1f min %5.1f | GM nom %5.1f min %5.1f'%res[-1])
json.dump(res,open('bode_tol.json','w'))
