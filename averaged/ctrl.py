# Python port of the controller in matlab/Charger_Bench.m (main + po_step + dpso_core + cs_core + scan_core).
# Faithful line-by-line translation; alg 2/3/4 (legacy cores, not used in the paper) not included.
import numpy as np
from scipy.special import gamma as _gamma

VB = 200.0; DMAXB = 0.95; HS = 4

class Rng:
    def __init__(self, seed): self.r = np.random.default_rng(int(seed))
    def rand(self): return self.r.random()
    def randn(self): return self.r.standard_normal()

def levy_signed(rng):
    beta = 1.5
    sig = (_gamma(1+beta)*np.sin(np.pi*beta/2)/(_gamma((1+beta)/2)*beta*2**((beta-1)/2)))**(1/beta)
    u = rng.randn()*sig; v = rng.randn()
    return u/abs(v)**(1/beta)

class DPSO:
    def __init__(s, TS):
        s.TS = TS; s.init = False
    def step(s, Vpv, Ipv, useRA, xmaxDyn, restart, useR2, useHT, xWarm, opt):
        NP=3; W=0.4; VMAX=0.075; VSTOP=0.01; DSPREAD=0.45*37.2/200; IT_MAX=20
        XLO=0.33; XHI0=0.90; TSC=s.TS
        SETTLE=round(650e-6*HS/TSC); LONG_SETTLE=round(1300e-6*HS/TSC)
        WARMUP=round(20e-3/TSC); NAVG=max(1,round(200e-6/TSC))
        HC_S0=0.01; HC_SMIN=0.0025; NWIN=20; PTHR=0.10
        DXR=opt.get('dVR',2.0)/200
        useAR = opt.get('AR',False)
        xhi = XHI0
        if useRA or useAR: xhi = max(XLO+0.05, min(XHI0, xmaxDyn))
        xr = max(XLO+0.05, min(XHI0, opt.get('xR2', xmaxDyn)))
        if not s.init:
            s.init=True; s.st=0; s.cnt=0; s.hold=0; s.xc=0.4; s.warmed=False
            s.X=np.zeros(NP); s.V=np.zeros(NP); s.Pb=-np.ones(NP); s.PbX=np.zeros(NP); s.Gb=-1; s.GbX=0.4
            s.k=0; s.it=0; s.Pm=np.zeros(NP); s.Vm=np.zeros(NP); s.Im=np.zeros(NP)
            s.pAcc=0; s.vAcc=0; s.iAcc=0; s.nAcc=0
            s.hcS=HC_S0; s.hcDir=1; s.hcP=0; s.Pw=-1; s.nw=0; s.rst=0
            s.xhiS=xr; s.atB=False; s.lostI=False
        if restart:
            s.cnt=0; s.hold=min(s.hold,SETTLE); s.pAcc=0; s.vAcc=0; s.iAcc=0; s.nAcc=0
        if s.cnt < s.hold:
            s.cnt += 1
            if s.cnt > s.hold-NAVG:
                s.pAcc+=Vpv*Ipv; s.vAcc+=Vpv; s.iAcc+=Ipv; s.nAcc+=1
            if useRA and s.st==3:
                if s.atB and useHT: s.xc = xhi
                else:
                    if (not s.atB) and s.xc > xhi+1e-6: s.lostI=True
                    s.xc = min(s.xc, xhi)
            elif useAR and s.st==3:
                s.xc = min(s.xc, xhi)
            return s.xc, phase_of(s.st), s.rst, s.xc
        P=Vpv*Ipv; Vme=Vpv; Ime=Ipv
        if s.nAcc>0: P=s.pAcc/s.nAcc; Vme=s.vAcc/s.nAcc; Ime=s.iAcc/s.nAcc
        s.pAcc=0; s.vAcc=0; s.iAcc=0; s.nAcc=0
        xn=s.xc; goGlobal=False
        if s.st==0:
            if not s.warmed:
                s.warmed=True; s.xc=xWarm; s.cnt=1; s.hold=WARMUP
                return s.xc, 0, s.rst, s.xc
            goGlobal=True
        elif s.st==1:
            k=s.k
            s.Pm[k]=P; s.Vm[k]=Vme; s.Im[k]=Ime; s.Pb[k]=P; s.PbX[k]=s.X[k]
            if P>s.Gb: s.Gb=P; s.GbX=s.X[k]
            if k<NP-1:
                s.k+=1; xn=s.X[s.k]
            else:
                shaded=(s.Im[2]-s.Im[1])/max(s.Im[2],1e-3)>=0.1 and (s.Vm[0]-s.Vm[1])/max(s.Vm[0],1)>=0.2
                if shaded:
                    s.st=2; s.it=1; s.k=0
                    s.X[0],s.V[0]=dpso_step(s.X[0],s.V[0],s.PbX[0],s.GbX,s.Gb,s.it,XLO,xhi,W,VMAX); xn=s.X[0]
                else:
                    s.st=3; xn=s.GbX; s.hcS=HC_S0; s.hcDir=1; s.hcP=s.Gb; s.Pw=(-1 if opt.get('RSkip') else s.Gb); s.nw=0; s.xhiS=xr
        elif s.st==2:
            k=s.k
            if P>s.Pb[k]: s.Pb[k]=P; s.PbX[k]=s.X[k]
            if P>s.Gb: s.Gb=P; s.GbX=s.X[k]
            if k<NP-1:
                s.k+=1; s.X[s.k],s.V[s.k]=dpso_step(s.X[s.k],s.V[s.k],s.PbX[s.k],s.GbX,s.Gb,s.it,XLO,xhi,W,VMAX); xn=s.X[s.k]
            else:
                s.it+=1
                doStop = s.it>=IT_MAX
                for i in range(NP):
                    for j in range(NP):
                        if i!=j and abs(s.X[i]-s.X[j])<DSPREAD and abs(s.V[i])<VSTOP: doStop=True
                if doStop:
                    s.st=3; xn=s.GbX; s.hcS=HC_S0; s.hcDir=1; s.hcP=s.Gb; s.Pw=(-1 if opt.get('RSkip') else s.Gb); s.nw=0; s.xhiS=xr
                else:
                    s.k=0; s.X[0],s.V[0]=dpso_step(s.X[0],s.V[0],s.PbX[0],s.GbX,s.Gb,s.it,XLO,xhi,W,VMAX); xn=s.X[0]
        else:
            if useRA and useR2 and s.atB and abs(xr-s.xhiS)>DXR:
                s.rst+=1; goGlobal=True; s.r2=getattr(s,'r2',0)+1
            elif useRA and useHT and s.lostI:
                s.rst+=1; goGlobal=True
            s.lostI=False
            s.nw+=1
            if s.nw>=NWIN:
                if s.Pw>0 and abs(P-s.Pw)/s.Pw>PTHR:
                    s.rst+=1; goGlobal=True
                s.Pw=P; s.nw=0
            if not goGlobal:
                if P<s.hcP: s.hcDir=-s.hcDir; s.hcS=max(HC_SMIN,s.hcS/2)
                s.hcP=P
                xn=min(xhi,max(XLO,s.xc+s.hcDir*s.hcS))
        if goGlobal:
            s.st=1; s.k=0; s.Gb=-1; s.it=0; s.xhiS=xr
            for i in range(NP):
                if useRA: s.X[i]=XLO+(i+1)/NP*(xhi-XLO)
                else: s.X[i]=XLO+(i+0.5)/NP*(xhi-XLO)
                s.V[i]=0; s.Pb[i]=-1
            xn=s.X[0]
        if useRA or useAR: xn=min(xn,xhi)
        s.atB = useRA and s.st==3 and xn>=xhi-1e-6
        s.hold=SETTLE
        if abs(xn-s.xc)>0.05: s.hold=LONG_SETTLE
        s.xc=xn; s.cnt=1
        return s.xc, phase_of(s.st), s.rst, s.xc

def dpso_step(Xi,Vi,PbXi,GbX,Gb,it,XLO,xhi,W,VMAX):
    vi=W*Vi+(PbXi-Xi)
    if Gb>=0: vi=vi+(GbX-Xi)
    vi=min(VMAX,max(-VMAX,vi))
    if abs(vi)<1e-4: vi=0.01*(2*(it%2)-1)
    return min(xhi,max(XLO,Xi+vi)), vi

def phase_of(st):
    return {1:1,2:5,3:6}.get(st,0)

def ra_draw(rng,XLO,XHI0,xhi,occB):
    xs=min(xhi,XLO+rng.rand()*(XHI0-XLO))
    if xs>=xhi-1e-9 and xhi<XHI0-1e-9:
        if occB: xs=XLO+rng.rand()*(xhi-XLO)
    return xs

class CS:
    def __init__(s, TS, seed):
        s.TS=TS; s.init=False; s.seed=seed
    def step(s, Vpv, Ipv, useRA, xmaxDyn, restart, useR2, useHT, xWarm, opt):
        NN=3; KAPPA=0.05; PA=0.25; DCONV=0.01; GMAX=25
        XLO=0.33; XHI0=0.90; TSC=s.TS
        SETTLE=round(650e-6*HS/TSC); LONG_SETTLE=round(1300e-6*HS/TSC)
        WARMUP=round(20e-3/TSC); NAVG=max(1,round(200e-6/TSC))
        NWIN=20; PTHR=0.10; DXR=opt.get('dVR',2.0)/200
        useAR=opt.get('AR',False)
        xhi=XHI0
        if useRA or useAR: xhi=max(XLO+0.05,min(XHI0,xmaxDyn))
        xr=max(XLO+0.05,min(XHI0,opt.get('xR2',xmaxDyn)))
        XHIdraw = xhi if useAR else XHI0      # AR: sample only inside reachable range
        if not s.init:
            s.init=True; s.rng=Rng(s.seed*1000+7)
            s.st=0; s.cnt=0; s.hold=0; s.xc=0.4; s.warmed=False
            s.X=np.zeros(NN); s.P=-np.ones(NN); s.Xn=np.zeros(NN); s.k=0; s.g=0; s.Gb=-1; s.GbX=0.4
            s.pAcc=0; s.nAcc=0; s.Pw=-1; s.nw=0; s.rst=0; s.xhiS=xr; s.atB=False; s.lostI=False
            s.mode=1; s.kA=0; s.dropI=np.zeros(NN,bool); s.bndM=False
        rng=s.rng
        if restart: s.cnt=0; s.hold=min(s.hold,SETTLE); s.pAcc=0; s.nAcc=0
        if s.cnt<s.hold:
            s.cnt+=1
            if s.cnt>s.hold-NAVG: s.pAcc+=Vpv*Ipv; s.nAcc+=1
            if useRA and s.st==3:
                if s.atB and useHT: s.xc=xhi
                else:
                    if (not s.atB) and s.xc>xhi+1e-6: s.lostI=True
                    s.xc=min(s.xc,xhi)
            elif useAR and s.st==3:
                s.xc=min(s.xc,xhi)
            return s.xc, phase_of(s.st), s.rst, s.xc
        Pm=Vpv*Ipv
        if s.nAcc>0: Pm=s.pAcc/s.nAcc
        s.pAcc=0; s.nAcc=0
        xn=s.xc; goGlobal=False
        if s.st==0:
            if not s.warmed:
                s.warmed=True; s.xc=xWarm; s.cnt=1; s.hold=WARMUP
                return s.xc,0,s.rst,s.xc
            goGlobal=True
        elif s.st==1:
            if s.bndM:
                s.bndM=False; s.Gb=Pm; s.GbX=s.xc; xn=s.X[0]; s.k=0
                s.hold=SETTLE
                if abs(xn-s.xc)>0.05: s.hold=LONG_SETTLE
                s.xc=xn; s.cnt=1
                return s.xc, phase_of(s.st), s.rst, s.xc
            k=s.k; s.P[k]=Pm
            if Pm>s.Gb: s.Gb=Pm; s.GbX=s.X[k]
            if k<NN-1: s.k+=1; xn=s.X[s.k]
            else:
                s.st=2; s.g=0; s.k=0; s.mode=1; s.Xn=s.X.copy()
                for i in range(NN):
                    s.Xn[i]=min(xhi,max(XLO,s.X[i]+KAPPA*(s.GbX-s.X[i])*levy_signed(rng)))
                xn=s.Xn[0]
        elif s.st==2:
            endGen=False
            if s.mode==1:
                k=s.k
                if Pm>s.P[k]: s.P[k]=Pm; s.X[k]=s.Xn[k]
                if Pm>s.Gb: s.Gb=Pm; s.GbX=s.Xn[k]
                s.k+=1
                if s.k<=NN-1: xn=s.Xn[s.k]
                elif (s.X.max()-s.X.min()<DCONV) or s.g+1>=GMAX:
                    s.st=3; xn=s.GbX; s.Pw=s.Gb; s.nw=0; s.xhiS=xr
                else:
                    nA=0
                    for i in range(NN):
                        s.dropI[i] = abs(s.X[i]-s.GbX)>1e-9 and rng.rand()<PA
                        if s.dropI[i]:
                            nA+=1; occB=False
                            for j in range(NN):
                                if j!=i:
                                    xj=s.X[j]
                                    if j<i and s.dropI[j]: xj=s.Xn[j]
                                    if xj>=xhi-1e-9: occB=True
                            s.Xn[i]=ra_draw(rng,XLO,XHIdraw,xhi,occB) if not useAR else XLO+rng.rand()*(xhi-XLO)
                    if nA>0:
                        s.kA=int(np.argmax(s.dropI)); s.mode=2; xn=s.Xn[s.kA]
                    else: endGen=True
            else:
                kA=s.kA
                s.X[kA]=s.Xn[kA]; s.P[kA]=Pm; s.dropI[kA]=False
                if Pm>s.Gb: s.Gb=Pm; s.GbX=s.X[kA]
                if s.dropI.any():
                    s.kA=int(np.argmax(s.dropI)); xn=s.Xn[s.kA]
                else: endGen=True
            if endGen:
                s.g+=1
                if s.g>=GMAX:
                    s.st=3; xn=s.GbX; s.Pw=s.Gb; s.nw=0; s.xhiS=xr
                else:
                    s.mode=1; s.k=0
                    for i in range(NN):
                        s.Xn[i]=min(xhi,max(XLO,s.X[i]+KAPPA*(s.GbX-s.X[i])*levy_signed(rng)))
                    xn=s.Xn[0]
        else:
            if useRA and useR2 and s.atB and abs(xr-s.xhiS)>DXR:
                s.rst+=1; goGlobal=True; s.r2=getattr(s,'r2',0)+1
            elif useRA and useHT and s.lostI:
                s.rst+=1; goGlobal=True
            s.lostI=False
            s.nw+=1
            if (not goGlobal) and s.nw>=NWIN:
                if s.Pw>0 and abs(Pm-s.Pw)/s.Pw>PTHR:
                    s.rst+=1; goGlobal=True
                s.Pw=Pm; s.nw=0
            if not goGlobal: xn=s.xc
        if goGlobal:
            s.st=1; s.k=0; s.Gb=-1; s.g=0; s.mode=1; s.xhiS=xr
            for i in range(NN):
                occB=any(s.X[j]>=xhi-1e-9 for j in range(i))
                s.X[i]=ra_draw(rng,XLO,XHI0,xhi,occB) if not useAR else XLO+rng.rand()*(xhi-XLO)
                s.P[i]=-1
            xn=s.X[0]
            if useRA: s.bndM=True; xn=xhi
        if useRA or useAR: xn=min(xn,xhi)
        s.atB = useRA and s.st==3 and xn>=xhi-DXR
        if s.atB: xn=xhi
        s.hold=SETTLE
        if abs(xn-s.xc)>0.05: s.hold=LONG_SETTLE
        s.xc=xn; s.cnt=1
        return s.xc, phase_of(s.st), s.rst, s.xc


SENS = [  # Dhat, Voff, dVc, dVR, tauR2, ocvR2   (index = floor(abl/4096))
    (0.95, 0.0, 0.5, 2.0, 5e-4, 0), (0.93, 0.0, 0.5, 2.0, 5e-4, 0), (0.97, 0.0, 0.5, 2.0, 5e-4, 0),
    (0.95, 1.0, 0.5, 2.0, 5e-4, 0), (0.95, -1.0, 0.5, 2.0, 5e-4, 0), (0.95, 0.0, 0.0, 2.0, 5e-4, 0),
    (0.95, 0.0, 1.5, 2.0, 5e-4, 0), (0.95, 0.0, 0.5, 1.0, 5e-4, 0), (0.95, 0.0, 0.5, 4.0, 5e-4, 0),
    (0.95, 0.0, 0.5, 2.0, 50e-3, 0), (0.95, 0.0, 0.5, 2.0, 5e-4, 1)]

class Controller:
    """alg 1 P&O, 5 DPSO, 6 RA-DPSO, 7 CS, 8 RA-CS, 9 SCAN, 10 RA-SCAN, -2 sweep.
    abl bits: 1 PI right-flank limiter, 2 fast left-flank reaction, 16 no RA, 32 conventional limiter,
    64 no R2, 128 no HT (hold tracking), 256 adaptive range (baselines), 512 right-flank soft start, 1024 20-kHz mode;
    floor(abl/4096) = sensitivity row."""
    def __init__(s, alg, abl=0, seed=1, Ilim=10.0, Vcv=54.6, TS=5e-6, opt=None):
        s.alg=alg; s.abl=int(abl); s.seed=seed; s.Ilim=Ilim; s.Vcv=Vcv; s.TS=TS
        si = s.abl//4096; s.sens = SENS[si]
        s.opt = dict(opt or {}); s.opt.setdefault('dVR', s.sens[3]); s.opt['AR'] = (s.abl & 256) != 0; s.opt['RSkip'] = (s.abl & 8) != 0   # DPSO restart rule: first comparison after convergence skipped
        s.dpso=DPSO(TS); s.cs=CS(TS,seed); s.scan=SCAN(TS)
        s.init=False
    def __call__(s, Vpv, Ipv, Vbat, Ibat):
        TS=s.TS; TAU=5e-4; aF=TS/TAU
        NPWM=max(1,round(50e-6/TS)); KP=s.opt.get('KP',12); KI=1000
        LIM_PER=round(650e-6*HS/TS)
        kL=0.05; kR=0.01; dxMax=0.01; ETA=0.97
        SR=10e3; dxSR=SR*NPWM*TS/VB
        kC=0.02; NEX=5; T_ENT=round(3e-3/TS); T_FAST=max(1,round(0.3e-3/TS))
        KPL=0.025; KIL=50.0                         # PI right-flank limiter (x per unit error; 1/s)
        alg=s.alg; abl=s.abl; opt=s.opt
        Dhat, Voff, dVc, dVR, tauR2, ocvR2 = s.sens
        Vbat = Vbat + Voff                          # battery-voltage sensor offset
        if not s.init:
            s.init=True; s.VbF=Vbat; s.IbF=0; s.IpF=0; s.lim=False; s.xL=0.5; s.xPk=0.5
            s.side=1; s.cnt=0; s.xLast=0.4; s.olCnt=0
            s.poX=0.4; s.poP=0; s.poDir=1; s.poCnt=0; s.poW=False
            s.vAcc=0; s.nAcc=0; s.Dcmd=0; s.integ=0; s.restartCore=False
            s.lastPh=0; s.lastRst=0; s.lastXh=0.4
            s.xTrim=0; s.xCmd=0.4; s.slewing=False; s.xFFl=0.4
            s.Icv=10; s.exCnt=0; s.exhaust=False; s.entCnt=0
            s.fromL=False; s.Ip0=0; s.kT=0; s.ssDone=False
            s.VbS=Vbat; s.fCnt=0; s.piMode=False; s.xI=0.4; s.fInh=False
        s.VbF+=aF*(Vbat-s.VbF); s.IbF+=aF*(Ibat-s.IbF); s.IpF+=aF*(Ipv-s.IpF)
        s.VbS+=(TS/tauR2)*(Vbat-s.VbS)
        VpvMin=max(s.VbF,1)/Dhat+1.0
        VpvMinC=max(s.VbF,1)/Dhat-dVc
        xmaxC=1-VpvMinC/VB; xmax=1-VpvMin/VB
        Rhat = 0.15 if s.Vcv > 70 else 0.0975
        VbR2 = s.VbF - Rhat*s.IbF if ocvR2 else s.VbS        # boundary signal used by R2
        xR2 = 1-(max(VbR2,1)/Dhat-dVc)/VB
        opt['xR2'] = xR2
        propC = alg in (6,8,10)
        rightLim = (not propC) or (abl & 32)!=0
        raOn = (abl & 16)==0
        useR2 = raOn and (abl & 64)==0
        useHT = raOn and (abl & 128)==0
        piLim = rightLim and (abl & 1)!=0
        fast = (not rightLim) and (abl & 2)!=0
        ssR = rightLim and (abl & 512)!=0
        xWarm=0.4
        if (not rightLim) and propC: xWarm=xmaxC
        if ssR: xWarm=1-150/VB                      # right-flank soft start: warm up near Voc
        s.vAcc+=Vpv; s.nAcc+=1
        phaseOut=s.lastPh; modeOut=0
        if alg==-2:
            s.olCnt+=1
            Vsw=150-130*min(1,max(0,(s.olCnt*TS-0.03)/0.7)); x=1-Vsw/VB
        else:
            IlimS=max(s.Ilim,1e-3)
            e=max(s.IbF/IlimS-1,(s.VbF-s.Vcv)/0.5)
            wasLim=s.lim
            if e>0.02: s.entCnt+=1
            else: s.entCnt=0
            if s.lastPh in (1,5): s.fInh=False                         # a new search of the core re-enables fast entry
            if fast and (not s.fInh) and s.lastPh==6 and e>-0.05: s.fCnt+=1     # near the limit while the core holds
            else: s.fCnt=0
            s.kT+=1
            if (not s.ssDone) and s.kT>=round((2e-3 if fast else 19.5e-3)/TS) and (((not rightLim) and propC) or ssR):
                s.ssDone=True
                if not s.lim:
                    s.lim=True; s.cnt=0; s.xTrim=0; s.slewing=False; s.exCnt=0; s.exhaust=False; s.fromL=False
                    s.entCnt=0; s.Icv=IlimS; s.xPk=1-119/VB
                    if ssR:                                       # right flank from ~Voc down to 119 V
                        s.side=-1; s.xL=1-150/VB; s.xI=s.xL; s.piMode=piLim
                    else:
                        s.side=1; s.xL=xmax; s.xCmd=xmax; s.xFFl=xmax
            if (not s.lim) and (s.entCnt>=T_ENT or (fast and s.fCnt>=T_FAST)):
                fastE = s.entCnt < T_ENT
                s.entCnt=0; s.fCnt=0
                s.lim=True; s.cnt=0; s.xTrim=0; s.slewing=False; s.exCnt=0; s.exhaust=False; s.fromL=False
                if alg>=2 and s.lastPh in (3,9,6,4): s.xPk=s.lastXh
                else: s.xPk=s.xLast
                s.xL=s.xLast; s.xCmd=s.xLast; s.xFFl=s.xPk
                s.Icv=min(IlimS,max(s.IbF,0.2))
                if fastE: s.Icv=IlimS; s.xPk=min(s.xPk,1-119/VB)    # held point is no longer the peak after a rise
                if rightLim:
                    if piLim:                                     # min-select PI: no jump, start at held point
                        s.side=-1; s.piMode=True; s.xI=s.xLast; s.xPk=s.xLast
                    else:
                        s.side=-1; s.xL=1-150/VB; s.fromL=True; s.Ip0=max(s.IpF,1)
                else: s.side=1
            if s.lim:
                Iset=max(0.2,min(IlimS,s.Icv)); eI=s.IbF/Iset-1
                s.cnt+=1
                if s.side<0 and s.piMode:
                    if s.nAcc>=NPWM:                             # PI on the constraint error, every PWM period
                        xIn = s.xI - KIL*NPWM*TS*e
                        s.xI = min(s.xPk, max(0.15, xIn))
                        s.xL = min(s.xPk, max(0.15, s.xI - KPL*e))
                        if s.xI>=s.xPk-1e-9 and e<-0.02: s.lim=False; s.piMode=False
                elif s.cnt>=LIM_PER:
                    s.cnt=0
                    if s.side>0:
                        s.exCnt = s.exCnt+1 if s.exhaust else 0
                        s.Icv=min(IlimS,max(0,s.Icv-IlimS*(kC*6 if fast else kC)*(s.VbF-s.Vcv)/0.5))
                        xUn=s.xFFl+s.xTrim
                        if (not s.slewing) and not (xUn>=xmax and eI>0) and not (xUn<=s.xPk and eI<0):
                            s.xTrim=min(1,max(-1,s.xTrim+min(dxMax,max(-dxMax,kL*eI))))
                    else:
                        if s.fromL and s.IpF<0.15*s.Ip0 and e<-0.1: s.xL+=2*dxMax
                        else: s.xL-=min(dxMax,max(-dxMax,kR*e))
                        s.xL=min(xmax,max(s.xL,0.15))
                        if s.xL>=s.xPk and e<-0.02: s.lim=False
                if s.lim and s.side>0 and s.nAcc>=NPWM:
                    if s.IpF>0.5:
                        Vt=s.VbF*Iset/(ETA*s.IpF)
                        xFF=1-min(VB*(1-s.xPk),max(VpvMin,Vt))/VB
                    else: xFF=s.xPk
                    s.xFFl=xFF
                    xDes=min(xmax,max(s.xPk,xFF+s.xTrim))
                    if fast and xDes>s.xCmd:                     # power-reducing direction: no slew limit
                        s.slewing=False; s.xCmd=xDes
                    else:
                        s.slewing=abs(xDes-s.xCmd)>dxSR
                        s.xCmd+=min(dxSR,max(-dxSR,xDes-s.xCmd))
                    s.exhaust=(not s.slewing) and s.xCmd>=xmax-1e-6 and e>0.02
                    cvNoLeft = s.VbF>=s.Vcv-0.1 and s.Icv<IlimS-1e-6 and s.IpF>0.5 and VpvMin*s.IpF*ETA>1.05*s.VbF*Iset
                    if fast and s.xCmd>=xmax-1e-6 and s.IpF>0.5 and Vpv>=VpvMin-1.5 and Vpv*s.IpF*ETA>1.05*s.VbF*Iset:
                        cvNoLeft = True                              # no left-flank solution at all (CC too)
                    if (s.exhaust and s.exCnt>=NEX) or cvNoLeft:
                        s.side=-1; s.xL=1-150/VB; s.fromL=True; s.Ip0=max(s.IpF,1)
                    elif xDes<=s.xPk+1e-6 and s.xCmd<=s.xPk+1e-6 and eI<-0.02 and s.Icv>=IlimS-1e-6:
                        s.lim=False
                    else: s.xL=s.xCmd
                if wasLim and not s.lim:
                    s.restartCore=True; s.poCnt=0; s.poP=-1e9; s.fInh=True
                x=s.xL; phaseOut=20; modeOut = 1 if s.side>0 else 2
            else:
                if alg==1:
                    x=s.po_step(Vpv*Ipv,xmax if not opt['AR'] else xmaxC); ph=0; rs=0; xh=x
                elif alg in (5,6):
                    x,ph,rs,xh=s.dpso.step(Vpv,Ipv,alg==6 and raOn,xmaxC,s.restartCore,useR2,useHT,xWarm,opt)
                elif alg in (7,8):
                    x,ph,rs,xh=s.cs.step(Vpv,Ipv,alg==8 and raOn,xmaxC,s.restartCore,useR2,useHT,xWarm,opt)
                else:
                    x,ph,rs,xh=s.scan.step(Vpv,Ipv,alg==10 and raOn,xmaxC,s.restartCore,useR2,useHT,xWarm,opt)
                s.restartCore=False
                s.xLast=x; s.lastPh=ph; s.lastRst=rs; s.lastXh=xh; phaseOut=ph
        Vref=VB*(1-x)
        if s.nAcc>=NPWM:
            vAvg=s.vAcc/s.nAcc; s.vAcc=0; s.nAcc=0
            Dff=s.VbF/max(Vref,1); ev=(vAvg-Vref)/max(Vref,1)
            Dnew=Dff*(1+KP*ev+s.integ)
            Dsat=min(DMAXB,max(0,Dnew))
            if Dsat==Dnew or (Dnew>DMAXB and ev<0) or (Dnew<0 and ev>0):
                s.integ+=KI*NPWM*TS*ev
            s.Dcmd=Dsat
        return s.Dcmd, phaseOut, modeOut, Vref, VpvMin
    def po_step(s,P,xmax):
        s.poCnt+=1
        if not s.poW:
            if s.poCnt>=round(20e-3/s.TS): s.poW=True; s.poCnt=0; s.poP=P
        else:
            if s.poCnt>=round(650e-6*HS/s.TS):
                s.poCnt=0
                if P<s.poP: s.poDir=-s.poDir
                s.poP=P
                s.poX=min(xmax,max(0.15,s.poX+s.poDir*0.005))
        return s.poX


class SCAN:
    """Third core: fixed-point scan at k*0.8*Voc,module (Patel & Agarwal, IEEE TIE 2008) + P&O.
    RA: points projected by R1 (boundary always a candidate, duplicates merged), P&O steps projected,
    R2 restarts the scan when the P&O point stays within dVR of the boundary and the boundary has moved > dVR."""
    def __init__(s, TS): s.TS = TS; s.init = False
    def step(s, Vpv, Ipv, useRA, xmaxDyn, restart, useR2, useHT, xWarm, opt):
        XLO = 0.33; XHI0 = 0.90; TSC = s.TS
        SETTLE = round(650e-6*HS/TSC); LONG_SETTLE = round(1300e-6*HS/TSC)
        WARMUP = round(20e-3/TSC); NAVG = max(1, round(200e-6/TSC))
        NWIN = 20; PTHR = 0.10; DXR = opt.get('dVR', 2.0)/200; DPO = 0.005
        useAR = opt.get('AR', False)
        xhi = XHI0
        if useRA or useAR: xhi = max(XLO+0.05, min(XHI0, xmaxDyn))
        xr = max(XLO+0.05, min(XHI0, opt.get('xR2', xmaxDyn)))
        if not s.init:
            s.init = True; s.st = 0; s.cnt = 0; s.hold = 0; s.xc = 0.4; s.warmed = False
            s.Xs = np.zeros(5); s.ns = 0; s.k = 0; s.Gb = -1; s.GbX = 0.4
            s.pAcc = 0; s.nAcc = 0; s.Pw = -1; s.nw = 0; s.rst = 0; s.xhiS = xr; s.atB = False
            s.poP = 0; s.poDir = 1
        if restart: s.cnt = 0; s.hold = min(s.hold, SETTLE); s.pAcc = 0; s.nAcc = 0
        if s.cnt < s.hold:
            s.cnt += 1
            if s.cnt > s.hold-NAVG: s.pAcc += Vpv*Ipv; s.nAcc += 1
            if (useRA or useAR) and s.st == 3:
                s.xc = min(s.xc, xhi)
            return s.xc, phase_of(s.st), s.rst, s.xc
        Pm = Vpv*Ipv
        if s.nAcc > 0: Pm = s.pAcc/s.nAcc
        s.pAcc = 0; s.nAcc = 0
        xn = s.xc; goGlobal = False
        if s.st == 0:
            if not s.warmed:
                s.warmed = True; s.xc = xWarm; s.cnt = 1; s.hold = WARMUP
                return s.xc, 0, s.rst, s.xc
            goGlobal = True
        elif s.st == 1:                                   # scan points
            if Pm > s.Gb: s.Gb = Pm; s.GbX = s.Xs[s.k]
            if s.k < s.ns-1:
                s.k += 1; xn = s.Xs[s.k]
            else:
                s.st = 3; xn = s.GbX; s.poP = s.Gb; s.poDir = 1; s.Pw = s.Gb; s.nw = 0; s.xhiS = xr
        else:                                             # P&O around the best point
            if useRA and useR2 and s.atB and abs(xr-s.xhiS) > DXR:
                s.rst += 1; goGlobal = True; s.r2 = getattr(s,'r2',0)+1
            s.nw += 1
            if (not goGlobal) and s.nw >= NWIN:
                if s.Pw > 0 and abs(Pm-s.Pw)/s.Pw > PTHR:
                    s.rst += 1; goGlobal = True
                s.Pw = Pm; s.nw = 0
            if not goGlobal:
                if Pm < s.poP: s.poDir = -s.poDir
                s.poP = Pm
                xn = min(xhi, max(XLO, s.xc + s.poDir*DPO))
        if goGlobal:
            s.st = 1; s.k = 0; s.Gb = -1; s.xhiS = xr
            pts = [1 - kk*0.8*37.2/VB for kk in (4, 3, 2, 1)]          # 119.0, 89.3, 59.5, 29.8 V
            if useRA:
                pts = [p for p in pts if p < xhi - DXR] + [xhi]       # R1: projected, boundary once
            elif useAR:
                pts = [p for p in pts if p <= xhi]
                if not pts: pts = [xhi]
            s.ns = len(pts); s.Xs[:s.ns] = pts; xn = s.Xs[0]
        if useRA or useAR: xn = min(xn, xhi)
        nowB = useRA and s.st == 3 and xn >= xhi - DXR
        if nowB and not s.atB: s.xhiS = xr                       # R2 reference = start of the boundary hold
        s.atB = nowB
        s.hold = SETTLE
        if abs(xn-s.xc) > 0.05: s.hold = LONG_SETTLE
        s.xc = xn; s.cnt = 1
        return s.xc, phase_of(s.st), s.rst, s.xc
