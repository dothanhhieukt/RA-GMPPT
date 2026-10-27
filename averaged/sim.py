import numpy as np, math
from plant import PVString, FastPV, battery, ocv, L, CIN, RD, CD, VF_DB, VF_DF, RON
from ctrl import Controller
CASES = {'C1':([1000]*4,25),'C2':([1000,1000,600,600],25),'C3':([1000,700,700,300],25),'C4':([1000,850,500,300],25),
         'C5':([1000,600,600,400],25),'C6':([900,700,400,200],25),'C7':([1000,850,500,300],50),'C8':([1000,850,500,300],65)}
_pvcache = {}
def pvfor(case):
    if case not in _pvcache:
        G, T = CASES[case]; _pvcache[case] = FastPV(PVString(G, T))
    return _pvcache[case]

def run(alg, abl=0, seed=1, pack=72, caseA='C4', caseB=None, soc0=0.6, Qscale=1.0, tsw=np.inf, Tend=0.8,
        Ilim=10.0, opt=None, TS=5e-6, dec=10, nsub=2, Gprof=None, noise=(0.05,0.005), cdec=None):
    caseB = caseB or caseA
    b = battery(pack, soc0, Qscale)
    if cdec is None: cdec = 10 if (int(abl) & 1024) else 1
    ctrl = Controller(alg, abl, seed, Ilim, b['Vcv'], TS*cdec, opt)
    D = 0.0; ph = 0; mode = 0; Vref = 0.0
    pvA, pvB = pvfor(caseA), pvfor(caseB)
    Ns = b['Ns']; Rb = b['R']; Q = b['Q']; soc = soc0
    v = 0.0; vd = 0.0; iL = 0.0
    nz = np.random.default_rng(1000+int(seed)); sv, si = noise
    NZ = nz.standard_normal((2, 4096)).tolist()
    N = int(round(Tend/TS)); h = TS/nsub
    out = {k: [] for k in ('t','V','P','Vref','Vb','Ib','D','phase','mode','SOC')}
    from plant import SOCtab, OCVtab
    def OCV(s): return Ns*float(np.interp(s, SOCtab, OCVtab))
    oc = OCV(soc)
    for n in range(N):
        t = n*TS
        pv = pvA if t < tsw else pvB
        ipv, _ = pv.i_g(v)
        vb = oc + Rb*iL
        j = n & 4095
        if n % cdec == cdec-1:
            D, ph, mode, Vref, Vmin = ctrl(v + sv*NZ[0][j], ipv + si*NZ[1][j], vb, iL)
        if j == 4095: NZ = nz.standard_normal((2, 4096)).tolist()
        for _ in range(nsub):
            ipv, g = pv.i_g(v)
            f = ipv - D*iL - (v - vd)/RD
            v = v + h*f/(CIN - h*(g - 1/RD))
            vd += h*(v - vd)/(RD*CD)
            vL = D*(v - VF_DB - RON*iL) - (1-D)*VF_DF - (oc + Rb*iL)
            iL = iL + h*vL/L
            if iL < 0: iL = 0.0
            soc += h*iL/(3600*Q)
        if n % 200 == 0: oc = OCV(soc)
        if n % dec == 0:
            out['t'].append(t); out['V'].append(v); out['P'].append(v*ipv); out['Vref'].append(Vref)
            out['Vb'].append(oc + Rb*iL); out['Ib'].append(iL); out['D'].append(D); out['phase'].append(ph)
            out['mode'].append(mode); out['SOC'].append(soc)
    return {k: np.array(x) for k, x in out.items()}
