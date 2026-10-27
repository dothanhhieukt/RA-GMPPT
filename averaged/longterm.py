# One-hour charging run on the averaged model: EN 50530-like irradiance ramps, moving partial shading,
# irradiance-dependent module temperature, real battery capacity. Controller at 20 kHz (TS = 50 us).
import numpy as np, sys, json, math, time
from plant import PVString, FastPV, battery, SOCtab, OCVtab, L, CIN, RD, CD, VF_DB, VF_DF, RON
from ctrl import Controller

TEND = 3600.0
def profile():
    """returns arrays (t, Gsun) at 1 s and the shading pattern schedule"""
    # EN 50530-style ramps between 30 % and 100 % irradiance with slopes 10..100 W/m2/s and holds
    seg = []
    G = 300.0
    for slope, hold in [(10, 60), (20, 60), (50, 30), (100, 30), (10, 120), (30, 60), (5, 120), (100, 20), (50, 60), (20, 60)]:
        for target in (1000.0, 300.0):
            n = int(abs(target-G)/slope)
            seg += list(np.linspace(G, target, n+1)[1:]); G = target
            seg += [G]*hold
    g = np.array(seg)
    while len(g) < TEND+1: g = np.concatenate([g, g])
    return g[:int(TEND)+1]
# shading factors per module (moving shadow), changed every 300 s
PATTERNS = [[1,1,1,1], [1,1,.6,.6], [1,.85,.5,.3], [1,.7,.7,.3], [.9,.7,.4,.2], [1,.85,.5,.3], [1,1,1,.4],
            [1,.6,.6,.4], [1,.85,.5,.3], [1,.9,.5,.5], [1,1,.6,.6], [1,1,1,1]]

_cache = {}
def pv_at(Gs, T):
    key = tuple(int(round(g/10))*10 for g in Gs) + (round(T*2)/2,)
    if key not in _cache:
        if len(_cache) > 400: _cache.clear()
        _cache[key] = FastPV(PVString(list(key[:4]), key[4], n=3000))
    return _cache[key]

def run(alg, abl, seed=1, pack=72, soc0=0.30, Ilim=10.0, dt=25e-6, cdec=2, noise=(0.05, 0.005), tend=TEND, log=None):
    b = battery(pack, soc0); Ns = b['Ns']; Rb = b['R']; Q = b['Q']; soc = soc0
    TS = dt*cdec
    ctrl = Controller(alg, abl, seed, Ilim, b['Vcv'], TS)
    gsun = profile()
    v = vd = 0.0; iL = 0.0; D = 0.0; Vref = 0.0; ph = 0
    nz = np.random.default_rng(1000+seed)
    oc = Ns*float(np.interp(soc, SOCtab, OCVtab))
    E = 0.0; Epv = 0.0; tUnr = 0.0; tB = 0.0; nOS = 0
    per = []        # per-second records: t, Ppv, Pbat, Vbat, Vpv, Vref
    nsec = int(round(1.0/dt)); N = int(round(tend/dt))
    accP = accPb = accVb = accV = 0.0
    lastPh = 0; nScan = 0
    for n in range(N):
        if n % nsec == 0:
            sec = n//nsec
            Gs = [gsun[sec]*f for f in PATTERNS[(sec//300) % len(PATTERNS)]]
            T = 25 + 0.025*float(np.mean(Gs))
            pv = pv_at(Gs, T)
            NZ = nz.standard_normal((2, nsec//cdec+2)); jz = 0
            oc = Ns*float(np.interp(soc, SOCtab, OCVtab))
            if n > 0:
                per.append((sec, accP/nsec, accPb/nsec, accVb/nsec, accV/nsec, Gs, T))
                accP = accPb = accVb = accV = 0.0
        ipv, g = pv.i_g(v)
        if n % cdec == cdec-1:
            D, ph, mode, Vref, Vmin = ctrl(v + noise[0]*NZ[0][jz], ipv + noise[1]*NZ[1][jz], oc + Rb*iL, iL); jz += 1
            if ph == 1 and lastPh == 6: nScan += 1
            lastPh = ph
            vbP = oc + Rb*iL
            if Vref < vbP/0.95 - 1.5: tUnr += TS
            if abs(v - Vmin + 1.0) < 1.0 and ph == 6: tB += TS     # operating at the boundary
        f = ipv - D*iL - (v - vd)/RD
        v = v + dt*f/(CIN - dt*(g - 1/RD))
        vd += dt*(v - vd)/(RD*CD)
        vL = D*(v - VF_DB - RON*iL) - (1-D)*VF_DF - (oc + Rb*iL)
        iL = iL + dt*vL/L
        if iL < 0: iL = 0.0
        soc += dt*iL/(3600*Q)
        pb = (oc + Rb*iL)*iL; pp = v*ipv
        E += pb*dt; Epv += pp*dt
        accP += pp; accPb += pb; accVb += oc + Rb*iL; accV += v
    r2 = getattr(ctrl.dpso, 'r2', 0) + getattr(ctrl.cs, 'r2', 0) + getattr(ctrl.scan, 'r2', 0)
    return dict(alg=alg, abl=abl, seed=seed, E=E, Epv=Epv, socEnd=soc, tUnr=tUnr, tBnd=tB, nScan=nScan, nR2=r2), per

def reach_energy(per, pack):
    """energy of the reachable maximum along the run (static curve, boundary from the simulated V_bat)"""
    Er = 0.0; Eg = 0.0
    for sec, P, Pb, Vb, Vpv, Gs, T in per:
        pv = pv_at(Gs, T).pvs
        V = np.linspace(1, pv.V.max(), 3000); Pc = V*pv.i(V)
        Ib = Pb/max(Vb, 1); vb0 = (Vb + 0.025)/0.95 + 0.5 + 0.02*Ib
        m = V >= vb0
        Pr = max(Pc[m].max() if m.any() else 0, float(np.interp(vb0, V, Pc)))
        Ilim_P = 10*Vb/0.99
        Er += min(Pr, Ilim_P); Eg += min(Pc.max(), Ilim_P)
    return Er, Eg

if __name__ == '__main__':
    alg, abl = int(sys.argv[1]), int(sys.argv[2]); tend = float(sys.argv[3]) if len(sys.argv) > 3 else TEND
    t0 = time.time()
    R, per = run(alg, abl, tend=tend)
    Er, Eg = reach_energy(per, 72)
    R.update(Ereach=Er, Egmpp=Eg, wall=time.time()-t0, tend=tend)
    json.dump(R, open(f'lt_{alg}_{abl}_{int(tend)}.json', 'w'))
    np.save(f'lt_{alg}_{abl}_{int(tend)}_per.npy', np.array([(s, P, Pb, Vb, Vpv, T) for s, P, Pb, Vb, Vpv, Gs, T in per]))
    print(R)
