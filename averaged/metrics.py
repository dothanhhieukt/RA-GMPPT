import numpy as np, sim, functools
from plant import battery

def movmean(x, n):
    if n <= 1: return x.copy()
    c = np.cumsum(np.insert(x, 0, 0.0)); y = np.empty_like(x)
    h1 = (n-1)//2; h2 = n - 1 - h1
    N = len(x)
    lo = np.clip(np.arange(N)-h2, 0, N); hi = np.clip(np.arange(N)+h1+1, 0, N)
    return (c[hi]-c[lo])/(hi-lo)

@functools.lru_cache(maxsize=None)
def sweep(pack, case, soc):
    d = sim.run(-2, 0, 1, pack, case, soc0=soc, Tend=0.75, noise=(0, 0), dec=1)
    t, V, P = d['t'], d['V'], d['P']
    dt = t[1]-t[0]; n1 = max(1, round(1e-3/dt))
    V1 = movmean(V, n1); P1 = movmean(P, n1)
    ir = t >= 0.035; vv = V1[ir]; pp = P1[ir]
    Vg = np.arange(np.floor(vv.min()*4)/4, np.ceil(vv.max()*4)/4+0.25, 0.25)
    bi = np.clip(np.round((vv-Vg[0])/0.25).astype(int), 0, len(Vg)-1)
    Pg = np.full(len(Vg), np.nan)
    for b, p in zip(bi, pp):
        if np.isnan(Pg[b]) or p > Pg[b]: Pg[b] = p
    ib = np.nanargmax(Pg); Vbnd = vv.min(); Pbnd = np.median(pp[vv <= Vbnd+0.2])
    return dict(Preach=Pg[ib], Vreach=Vg[ib], Vbnd=Vbnd, Pbnd=Pbnd, Vg=Vg, Pg=Pg)

def metrics(d, pack, caseB, socRef, Ilim=10.0, tsw=np.inf, group='static', Tend=None):
    rf = sweep(pack, caseB, round(socRef, 2))
    t = d['t']; dt = t[1]-t[0]; Tend = Tend or t[-1]+dt
    tStart = 0 if np.isinf(tsw) else tsw
    V, P, Vr, Vb, Ib, D, ph, md = d['V'], d['P'], d['Vref'], d['Vb'], d['Ib'], d['D'], d['phase'], d['mode']
    n1 = max(1, round(1e-3/dt)); n5 = max(1, round(50e-6/dt))
    P1 = movmean(P, n1); Vb1 = movmean(Vb, n1); Ib1 = movmean(Ib, n1); Ib50 = movmean(Ib, n5)
    iw = t >= Tend - 0.05; ie = t >= tStart
    M = dict(Preach=rf['Preach'])
    M['Vss'] = np.median(V[iw]); M['Pss'] = np.median(P1[iw]); M['Ibss'] = Ib[iw].mean(); Vbss = Vb[iw].mean()
    M['isCC'] = rf['Preach']*0.96 > 1.01*Vbss*Ilim
    M['Ebat'] = np.trapezoid((Vb*Ib)[ie], t[ie]); M['Ipk'] = Ib1[ie].max()
    Vcv = battery(pack, 0.5)['Vcv']
    M['VbErr'] = Vbss - Vcv
    # dVcv during CV phase
    ok = ie & (np.abs(Ib1-Ilim) <= 0.02*Ilim); w = max(1, round(10e-3/dt))
    cs = np.cumsum(np.insert(ok.astype(int), 0, 0)); kk = np.where(cs[w:]-cs[:-w] == w)[0]
    M['tset'] = (t[kk[0]]-tStart) if len(kk) else np.nan
    oi = np.where(ie & ~ok)[0]   # settled: end of the last excursion outside +-2 % (1 ms avg)
    if not ok[-1]: M['tset2'] = np.nan
    elif len(oi) == 0: M['tset2'] = 0.0
    else: M['tset2'] = t[oi[-1]+1] - tStart
    kCV = None
    if len(kk):
        c = np.where((np.arange(len(t)) >= kk[0]) & (Vb1 >= Vcv-0.02))[0]
        if len(c): kCV = c[0]
    M['dVcv'] = (Vb1[ie].max() if kCV is None else Vb1[kCV:].max()) - Vcv
    if group == 'cv':
        M['success'] = abs(M['VbErr']) <= 0.2 and M['dVcv'] <= 0.5; M['IbErr'] = M['Ibss']-Ilim
    elif M['isCC']:
        M['IbErr'] = M['Ibss']-Ilim; M['success'] = abs(M['IbErr']) <= 0.02*Ilim
    else:
        M['eff'] = M['Pss']/rf['Preach']
        M['etaE'] = np.trapezoid(P[ie], t[ie])/(rf['Preach']*(t[-1]-tStart))
        M['success'] = M['eff'] >= 0.99
    inB = np.abs(P1-M['Pss']) <= 0.02*max(M['Pss'], 1)
    w = max(1, round(0.02/dt)); cs = np.cumsum(np.insert(inB.astype(int), 0, 0))
    k0 = np.argmax(ie); good = np.where(cs[w:]-cs[:-w] == w)[0]; good = good[good >= k0]
    M['tconv'] = t[good[0]]-tStart if len(good) else np.nan
    M['tUnreach'] = dt*np.sum(ie & (Vr < Vb1/0.95-1.5))
    M['tDsat'] = dt*np.sum(ie & (D >= 0.949))
    M['tOver'] = dt*np.sum(ie & (Ib50 > 1.05*Ilim))
    M['Qover'] = np.trapezoid(np.maximum(Ib50[ie]-Ilim, 0), t[ie])
    pe = ph[ie]; M['nRescan'] = int(np.sum((pe[1:] == 1) & (pe[:-1] == 6)))
    return M
