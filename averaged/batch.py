import numpy as np, sim, metrics, multiprocessing as mp, sys, json, itertools
def job(a):
    name, alg, abl, seed, pack, cA, cB, soc0, Q, tsw, Tend, socRef, Ilim, group, kw = a
    d = sim.run(alg, abl, seed, pack, cA, cB, soc0, Q, tsw, Tend, Ilim=Ilim, **kw)
    M = metrics.metrics(d, pack, cB, socRef, Ilim, tsw, group)
    M.update(name=name, alg=alg, abl=abl, seed=seed, scen=f'{pack}-{cA}{">"+cB if cB!=cA else ""}-S{round(100*soc0):02d}' + (f'-I{Ilim:g}' if Ilim != 10 else '') + (f'-Q{Q:g}' if Q < 1 else ''))
    return M
def runall(jobs, nproc=2):
    with mp.Pool(nproc) as p: return p.map(job, jobs, chunksize=1)
def J(name, alg, abl, seed, pack, case, soc, Tend=0.8, Ilim=10.0, caseB=None, tsw=np.inf, Q=1.0, socRef=None, group='static', **kw):
    return (name, alg, abl, seed, pack, case, caseB or case, soc, Q, tsw, Tend, socRef if socRef is not None else soc, Ilim, group, kw)
