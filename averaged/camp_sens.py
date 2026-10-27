import batch, pandas as pd, numpy as np, os, multiprocessing as mp
S=[(72,'C4',.6,{}),(72,'C7',.95,{}),(72,'C8',.6,{}),(48,'C7',.6,{}),(72,'C8',.6,dict(Q=1.8e-4,socRef=0.88,Tend=1.12,group='vrise'))]
jobs=[]
for si in range(11):
    for p,c,s,kw in S:
        for r in range(1,31): jobs.append(batch.J(f'RA-CS-s{si}',8,4096*si,r,p,c,s,**kw))
        for r in range(1,11): jobs.append(batch.J(f'RA-DPSO-s{si}',6,4096*si,r,p,c,s,**kw))
out='camp_sens.csv'
done=set()
if os.path.exists(out):
    D=pd.read_csv(out); done=set(zip(D.name,D.scen,D.seed))
def key(j):
    name,alg,abl,seed,pack,cA,cB,soc0,Q,tsw,Tend,socRef,Ilim,group,kw=j
    return (name,f'{pack}-{cA}-S{round(100*soc0):02d}'+(f'-Q{Q:g}' if Q<1 else ''),seed)
jobs=[j for j in jobs if key(j) not in done]
print(len(jobs),flush=True)
with mp.Pool(2) as pool:
    buf=[]
    for i,M in enumerate(pool.imap_unordered(batch.job,jobs,chunksize=4)):
        buf.append(M)
        if len(buf)>=40 or i==len(jobs)-1:
            pd.DataFrame(buf).to_csv(out,mode='a',header=not os.path.exists(out),index=False); buf=[]; print(i+1,flush=True)
