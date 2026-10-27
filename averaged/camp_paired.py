import batch, pandas as pd, numpy as np, sys, os
MP=[(48,'C4',.6),(48,'C7',.2),(48,'C7',.6),(48,'C8',.6),(72,'C4',.6),(72,'C6',.6),(72,'C7',.6),(72,'C7',.95),(72,'C8',.2),(72,'C8',.6)]
V=dict(Q=1.8e-4,socRef=0.88,Tend=1.12,group='vrise')
part=sys.argv[1]
jobs=[]
if part=='static':
    for n,a,b in [('CS',7,0),('RA-CS',8,0),('CS-AR',7,256),('RA-CS-noRA',8,16),('RA-CS-R1',8,192)]:
        for p,c,s in MP:
            for r in range(1,31): jobs.append(batch.J(n,a,b,r,p,c,s))
    for n,a,b in [('P&O',1,0),('DPSO',5,0),('DPSO-AR',5,256),('RA-DPSO',6,0),('SCAN',9,0),('SCAN-AR',9,256),('RA-SCAN',10,0)]:
        for p,c,s in MP:
            for r in range(1,11): jobs.append(batch.J(n,a,b,r,p,c,s))
elif part=='vrise':
    for n,a,b,N in [('CS',7,0,30),('RA-CS',8,0,30),('RA-CS-R1',8,192,30),('CS-AR',7,256,30),('DPSO',5,0,10),('RA-DPSO',6,0,10),('RA-DPSO-R1',6,192,10),('DPSO-AR',5,256,10),('SCAN',9,0,10),('RA-SCAN',10,0,10),('RA-SCAN-R1',10,192,10)]:
        for r in range(1,N+1): jobs.append(batch.J(n,a,b,r,72,'C8',.6,**V))
elif part=='dyn':
    for n,a,b,N in [('CS',7,0,30),('RA-CS',8,0,30),('RA-CS-F',8,2,30),('DPSO',5,0,10),('RA-DPSO',6,0,10),('RA-DPSO-F',6,2,10),('P&O',1,0,10),('DPSO-PI-SS',5,513,10)]:
        for c1,c2 in [('C4','C1'),('C1','C4'),('C2','C4'),('C4','C2')]:
            for r in range(1,N+1): jobs.append(batch.J(n,a,b,r,72,c1,.6,Tend=1.1,caseB=c2,tsw=0.5,group='dynamic'))
out=f'camp_{part}.csv'
done=set()
if os.path.exists(out):
    D=pd.read_csv(out); done=set(zip(D.name,D.scen,D.seed))
def key(j):
    name,alg,abl,seed,pack,cA,cB,soc0,Q,tsw,Tend,socRef,Ilim,group,kw=j
    sc=f'{pack}-{cA}{">"+cB if cB!=cA else ""}-S{round(100*soc0):02d}'+(f'-Q{Q:g}' if Q<1 else '')
    return (name,sc,seed)
jobs=[j for j in jobs if key(j) not in done]
print(part,len(jobs),'to run',flush=True)
import multiprocessing as mp
with mp.Pool(2) as pool:
    buf=[]
    for i,M in enumerate(pool.imap_unordered(batch.job,jobs,chunksize=2)):
        buf.append(M)
        if len(buf)>=40 or i==len(jobs)-1:
            pd.DataFrame(buf).to_csv(out,mode='a',header=not os.path.exists(out),index=False); buf=[]
            print(i+1,flush=True)
