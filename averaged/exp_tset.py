import batch, pandas as pd
S=[(48,'C1',.6,10),(72,'C1',.6,10),(48,'C3',.2,10),(48,'C2',.6,10),(48,'C4',.6,6),(72,'C4',.6,4),(72,'C6',.6,3.5)]
A=[('DPSO',5,0),('DPSO-SS',5,512),('DPSO-PI',5,1),('DPSO-PI-SS',5,513),('RA-DPSO',6,0),('RA-DPSO-fast',6,2)]
jobs=[batch.J(n,a,b,1,p,c,s,Ilim=I) for n,a,b in A for p,c,s,I in S]
for n,a,b,N in [('CS',7,0,30),('RA-CS',8,0,30),('RA-CS-F',8,2,30),('DPSO',5,0,10),('RA-DPSO',6,0,10),('RA-DPSO-F',6,2,10),('P&O',1,0,10),('DPSO-PI-SS',5,513,10)]:
    for r in range(1,N+1): jobs.append(batch.J(n,a,b,r,72,'C4',.6,Tend=1.1,caseB='C1',tsw=0.5,group='dynamic'))
if __name__=='__main__':
    R=pd.DataFrame(batch.runall(jobs)); R.to_csv('t_tset.csv')
