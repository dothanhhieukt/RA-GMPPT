import batch, sim, pandas as pd, sys
SM=[(48,'C4',.6),(48,'C7',.2),(48,'C7',.6),(48,'C8',.6),(72,'C4',.6),(72,'C6',.6),(72,'C7',.6),(72,'C7',.95),(72,'C8',.2),(72,'C8',.6)]
J=[]
for n,a,b,N in [('DPSO',5,0,3),('RA-DPSO',6,0,3),('CS',7,0,5),('RA-CS',8,0,5)]:
    for r in range(1,N+1):
        for p,c,s in SM: J.append(batch.J(n,a,b,r,p,c,s))
for n,a,b in [('RA-DPSO',6,0),('RA-DPSO-F',6,2),('DPSO',5,0)]:
    for r in (1,2,3): J.append(batch.J(n,a,b,r,72,'C4',.6,Tend=1.1,caseB='C1',tsw=0.5,group='dynamic'))
    J.append(batch.J(n,a,b,1,48,'C1',.6))
if __name__=='__main__':
    tag=sys.argv[1]
    if tag=='new': sim.RD=1.0; sim.CD=470e-6
    R=pd.DataFrame(batch.runall(J,2)); R.to_csv(f't_damp_{tag}.csv')
