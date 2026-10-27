import batch, pandas as pd
SM=[(48,'C4',.6),(48,'C7',.2),(48,'C7',.6),(48,'C8',.6),(72,'C4',.6),(72,'C6',.6),(72,'C7',.6),(72,'C7',.95),(72,'C8',.2),(72,'C8',.6)]
J=[]
for n,a,b in [('DPSO',5,0),('RA-DPSO',6,0),('DPSO-AR',5,256)]:
    for r in (1,2,3):
        for p,c,s in SM: J.append(batch.J(n,a,b,r,p,c,s)); J.append(batch.J(n+'*',a,b|8,r,p,c,s))
for n,a,b in [('DPSO',5,0),('RA-DPSO',6,0),('RA-DPSO-F',6,2)]:
    for r in (1,2,3):
        for c1,c2 in [('C4','C1'),('C1','C4'),('C2','C4'),('C4','C2')]:
            for nn,bb in [(n,b),(n+'*',b|8)]: J.append(batch.J(nn,a,bb,r,72,c1,.6,Tend=1.1,caseB=c2,tsw=0.5,group='dynamic'))
V=dict(Q=1.8e-4,socRef=0.88,Tend=1.12,group='vrise')
for n,a,b in [('DPSO',5,0),('RA-DPSO',6,0),('DPSO-AR',5,256)]:
    for r in (1,2,3):
        for nn,bb in [(n,b),(n+'*',b|8)]: J.append(batch.J(nn,a,bb,r,72,'C8',.6,**V))
if __name__=='__main__':
    R=pd.DataFrame(batch.runall(J,2)); R.to_csv('t_rskip.csv')
