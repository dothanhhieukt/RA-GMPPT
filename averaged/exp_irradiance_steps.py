import batch, pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
S=[('C4','C1'),('C1','C4'),('C2','C4'),('C4','C2')]
A=[('DPSO',5,0,1),('DPSO-PI-SS',5,513,1),('RA-DPSO',6,0,1),('RA-DPSO-fast',6,2,1)]+[('CS',7,0,r) for r in range(1,6)]+[('RA-CS',8,0,r) for r in range(1,6)]+[('RA-CS-fast',8,2,r) for r in range(1,6)]
jobs=[batch.J(n,a,b,r,72,c1,.6,Tend=1.1,caseB=c2,tsw=0.5) for n,a,b,r in A for c1,c2 in S]
R=pd.DataFrame(batch.runall(jobs)); R.to_csv('t3.csv')
for m in ['Ipk','tOver','Qover','Ebat','Pss']:
    print(m); print(R.pivot_table(index='scen',columns='name',values=m,aggfunc='mean').round(3))
