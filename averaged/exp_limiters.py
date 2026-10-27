import batch, pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
S=[(48,'C1',.6,10),(72,'C1',.6,10),(48,'C3',.2,10),(48,'C2',.6,10),(48,'C4',.6,6),(72,'C4',.6,4),(72,'C6',.6,3.5)]
A=[('DPSO',5,0),('DPSO-SS',5,512),('DPSO-PI',5,1),('DPSO-PI-SS',5,513),('RA-DPSO',6,0),('RA-DPSO-fast',6,2),('RA-DPSO-conv',6,32),('RA-DPSO-conv-SS',6,544),('RA-DPSO-PI-SS',6,545)]
jobs=[batch.J(n,a,b,1,p,c,s,Ilim=I) for n,a,b in A for p,c,s,I in S]
R=pd.DataFrame(batch.runall(jobs)); R.to_csv('t2.csv')
for m in ['Ipk','tOver','tset','Ebat']:
    print(m); print(R.pivot_table(index='scen',columns='name',values=m).round(3))
