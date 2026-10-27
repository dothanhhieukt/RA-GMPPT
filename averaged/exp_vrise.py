import batch, pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
V=dict(Q=1.8e-4,socRef=0.88,Tend=1.12)
A=[('DPSO',5,0,1),('RA-DPSO',6,0,1),('RA-DPSO-R1',6,192,1),('SCAN',9,0,1),('RA-SCAN',10,0,1),('RA-SCAN-R1',10,192,1),('DPSO-AR',5,256,1)]+[(n,a,b,r) for n,a,b in [('CS',7,0),('RA-CS',8,0),('RA-CS-R1',8,192)] for r in range(1,11)]
jobs=[batch.J(n,a,b,r,72,'C8',.6,group='vrise',**V) for n,a,b,r in A]
R=pd.DataFrame(batch.runall(jobs)); R.to_csv('t5.csv')
print(R.groupby('name')[['eff','Ebat','nRescan','Vss']].agg(['mean','min','max']).round(3))
