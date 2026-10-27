import batch, pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
S=[(72,'C4',.6),(72,'C6',.6),(72,'C7',.95),(72,'C8',.6),(48,'C4',.6),(48,'C7',.6),(48,'C8',.6)]
A=[('RA-DPSO',6,0),('RA-DPSO-fast',6,2)]
jobs=[batch.J(n,a,b,1,p,c,s) for n,a,b in A for p,c,s in S]
jobs+=[batch.J(n,a,b,1,48,'C1',.93,Q=1e-3,socRef=.99,group='cv') for n,a,b in A+[('DPSO',5,0),('DPSO-PI-SS',5,513)]]
R=pd.DataFrame(batch.runall(jobs))
print(R[['name','scen','eff','Ebat','Ipk','tOver','dVcv','VbErr','success']].round(3).to_string())
