import batch, pandas as pd
pd.set_option('display.width',250)
S=[(72,'C4',.6),(72,'C6',.6),(72,'C7',.6),(72,'C7',.95),(72,'C8',.6),(72,'C8',.2),(48,'C4',.6),(48,'C7',.6),(48,'C8',.6),(48,'C7',.2),(48,'C1',.6),(72,'C1',.6)]
A=[('DPSO',5,0),('DPSO-20k',5,1024),('RA-DPSO',6,0),('RA-DPSO-20k',6,1024)]
R=pd.DataFrame(batch.runall([batch.J(n,a,b,1,p,c,s) for n,a,b in A for p,c,s in S]))
for m in ['eff','Ebat','Ipk','IbErr']: print(m); print(R.pivot_table(index='scen',columns='name',values=m).round(3))
