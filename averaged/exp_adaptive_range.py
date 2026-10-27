import batch, pandas as pd
pd.set_option('display.width',200)
S=[(72,'C4',.6),(72,'C6',.6),(72,'C7',.6),(72,'C7',.95),(72,'C8',.6),(72,'C8',.2),(48,'C4',.6),(48,'C7',.6),(48,'C8',.6),(48,'C7',.2)]
A=[('DPSO',5,0),('DPSO-AR',5,256),('RA-DPSO',6,0),('SCAN',9,0),('SCAN-AR',9,256),('RA-SCAN',10,0),('P&O',1,0)]
jobs=[batch.J(n,a,b,1,p,c,s) for n,a,b in A for p,c,s in S]
R=pd.DataFrame(batch.runall(jobs))
print(R.pivot_table(index='scen',columns='name',values='eff').round(3))
print(R.pivot_table(index='scen',columns='name',values='Ebat').round(1))
R.to_csv('t1.csv')
