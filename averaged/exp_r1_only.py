import batch, pandas as pd
S6=[(48,'C7',.6),(48,'C8',.6),(72,'C4',.6),(72,'C7',.6),(72,'C7',.95),(72,'C8',.6)]
R=pd.DataFrame(batch.runall([batch.J('RA-DPSO-R1',6,192,r,p,c,s) for p,c,s in S6 for r in range(1,11)]))
g=R.groupby('scen'); print(g.eff.mean().mean(), R.eff.min(), g.tUnreach.mean().mean()*1e3, g.Ebat.mean().mean())
