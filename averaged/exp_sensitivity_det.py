import batch, pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
S=[(72,'C4',.6),(72,'C7',.95),(72,'C8',.6),(48,'C7',.6)]
jobs=[]
for si in range(11):
    for p,c,s in S: jobs.append(batch.J(f'sens{si}',6,4096*si,1,p,c,s))
    jobs.append(batch.J(f'sens{si}',6,4096*si,1,72,'C8',.6,Q=1.8e-4,socRef=.88,Tend=1.12,group='vrise'))
R=pd.DataFrame(batch.runall(jobs)); R.to_csv('t6.csv')
print(R.pivot_table(index='name',columns='scen',values='eff').round(3))
print(R.pivot_table(index='name',columns='scen',values='Ebat').round(1))
print(R.pivot_table(index='name',columns='scen',values='nRescan').round(1))
