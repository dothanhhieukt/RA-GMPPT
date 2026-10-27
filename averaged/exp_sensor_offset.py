import batch, pandas as pd
J=[]
for r in (1,2,3):
    for p,c,s in [(72,'C7',.95),(72,'C8',.6)]:
        for nm,abl in [('D093',4096),('D093*',4096|8),('Voff1',3*4096),('Voff1*',3*4096|8)]:
            J.append(batch.J(nm,6,abl,r,p,c,s))
if __name__=='__main__':
    R=pd.DataFrame(batch.runall(J,2)); R.to_csv('t_sensr.csv')
    print(R.groupby(['name','scen'])[['Ebat','nRescan']].mean())
