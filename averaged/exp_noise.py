import batch, pandas as pd
SM=[(48,'C4',.6),(48,'C7',.2),(48,'C7',.6),(48,'C8',.6),(72,'C4',.6),(72,'C6',.6),(72,'C7',.6),(72,'C7',.95),(72,'C8',.2),(72,'C8',.6)]
J=[]
for f in (0.5,2.0):
    for n,a,b,N in [('DPSO',5,0,3),('RA-DPSO',6,0,3),('DPSO-AR',5,256,3),('CS',7,0,10),('RA-CS',8,0,10)]:
        for r in range(1,N+1):
            for p,c,s in SM: J.append(batch.J(f'{n}|{f}',a,b,r,p,c,s,noise=(0.05*f,0.005*f)))
if __name__=='__main__':
    R=pd.DataFrame(batch.runall(J,2)); R.to_csv('t_noise.csv')
