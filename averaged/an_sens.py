import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
D=pd.read_csv('camp_sens.csv')
D['row']=D.name.str.extract(r's(\d+)$').astype(int); D['core']=D.name.str.extract(r'^(RA-\w+)-s')
lab={0:'nominal',1:'D̂max 0.93',2:'D̂max 0.97',3:'Vbat offset +1 V',4:'Vbat offset −1 V',5:'ΔVc 0',6:'ΔVc 1.5 V',7:'ΔVR 1 V',8:'ΔVR 4 V',9:'R2 filter 50 ms',10:'R2 on OCV estimate'}
for core in ['RA-DPSO','RA-CS']:
    X=D[D.core==core]
    E=X.pivot_table(index='row',columns='scen',values='Ebat',aggfunc='mean')
    eff=X.pivot_table(index='row',columns='scen',values='eff',aggfunc='mean')
    n=X.pivot_table(index='row',columns='scen',values='nRescan',aggfunc='mean')
    dE=100*(E/E.loc[0]-1)
    print(core); print('eff'); print(eff.round(3)); print('dE%'); print(dE.round(2)); print('nRescan'); print(n.round(2))
