import pandas as pd
import numpy as np
from valenti_ni56 import *
df=pd.read_csv('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata_wygoda_26dec.csv')


sug_beta={'Ib':0.79, 'Ic':0.82, 'Ic BL':0.53, 'IIb':0.76, 'Ib pec':0.79, 'Ia':1.6}
df_save = pd.DataFrame(columns=['name', 'Mni','Mni_e', 'Arnett','Tail'])
for index,sn in enumerate(df['name'].tolist()):
    print index
    # print sn,float(df[df['name'] == sn]['peakL'].tolist()[0].split(";")[0])
    Mni,mni_e=nickel_mass_khatami(float(df[df['name'] == sn]['peakt'].tolist()[0].split(";")[0]), 10**float(df[df['name'] == sn]['peakL'].tolist()[0].split(";")[0]), 10**(float(df[df['name'] == sn]['peakL'].tolist()[0].split(";")[1])), sug_beta[df[df['name'] == sn]['sn_type_full'].tolist()[0]])
    df_save.loc[index] = [sn,Mni,mni_e,float(df[df['name'] == sn]['arnett_ni_mass'].tolist()[0].split(";")[0]),float(df[df['name'] == sn]['tail_ni_mass'].tolist()[0].split(";")[0])]
print df_save
print np.sqrt(np.sum((df_save['Mni']-df_save['Tail'])**2))/28, np.sqrt(np.sum((df_save['Arnett']-df_save['Tail'])**2))/28
df_save.to_csv('./Data/khatamiMni.csv')