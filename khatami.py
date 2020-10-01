import pandas as pd
import numpy as np
df=pd.read_csv('./Data/khatamiMni.csv')
#res['khatami_mni']=np.round(np.array(df['Mni']),3)
#res['khatami_mni_e']=np.round(np.array(df['Mni_e']),3)
#np.savetxt("Data/table1_26dec.csv", ab, fmt="%10s & %12s & %6s & %10.1f (%10.1f) & %10.4f (%10.4f) & %10.2f (%10.2f) & %10.1f (%10.1f)", newline=' \\\\\n')
#np.savetxt("Data/table2_26dec.csv", res, fmt="%10s & $%6s$ & %10.2f (%10.2f) & %10.1f (%10.1f) & %10.3f (%10.3f) & %10.1f (%10.1f) & %10.2f (%10.2f) & %10.2f (%10.2f) & %10.3f (%10.3f) ", newline=' \\\\\n')

print np.round(np.mean(np.array(df[df['type']=='IIb']['Mni'])),2),np.round(np.median(np.array(df[df['type']=='IIb']['Mni'])),2),np.round(np.std(np.array(df[df['type']=='IIb']['Mni'])),2)
print np.round(np.mean(np.array(df[df['type']=='Ib']['Mni'])),2),np.round(np.median(np.array(df[df['type']=='Ib']['Mni'])),2),np.round(np.std(np.array(df[df['type']=='Ib']['Mni'])),2)
print np.round(np.mean(np.array(df[df['type']=='Ic']['Mni'])),2),np.round(np.median(np.array(df[df['type']=='Ic']['Mni'])),2),np.round(np.std(np.array(df[df['type']=='Ic']['Mni'])),2)

print np.round(np.mean(np.array(df[df['type']=='Ic BL']['Mni'])),2),np.round(np.median(np.array(df[df['type']=='Ic BL']['Mni'])),2),np.round(np.std(np.array(df[df['type']=='Ic BL']['Mni'])),2)

print np.round(np.mean(np.array(df['Mni'])),2),np.round(np.median(np.array(df['Mni'])),2),np.round(np.std(np.array(df['Mni'])),2)

    #,np.median(np.round(np.array(df['Mni']),3)), np.std(np.round(np.array(df['Mni']),3))
#print np.mean(np.round(np.array(df['Mni'][df['type']=='IIb']),3)),np.median(np.round(np.array(df['Mni']),3)), np.std(np.round(np.array(df['Mni']),3))
print np.mean(np.abs(np.array(df['Mni']-df['Tail']))/np.array(df['Tail']))
print np.mean(np.abs(np.array(df['Arnett']-df['Tail']))/np.array(df['Tail']))

print np.sqrt(np.mean(np.square(np.array(df['Mni']-df['Tail']))))
print np.sqrt(np.mean(np.square(np.array(df['Arnett']-df['Tail']))))