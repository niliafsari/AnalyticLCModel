import numpy as np
import pandas as pd

df_old=pd.read_csv("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata_wygoda_26dec.csv")
df_new=pd.read_csv("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata_wygoda_26dec_tpeak.csv")

df_old=df_old[['name','tail_ni_mass','beta_req']]
#df_old[['tail_ni_mass']]=df_old[['tail_ni_mass']]
new=df_old["tail_ni_mass"].str.split(";", n = 1, expand = True)
df_old['tail_ni_mass']=new[0]
new=df_old["beta_req"].str.split(";", n = 1, expand = True)
df_old['beta_req']=new[0]

df_new=df_new[['name','tail_ni_mass','beta_req']]
#df_new[['tail_ni_mass']]=df_new[['tail_ni_mass']]
new=df_new["tail_ni_mass"].str.split(";", n = 1, expand = True)
df_new['tail_ni_mass']=new[0]
new=df_new["beta_req"].str.split(";", n = 1, expand = True)
df_new['beta_req']=new[0]

df_old.to_csv('Data/old_df.csv', index=False)
df_new.to_csv('Data/new_df.csv', index=False)