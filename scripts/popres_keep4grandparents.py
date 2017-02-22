import pandas as pd
import numpy as np
from collections import Counter

p = pd.read_table("intermediate/POPRES_SINGLETAB.csv")
same_grandma = (p.COUNTRY_PGM == p.COUNTRY_MGM)
same_grandpa = (p.COUNTRY_PGF == p.COUNTRY_MGF)
same_maternal = (p.COUNTRY_MGF == p.COUNTRY_MGM)
same_self_mat = p.COUNTRY_SELF == p.COUNTRY_MGM

same_gf1 = np.logical_and(same_grandma, same_grandpa)
same_gf2 = np.logical_and(same_maternal, same_self_mat)
same = np.logical_and(same_gf1, same_gf2)

pp = p[same]
p = pp[['SUBJID', 'COUNTRY_SELF', 'PRIMARY_LANGUAGE']]

sf = np.logical_and(p.COUNTRY_SELF == 'Switzerland',
            p.PRIMARY_LANGUAGE == 'French')
si = np.logical_and(p.COUNTRY_SELF == 'Switzerland',
            p.PRIMARY_LANGUAGE == 'Italian')
sg = np.logical_and(p.COUNTRY_SELF == 'Switzerland',
            p.PRIMARY_LANGUAGE == 'German')

p.loc[sf,'COUNTRY_SELF'] = 'Swiss-French'
p.loc[sg,'COUNTRY_SELF'] = 'Swiss-German'
p.loc[si,'COUNTRY_SELF'] = 'Swiss-Italian'


p.columns = 'SID', 'Country', 'Language'

a = pd.read_csv("regions/location_simplified.csv")
#a = a[a.Source == 'Popres']
a = a.iloc[:, [0,1]]
a.columns = 'PID', 'Country'

g = pd.merge(p, a, how='left')

popid_remove = [803, 809]
for i in popid_remove:
    g = g[g.PID != i]

g = g[['SID', 'PID']]
g = g[g.PID.notnull()]
g.PID = g.PID.astype(int)

g['Source'] = 'POPRES'
g['Permission'] = 'Private'
g.to_csv("intermediate/sample_pop_popres.csv", index=False)


