"""script replaced by snakerule sample_pop_estonia"""

import pandas as pd

a = pd.read_excel("sources/Data_for_Ben_Meta.xlsx")
b = pd.read_csv("regions/location_simplified.csv")
a = a.fillna('')
a[a=='Belorussia'] = 'Belarus'
a[a=='Belarusia'] = 'Belarus'

b = b.fillna('')
s1 = set(a.columns)
s2 = set(b.columns)
indices = s1.intersection(s2)
indices.discard('Latitude')
indices.discard('Longitude')
indices.discard('Country of origin / collected in')

c = pd.merge(a, b, on=list(indices), how='left')
d = c[['Sample ID', 'ID']]

d.to_csv("intermediate/sample_pop_estonians.csv")


