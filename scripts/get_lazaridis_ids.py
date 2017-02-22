import numpy as np
import pandas as pd

t = pd.read_excel("raw/EuropeAllData/10_13_2014 Excel version of " +
                  "Table S9.4.xlsx")
pops = t['Verbose Population ID']

q = pd.read_table("raw/EuropeAllData/vdata.ind", sep=r"\s*", header=None)
q2 = q[[0, 2]]
q2['data'] = 'human_origins'
q2['x'] = 'private'

kp = np.in1d(q2[2], pops)

q2 = q2[kp]
q2.to_csv("raw/EuropeAllData/vdata_inds.txt",
          sep="\t", index=False, header=False)                 



