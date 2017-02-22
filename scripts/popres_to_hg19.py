import pandas as pd
import numpy as np

"""
this file compares the Reich data with the Human Origins annotation file from
Affymetrix. As all SNP appear to be on the forward strand, this does not seem
to be necessary.
"""

data = pd.read_table("tmp/popres0.bim",
                     sep=r"\s*", header=None)
data.columns = ['CHROM', 'Probe Set ID', 'X', 'POS', 'ALLELE1', 'ALLELE2']

autosomes = np.in1d(data['CHROM'], range(1, 23))
data = data[autosomes]

affy_file2 = pd.read_csv("chip/Mapping250K_Nsp.na32.annot.csv",
                         comment='#')
affy_file1 = pd.read_csv("chip/Mapping250K_Sty.na32.annot.csv",
                         comment='#')
affy_file = pd.concat((affy_file1, affy_file2))

m_data = pd.merge(data, affy_file, how='left')
                                                                          

m_data['to_flip'] = m_data['Strand'] == '-'

id_pos = zip(m_data['Chromosome'], m_data['Physical Position'])
new_ids = ["%s_%s" % i for i in id_pos]
m_data['ID2'] = new_ids

m_data[m_data['to_flip']].to_csv('raw/popres/to_flip.txt', sep='\t',
                                 columns=('ID2',), index=False, header=None)

opt = m_data[['Chromosome', 'ID2', 'X', 'Physical Position',
              'ALLELE1', 'ALLELE2']]
opt.to_csv("tmp/popres0.bim", index=False, header=None, sep="\t")
