import pandas as pd
import numpy as np

"""
this file compares the Hugo data with the Human Origins annotation file from
Affymetrix. As all SNP appear to be on the forward strand, this does not seem
to be necessary.
"""

data = pd.read_table("../hugo/hugo.bim", sep=r"\s*", header=None)
data.columns = ['CHROM', 'dbSNP RS ID', 'X', 'POS', 'ALLELE1', 'ALLELE2']

autosomes = np.in1d(data['CHROM'], range(1, 23))
data = data[autosomes]

affy_file = pd.read_csv("../chipfiles/Mapping50K_Xba240.na32.annot.csv",
                        comment='#')

m_data = pd.merge(data, affy_file, how='inner')
                                                                          
chr1 = pd.Series(m_data['CHROM'], dtype=str)
chr2 = pd.Series(m_data['Chromosome'], dtype=str)
                                                                          
m_data = m_data[chr1 == chr2]


m_data['to_flip'] = m_data['Strand'] == '-'

opt = m_data[['dbSNP RS ID', 'CHROM', 'POS', 
              'to_flip']]
opt.to_csv("hugo_clean.csv", index=False)
