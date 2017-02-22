import pandas as pd
import numpy as np

"""
this file compares the Reich data with the Human Origins annotation file from
Affymetrix. As all SNP appear to be on the forward strand, this does not seem
to be necessary.
"""

data = pd.read_table("tmp/lazaridis0.bim", sep=r"\s*", header=None)
data.columns = ['CHROM', 'Affy SNP ID', 'X', 'POS', 'ALLELE1', 'ALLELE2']

autosomes = np.in1d(data['CHROM'], range(1, 23))
data = data[autosomes]

affy_file = pd.read_csv("chip/Axiom_GW_HuOrigin.na35.annot.csv",
                        comment='#')
affy_file.drop_duplicates(subset=("dbSNP RS ID", 'Affy SNP ID', 'Ref Allele'),
                          inplace=True)

m_data = pd.merge(data, affy_file, how='inner')
type1 = m_data['ALLELE1'] == m_data['Ref Allele']
type2 = m_data['ALLELE2'] == m_data['Ref Allele']


type1e = np.logical_and(m_data['ALLELE1'] == m_data['Ref Allele'],
                        m_data['ALLELE2'] == m_data['Alt Allele'])
type2e = np.logical_and(m_data['ALLELE2'] == m_data['Ref Allele'],
                        m_data['ALLELE1'] == m_data['Alt Allele'])
                                                                          
weird = np.logical_not(np.logical_or(type1e, type2e))
m_data = m_data[np.logical_not(weird)]
                                                                          
alleles = ['ALLELE1', 'ALLELE2', 'Allele A', 'Allele B',
           'Ref Allele', 'Alt Allele']
                                                                          
chr1 = pd.Series(m_data['CHROM'], dtype=str)
chr2 = pd.Series(m_data['Chromosome'], dtype=str)
                                                                          
m_data = m_data[chr1 == chr2]


m_data['to_flip'] = m_data['Strand'] == '-'

row_no_id = m_data['dbSNP RS ID'] == '---'
no_id_pos = zip(m_data[row_no_id]['CHROM'], m_data[row_no_id]['POS'])
new_ids = ["%s_%s" % i for i in no_id_pos]
m_data['dbSNP RS ID'][row_no_id] = new_ids

opt = m_data[['Affy SNP ID', 'dbSNP RS ID', 'CHROM', 'POS', 'Ref Allele',
              'to_flip']]
opt.to_csv("reich_clean.csv", index=False)
