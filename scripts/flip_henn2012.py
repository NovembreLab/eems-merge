import pandas as pd
import numpy as np
import os

"""
this file compares the Reich data with the Human Origins annotation file from
Affymetrix. As all SNP appear to be on the forward strand, this does not seem
to be necessary.
"""

data = pd.read_table("qatari/NWAfrica_HM3_Qat.bim",
                     sep=r"\s*", header=None)
data.columns = ['CHROM', 'dbSNP RS ID', 'X', 'POS', 'ALLELE1', 'ALLELE2']

data['keep'] = True
rm1 = np.logical_and(data['ALLELE1'] == 'A', data['ALLELE2'] == 'T')
rm2 = np.logical_and(data['ALLELE1'] == 'T', data['ALLELE2'] == 'A')
rm3 = np.logical_and(data['ALLELE1'] == 'C', data['ALLELE2'] == 'G')
rm4 = np.logical_and(data['ALLELE1'] == 'G', data['ALLELE2'] == 'C')
data['keep'][rm1] = False
data['keep'][rm2] = False
data['keep'][rm3] = False
data['keep'][rm4] = False

data = data[data['keep']]

autosomes = np.in1d(data['CHROM'], range(1, 23))
data = data[autosomes]

affy_file = pd.read_csv("chip/GenomeWideSNP_6.na32.annot.csv",
                        comment='#')

affy_file.drop_duplicates(subset=("dbSNP RS ID"),
                          inplace=True)

m_data = pd.merge(data, affy_file, how='inner')
                                                                          
chr1 = pd.Series(m_data['CHROM'], dtype=str)
chr2 = pd.Series(m_data['Chromosome'], dtype=str)
                                                                          
m_data = m_data[chr1 == chr2]
m_data = m_data[m_data['Strand'] != '---']
keep_file = m_data[['dbSNP RS ID', 'Physical Position']]
keep_file.to_csv("qatari/keep_snp.txt", sep='\t', index=False, header=None)

m_data = pd.merge(data, affy_file, how='inner')
                                                                          
chr1 = pd.Series(m_data['CHROM'], dtype=str)
chr2 = pd.Series(m_data['Chromosome'], dtype=str)
                                                                          
m_data = m_data[chr1 == chr2]
m_data = m_data[m_data['Strand'] == '-']
flip_file = m_data['dbSNP RS ID']
flip_file.to_csv('qatari/flip.txt', sep='\t',
                 index=False, header=None)

s = 'plink --bfile %s --extract %s --update-map %s --make-bed --out %s --autosome'
s = s % ("qatari/NWAfrica_HM3_Qat", 'qatari/keep_snp.txt', 'qatari/keep_snp.txt', 'qatari/hg37')
s += ' --flip %s ' % 'qatari/flip.txt'
print(s)
os.system(s)
