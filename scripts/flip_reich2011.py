import pandas as pd
import os
import numpy as np

"""
this file compares the Reich data with the Human Origins annotation file from
Affymetrix. As all SNP appear to be on the forward strand, this does not seem
to be necessary.
"""

def run(bim='raw/reich2011/Denisova-SEAsia-Oceania.bim',
        chip='chip/GenomeWideSNP_6.na32.annot.csv',
        keep_file_name="tmp/reich2011.keep"):
    data = pd.read_table(bim, sep=r"\s*", header=None, engine='python')
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


    affy_file = pd.read_csv(chip, comment='#')
    affy_file.drop_duplicates(subset=("dbSNP RS ID"),
                              inplace=True)

    m_data = pd.merge(data, affy_file, how='inner')
                                                                              
    chr1 = pd.Series(m_data['CHROM'], dtype=str)
    chr2 = pd.Series(m_data['Chromosome'], dtype=str)
                                                                              
    m_data = m_data[chr1 == chr2]
    m_data = m_data[m_data['Strand'] != '---']
    m_data = m_data[m_data['Physical Position'] != '---']
    keep_file = m_data[['dbSNP RS ID', 'Physical Position']]
    keep_file.to_csv(keep_file_name, sep='\t',
                     index=False, header=None)

try:
    inbim = snakemake.input.bim
    chip = snakemake.input.chip
    keep = snakemake.output.keep
    outbed = snakemake.output.bed
    run(inbim, chip, keep)
    s =  'plink --bfile %s --extract %s --update-map %s --make-bed --out %s --autosome '
    s = s % (os.path.splitext(inbim)[0], keep, keep, os.path.splitext(outbed)[0])
    os.system(s)
except NameError as E:
    missing_var = E.args[0].split("'")[1] 
    if missing_var != 'snakemake': raise E
