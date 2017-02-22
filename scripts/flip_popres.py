import pandas as pd
import numpy as np
import os

"""
this file compares the Reich data with the Human Origins annotation file from
Affymetrix. As all SNP appear to be on the forward strand, this does not seem
to be necessary.
"""

def run(bimfile='raw/POPRES_Genotypes_QC1_v2.bim',
        chip1="chip/Mapping250K_Nsp.na32.annot.csv",
        chip2="chip/Mapping250K_Sty.na32.annot.csv",
        outname="raw/popres/keep_snp.txt",
        flip_file_name='raw/popres/flip.txt'):
    data = pd.read_table(bimfile, sep=r"\s*", header=None)
    data.columns = ['CHROM', 'Probe Set ID', 'X', 'POS', 'ALLELE1', 'ALLELE2']

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

    affy_file2 = pd.read_csv(chip1, comment='#', dtype='str')
    affy_file1 = pd.read_csv(chip2, comment='#', dtype='str')
    affy_file = pd.concat((affy_file1, affy_file2))
    affy_file.drop_duplicates(subset=("dbSNP RS ID", 'Probe Set ID'),
                              inplace=True)

    m_data = pd.merge(data, affy_file, how='inner')
                                                                              
    chr1 = pd.Series(m_data['CHROM'], dtype=str)
    chr2 = pd.Series(m_data['Chromosome'], dtype=str)
                                                                              
    m_data = m_data[chr1 == chr2]
    m_data = m_data[m_data['Strand'] != '---']
    m_data = m_data[m_data['Physical Position'] != '---']
    keep_file = m_data[['Probe Set ID', 'Physical Position']]
    keep_file.to_csv(outname, sep='\t', index=False, header=None)

    m_data = pd.merge(data, affy_file, how='inner')
                                                                              
    chr1 = pd.Series(m_data['CHROM'], dtype=str)
    chr2 = pd.Series(m_data['Chromosome'], dtype=str)
                                                                              
    m_data = m_data[chr1 == chr2]
    m_data = m_data[m_data['Strand'] == '-']
    flip_file = m_data['Probe Set ID']
    flip_file.to_csv(flip_file_name, sep='\t',
                     index=False, header=None)
try:
    bimfile = snakemake.input.bim
    chip1, chip2 = snakemake.input.chip
    keepname = snakemake.output.keep
    outbed = snakemake.output.bed
    flip_name = snakemake.output.flip
    print(chip1, chip2)
    run(bimfile, chip1, chip2, keepname, flip_name)
    s = 'plink --bfile %s --extract %s --update-map %s --make-bed --out %s --autosome'
    s = s % (os.path.splitext(bimfile)[0], keepname, keepname, os.path.splitext(outbed)[0])
    s += ' --flip %s ' % flip_name
    os.system(s)
except ValueError:
    print("error")


