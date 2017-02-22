import numpy as np
import os
import pandas as pd

#a = pd.read_table("raw/xing/affy6_344_raw_genotype_xing", index_col=0,
#                  skiprows=1, header=None)
def make_fam(inds):
    """
    Individual's family ID ('FID')
    Individual's within-family ID ('IID'; cannot be '0')
    Within-family ID of father ('0' if father isn't in dataset)
    Within-family ID of mother ('0' if mother isn't in dataset)
    Sex code ('1' = male, '2' = female, '0' = unknown)
    Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing
    data if case/control)
    """
    nrow = inds.shape[0]
    z = np.zeros(nrow, dtype='int')
    fam = np.vstack((inds, inds, z, z,z,z)).transpose()

    return fam

def run(infile="raw/affy6_344_raw_genotype_xing",
        chip="chip/GenomeWideSNP_6.na32.annot.csv",
        to_flip='tmp/to_flip.txt',
        out_tped='x',
        out_tfam='y'):
    a = pd.read_table(infile, usecols=[0], skiprows=2, header=None)
    #affy_file = pd.read_csv(chip, comment='#', engine="python")
    affy_file = pd.read_csv(chip, skiprows=21)
    print("loaded chip")
    affy_file.drop_duplicates(subset=("dbSNP RS ID"),
                              inplace=True)
    a.columns = ['Probe Set ID']
    b = pd.merge(a, affy_file, how='inner')
    b = b[b['Strand'] != '---']

    b['to_flip'] = b['Strand'] == '-'

    id_pos = zip(b['Chromosome'], b['Physical Position'])
    new_ids = ["%s_%s" % i for i in id_pos]
    b['ID2'] = new_ids

    b[b['to_flip']].to_csv(to_flip, sep='\t',
                                     columns=('ID2',), index=False, header=None)

    b['Z'] = 0
    d = zip(b['Chromosome'],b['ID2'],b['Z'], b['Physical Position'])
    d_str = ["%s\t%s\t%s\t%s" % i for i in d]

    al = zip(b['Allele A'], b['Allele B'])

    keys = b['Probe Set ID']

    bdict = dict((k, v) for (k, v) in zip(keys, d_str))
    allele_dict = dict((k, v) for (k, v) in zip(keys, al))

    with open(out_tped, 'w') as tped:
        with open(infile) as reader:
            for row in reader:
                if row.startswith("#"):
                    continue
                row = row.split()
                snp_id = row[0]
                if snp_id not in bdict:
                    continue
                a1, a2 = allele_dict[snp_id]
                front = bdict[snp_id]
                allele_string = front
                for gt in row[1:]:
                    if gt == '-1':
                        allele_string += ' 0 0'
                    elif gt == '0':
                        allele_string += ' %s %s' % (a1, a1)
                    elif gt == '1':
                        allele_string += ' %s %s' % (a1, a2)
                    elif gt == '2':
                        allele_string += ' %s %s' % (a2, a2)
                allele_string += '\n'
                tped.write(allele_string)

    inds = list(pd.read_table(infile, nrows=1,
                              skiprows=1))

    fam = make_fam(np.array(inds[1:]))
    np.savetxt(out_tfam, fam, '%s')

try:
    infile = snakemake.input.file
    chip = snakemake.input.chip
    tped = snakemake.output.tped
    tfam = snakemake.output.tfam
    flip = snakemake.output.flip
    outbed = snakemake.output.bed
    run(infile, chip, flip, tped, tfam)
    s =  'plink --tfile %s --flip %s --make-bed --out %s --autosome '
    s = s % (os.path.splitext(tped)[0], flip, os.path.splitext(outbed)[0])
    os.system(s)

except NameError as E:
    missing_var = E.args[0].split("'")[1] 
    if missing_var != 'snakemake': raise E





