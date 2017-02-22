import numpy as np
import pandas as pd


def make_fam(inds):
    """
    Individual's family ID ('FID') 
    Individual's within-family ID ('IID'; cannot be '0')
    Within-family ID of father ('0' if father isn't in dataset)
    Within-family ID of mother ('0' if mother isn't in dataset)
    Sex code ('1' = male, '2' = female, '0' = unknown)' ' '
    Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing
    data if case/control)
    """
    nrow = inds.shape[0]
    z = np.zeros(nrow, dtype='int')
    fam = np.vstack((inds, inds, z, z,z,z)).transpose()
    return(fam)


#delete stuff
def del_row(data, rows):
    return np.delete(data, rows, axis=0)

#tped = to_tped(chrom_l, rs_id, pos_l, alleles, data)
def to_tped(chrom, rs_id, pos, alleles, data):
    shape = data.shape[0], data.shape[1] * 2 + 4
    tped = np.ones(shape, dtype=object)
    tped[:, 0] = chrom
    tped[:, 1] = rs_id
    tped[:, 3] = pos
    alleles = np.array(alleles, object)
    a1, a2 = np.array([a.split("/") for a in alleles]).transpose()
    for i, row in enumerate(data):
        if i % 10000 == 0: print(i)
        for j, loc in enumerate(row):
            jj = 4 + 2 * j
            jjj = 5 + 2 * j
            if loc == '9':
                res = 0,0
            elif loc == '1':
                res = a1[i], a2[i]
            elif loc == '2':
                res = a2[i], a2[i]
            elif loc == '0':
                res = a1[i], a1[i]
            else:
                res = 0, 0
            tped[i, jj] = res[0]
            tped[i, jjj] = res[1]

    return tped


def run_all(raw_genotypes="raw/hugo/Genotypes_All.txt",
            lifted_in='supplementary/lifted.bed',
            unlifted_in='supplementary/unlifted.bed',
            out_tped='hugo.tped', out_fam='hugo.fam'):
    q = pd.read_table(raw_genotypes, header=None, dtype='str')
    a = np.array(q)
    print('loaded')

    inds = a[0, 5:]
    affy_id = a[1:, 0]
    rs_id = a[1:, 1]
    chrom = a[1:, 2]
    pos_hg36 = a[1:, 3]
    alleles = a[1:, 4]

    data = a[1:, 5:]


    lifted = np.genfromtxt(lifted_in, dtype=str)
    unlifted = np.genfromtxt(unlifted_in, dtype=str)


    l_hg36 = pos_hg36.tolist()
    unlifted_rows = [l_hg36.index(u) for u in unlifted[:,1]]

    pos_hg36 = del_row(pos_hg36, unlifted_rows)
    chrom_l = np.array([s[3:] for s in lifted[:,0]])        
    pos_l = np.array(lifted[:,1])

    rs_id = del_row(rs_id, unlifted_rows)
    affy_id = del_row(affy_id, unlifted_rows)
    alleles = del_row(alleles, unlifted_rows)
    data = del_row(data, unlifted_rows)



    print('making tped..')
    tped = to_tped(chrom_l, rs_id, pos_l, alleles, data)
    print('made tped')
    np.savetxt(out_tped, tped, fmt="%s")
    del tped


    fam = make_fam(inds)
    np.savetxt(out_fam, fam, fmt='%s')



try:
    run_all(snakemake.input.raw, snakemake.input.lifted,
            snakemake.input.unlifted,
            snakemake.output.tped, snakemake.output.fam)
except NameError:
    print('error')
    pass



