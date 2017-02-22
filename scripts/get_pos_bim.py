import pandas as pd


def unify_bim_id(data):
    data.columns = ['CHROM', 'ID', 'X', 'POS', 'ALLELE1', 'ALLELE2']
    s1 = data['CHROM']
    s2 = data['POS']
    s = zip(s1, s2)
    new_names = ["%s_%s" % i for i in s] 
    data['ID'] = new_names 
    return data

if __name__ == '__main__':
    import sys
    import os
    if len(sys.argv) > 2:
        data = pd.read_table(sys.argv[1] + ".bim", sep=r"\s*", header=None)
        data = unify_bim_id(data)
        data.to_csv(sys.argv[2] + ".bim", index=False, header=None, sep=" ")
#        os.system("ln %s.bed %s.bed" % (sys.argv[1], sys.argv[2]))
#        os.system("ln %s.fam %s.fam" % (sys.argv[1], sys.argv[2]))
