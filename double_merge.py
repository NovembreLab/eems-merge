import os


def plink_merge(in1, in2, out, tmp='tmp/merge'):
    """ uses plink to merge datasets

    this is hackish, but works: it runs plink three times to
        merge two datasets. 
        
        RUN1: naive merge. Will work if no SNP are flipped,
            otherwise the flip-file will be created.

        (if flipped SNP present)
        RUN2a : flip second input file
        RUN3a : merge in1 with flipped file

        (if no flipped SNP present)
        RUN2a : flip will fail, 
        RUN3a : merge in1 with flipped file
        
    """
    plink_merge_cmd = 'plink --bfile %s --bmerge %s --out %s --make-bed'
    plink_flip_cmd = 'plink --bfile %s --out %s --make-bed --flip %s-merge.missnp'

    os.system('rm %s*.missnp' % tmp)
    os.system('rm %s*.bim' % tmp)
    os.system('rm %s*.bed' % tmp)
    os.system('rm %s*.fam' % tmp)

    #find flipped
    cmd1 = plink_merge_cmd % (in1, in2, tmp)
    os.system(cmd1)

    #flip alleles
    cmd2 = plink_flip_cmd % (in2, tmp, tmp)
    os.system(cmd2)

    #merge
    cmd3 = plink_merge_cmd % (in1, tmp, out)
    os.system(cmd3)

