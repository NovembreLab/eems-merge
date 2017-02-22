#!/usr/bin/env python

import os, sys

convertf = "/home/peterb/old_data/programs/EIG6.0.1/bin/convertf"

s="""
genotypename:    %s.geno
snpname:         %s.snp
indivname:       %s.ind
outputformat:    PACKEDPED
genotypeoutname: %s.bed
snpoutname:      %s.bim
indivoutname:    %s.fam
"""

s2="""
genotypename:    %s.geno
snpname:         %s.snp
indivname:       %s.ind
outputformat:    PED
genotypeoutname: %s.ped
snpoutname:      %s.bim
indivoutname:    %s.fam
"""


def run(files_to_convert):


    for f in files_to_convert:
        with open( "tmp.par", "w") as handle:
            handle.write( s%( (f,)*6))

    os.system("%s -p tmp.par" % convertf)
    os.system("sed -i 's/X/23/; s/Y/24/' %s.bim" % f)
#    os.system("rm %s.snp %s.ind %s.geno"%(f,f,f))

try:
    run([snakemake.wildcards.name])
except NameError:
    if len(sys.argv) > 1:
        run(sys.argv[1:])

