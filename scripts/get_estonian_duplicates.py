import pandas as pd

#bim = snakemake.input.bim
inbim = 'data/Data_for_Ben.bim'
bim = pd.read_table(inbim, header=None)
dup_ids = bim[bim[3].duplicated()][1]

