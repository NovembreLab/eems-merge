import pandas as pd
inds = pd.read_table("meta/individuals2.txt", sep=r"\s+")
print inds.shape
inds.drop_duplicates(subset=['ID'], take_last=True, inplace=True)
print inds.shape
inds.to_csv("meta/individuals3.txt", sep="\t", index=False)
