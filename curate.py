import pandas as pd
import numpy as np

def read_fam_files(fam):
    f_content = [pd.read_table(f, header=None, sep="[\t ]+", engine='python') for f in fam]
    all_sample_ids =  pd.concat(f_content)[[0,1]]
    return all_sample_ids

def read_sample_files(sample_files):
    s_content = [pd.read_csv(f, engine='python') for f in sample_files]
    all_sample_ids =  pd.concat(s_content)
    null = pd.isnull(all_sample_ids.PID)
    all_sample_ids = all_sample_ids[np.logical_not(null)]
    all_sample_ids.PID = pd.Series(all_sample_ids.PID, dtype=int)
    return all_sample_ids
    

def read_location_data(fname):
    loc = pd.read_csv(fname)
    return loc

def read_dup_dict(d):
    return dict(np.loadtxt(d, dtype=int))
    

def deduplicate(data, dup_dict):
    new_ids =[]
    for pid in data.PID:
        pid = int(pid)
        if pid in dup_dict:
            new_ids.append(dup_dict[pid])
        else:
                new_ids.append(pid)
    data.PID = new_ids
    return data

def remove_sample_without_loc(data, loc):
    intersection = set.intersection(set(data.PID), set(loc.ID))
    data = data[data.PID.isin(intersection)]
    loc = loc[loc.ID.isin(intersection)]
    print(data.shape, loc.shape)
    return data, loc

try:
    fam_files = read_fam_files(snakemake.input.famfiles)
    sample_files = read_sample_files(snakemake.input.sample_files)
    location_data = read_location_data(snakemake.input.data)
    dup_dict = read_dup_dict(snakemake.input.dup_dict)
    deduplicate(sample_files, dup_dict)
    d, l = remove_sample_without_loc(sample_files, location_data)
    d.to_csv(snakemake.output.sample_out, index=None)
    l.to_csv(snakemake.output.loc_out, index=None)
except NameError:
    pass
