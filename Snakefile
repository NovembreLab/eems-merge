import pandas as pd
import numpy as np
import snakemake.utils
from snakemake.utils import R
import os

include: 'cleanup_snp.snake'
include: 'meta.snake'

base = lambda x : os.path.splitext(x)[0]


def plink_merge(in1, in2, out, tmp='tmp/merge',
        triallelic='merged/triallelic.txt'):
    """ uses plink to merge datasets

    this is hackish, but works: it runs plink three times to
        merge two datasets. 
        
        RUN1: naive merge. Will work if no SNP are flipped,
            otherwise the flip-file will be created.

        (if flipped SNP present)
        RUN2a : flip second input file
        RUN3a : merge in1 with flipped file

        (if no flipped SNP present)
        RUN2b : flip will fail, 
        RUN3b : merge in1 with flipped file, this will not 
                change stuff
    
        (if still not working)
        RUN4: as RUN2, but excluding triallelic stuff'
        RUN5: as RUN3
    """
    plink_merge_cmd = 'plink --bfile %s --bmerge %s --out %s --make-bed'
    plink_flip_cmd = 'plink --bfile %s --out %s --make-bed --flip %s-merge.missnp'
    plink_exclude_cmd = plink_flip_cmd + ' --exclude %s-merge.missnp'

    os.system('rm -f %s*.missnp' % tmp)
    os.system('rm -f %s*.bim' % tmp)
    os.system('rm -f %s*.bed' % tmp)
    os.system('rm -f %s*.fam' % tmp)

    #find flipped
    print('-------\n\n\n\n\n\n\n----- RUN 1')
    cmd1 = plink_merge_cmd % (in1, in2, tmp)
    os.system(cmd1)

    #flip alleles
    print('-------\n\n\n\n\n\n\n----- RUN 2')
    cmd2 = plink_flip_cmd % (in2, tmp, tmp)
    os.system(cmd2)

    #merge
    print('-------\n\n\n\n\n\n\n----- RUN 3')
    cmd3 = plink_merge_cmd % (in1, tmp, out)
    os.system(cmd3)
     
    if os.path.exists('%s-merge.missnp' % out):
        shell('cat %s-merge.missnp >> %s' % (out, triallelic))
        print('-------\n\n\n\n\n\n\n----- RUN 4')
        cmd4 = plink_exclude_cmd % (in2, tmp, tmp, out)
        os.system(cmd4)

        print('-------\n\n\n\n\n\n\n----- RUN 5')
        cmd3 = plink_merge_cmd % (in1, tmp, out)
        os.system(cmd3)
    
        os.remove('%s-merge.missnp' % out)
        

def unify_bim_id(data):
    data.columns = ['CHROM', 'ID', 'X', 'POS', 'ALLELE1', 'ALLELE2']
    s1 = data['CHROM']
    s2 = data['POS']
    s = zip(s1, s2)
    new_names = ["%s_%s" % i for i in s] 
    data['ID'] = new_names 
    return data

rule unique_estonian_meta:
    """ creates population raw file for estonian data 
        this file is one of the bases for regions/location_simplified.csv
    """
    input: "sources/Data_for_Ben_Meta.xlsx"
    output: "regions/estonia.csv"
    run:
        d = pd.read_excel(input[0])  
        d2 = (d[['Population',
                 'Group_Population',
                 'Region',
                 'State / Province / City',
                 'Country of origin / collected in',
                 'Continent',
                 'Sampling Location (Cambridgesamples)',
                 'Latitude',
                 'Longitude']])
        d2.drop_duplicates().to_csv(output[0], index=False)

rule niceify_lat_long:
    """sed 's/"[ ]*\([-0-9.]*\),[ ]*\([-0-9.]*\)",/\\1,\\2/' {input} >{output}"""
    input: "regions/location_simplified.csv"
    output: "regions/location_coords.csv"
    shell: " echo 3 "
 
rule estonian_studies:
    """ creates a file with all the unique study ID from
        the estonian data. currently not furhter used
    """
    input: "sources/Data_for_Ben_Meta.xlsx"
    output: "regions/estonian_studies.csv"
    run:
        d = pd.read_excel(input[0], sheetname='Sheet1')  
        d.source.drop_duplicates().to_csv(
            output[0], index=False)                                                  

rule uniqueify_hugo:
    input:
        "regions/location_coords.csv"
    output:
        "regions/location_hugo.csv"
    run:
        s = """d <- read.csv("{input[0]}")
            a <- d[d$Source == 'HUGO',]
            b <- aggregate(cbind(a$Latitude, a$Longitude),
                   list(a$Population), mean)
            names(b) <- c("Population", "Latitude", "Longitude")
            indices <- a[!duplicated(a$Population),c('ID', 'Population')]
            m <- merge(indices, b)[,c(2,1,3,4)]

            write.csv(m, "{output[0]}", row.names=F, quote=F)
            """
        print(s)
        R(s)

rule extract_hugo_loc:
    """
        from the HUGO html file, extracts the latitude and longitude coordinates
    """
    input: "sources/PASNP_Map.htm"
    output:
        temp('hugo_coords.txt'),
        temp('hugo_labels.txt'),
        "sources/hugo_meta.csv",
        temp('hugo_long.txt'),
        temp('hugo_lat.txt'),
    shell:
        """grep Latitude {input[0]}  | sed 's/.*\\[\\(.*\\)\].*/\\1/' > {output[0]}
        cut -f1 -d, {output[0]} > {output[3]}
        cut -f2 -d, {output[0]} > {output[4]}
        grep Latitude {input[0]}  | sed "s/.*map,'\\(.*\\)','.*/\\1/" |cut -f1 -d" "> {output[1]}
        paste {output[1]} {output[4]} {output[3]} -d, > {output[2]}"""

rule extract_stoneking_loc:
    input: "sources/Stoneking.pops.txt"
    output: "regions/Stoneking.pops.csv"
    run: 
        d = pd.read_table(input[0])
        d = d[['pop ID', 'sampling location']].drop_duplicates()
        d.to_csv(output[0], index=False)


rule download_xing:
    output:
        protected("raw/xing_sample_pop.txt")
    shell:
        "wget -O {output} http://jorde-lab.genetics.utah.edu/pub/affy6_xing2010/UID_344_Map.txt"


rule download_converge:
    output:
        protected("raw/converge_raw_data.txt.gz")
    shell:
        "wget -O {output} ftp://climb.genomics.cn/pub/10.5524/100001_101000/100155/MD_GWAS_SNPresults.txt.gz"

rule download_paschou:
    output:
        "data/paschou.zip"
    shell:
        "wget -O {output} "
        "http://drineas.org/Maritime_Route/RAW_DATA/PLINK_FILES/MARITIME_ROUTE.zip"
       

rule extract_paschou:
    input: "raw/paschou.zip"
    output:
        "raw/MARITIME_ROUTE.bed",
        "raw/MARITIME_ROUTE.bim",
        "raw/MARITIME_ROUTE.fam",
    shell: "unzip {input}"



rule deduplicate_master_table:
    input:
        table='regions/location_coords.csv',
        script='deduplicate_master_file.py'
    output:
        'regions/locations_deduplicated.csv',
        'duplicate_dict.txt'
    script:
        input.script

        
rule sample_pop_estonia:
    """ assigns estonian individuals id to population
        ids """
    input:
        "sources/Data_for_Ben_Meta.xlsx",
        'regions/location_coords.csv'
    output:
        'intermediate/sample_pop_estonians.csv'
    run:
        a = pd.read_excel(input[0])
        b = pd.read_csv(input[1])
        a = a.fillna('')
        b = b.fillna('')
        s1 = set(a.columns)
        s2 = set(b.columns)
        indices = s1.intersection(s2)
        indices.discard('Latitude')
        indices.discard('Longitude')
        #indices.discard('Country of origin / collected in')

        c = pd.merge(a, b, on=list(indices), how='left')
        d = c[['Sample ID', 'ID']]
        d.columns = 'SID', 'PID'
        d['Source'] = 'Estonians'
        d['Permission'] = 'Public'

        d.to_csv(output[0], index=False)

rule sample_pop_popres:
    input:
        "sources/POPRES_Phenotypes.txt"
    output:
        temp("intermediate/POPRES_SINGLETAB.csv"),
        'intermediate/sample_pop_popres.csv'
    shell:
        "sed  's/\t$//' {input} >{output[0]};"
        "python3 scripts/popres_keep4grandparents.py"

rule sample_pop_stoneking:
    """creates a table assigning individual labels to 
        population labels for the stoneking dataset"""
    input:
        "sources/Stoneking.pops.txt",
        "regions/location_coords.csv",
    output:
        'intermediate/sample_pop_stoneking.csv'
    script:
        "scripts/stoneking_merge.R"
    
rule sample_pop_paschou:
    input:
        loc='regions/location_coords.csv',
        samples='data/MARITIME_ROUTE.fam'
    output:
        'intermediate/sample_pop_paschou.csv'
    run:
        loc = pd.read_csv(input.loc)
        locp = loc[loc.Source=='Paschou2014']

        samp = pd.read_csv(input.samples, header=None,  sep=" ")
        samp = samp[[0,1]]
        samp.columns = ['Population', 'SampleID']

        easy_case = [(i.SampleID, i.Population, 
                locp[locp.Population==i.Population].ID.iloc[0]) 
            for  i in samp.itertuples() 
                if i.Population in list(locp.Population)]
        hard_case = [[i.SampleID, i.Population, 
                loc[loc.Population==i.Population].ID]
            for i in samp.itertuples()
                if i.Population not in list(locp.Population)]
        for i, h in enumerate(hard_case):
            if len(hard_case[i][2]) > 0:
                hard_case[i][2] = int(hard_case[i][2].iloc[0])
            else:
                if h[1] == 'Bedouin': hard_case[i][2] = 511
                if h[1] == 'Russia': hard_case[i][2] = 534
                if h[1] == 'Yemenite': hard_case[i][2] = 11

        cases = easy_case + hard_case
        cases = [(c[0], c[2]) for c in cases]
        cases = np.array(cases)
        src = np.repeat('PASCHOU', cases.shape[0])
        perm = np.repeat('Public', cases.shape[0])
        res = np.c_[cases, src, perm]

        np.savetxt(output[0], res, fmt="%s", delimiter=",",
                   header='SID,PID,Source,Permission', comments='')

rule sample_pop_hugo:
    input:
        loc='regions/location_hugo.csv',
        samples='data/hugo.fam'
    output:
        'intermediate/sample_pop_hugo.csv'
    run:
        loc = pd.read_csv(input.loc)
        samp = pd.read_csv(input.samples, header=None,  sep=" ")
        loc_dict = dict((l[1].Population, l[1].ID) for l in loc.iterrows())
        s = samp[0]
        res = [(s, loc_dict[s[:5]]) for s in samp[0] if s.count("-") == 4]
        res = np.array(res)
        src = np.repeat('HUGO', res.shape[0])
        perm = np.repeat('Private', res.shape[0])
        res = np.c_[res, src, perm]
        np.savetxt(output[0], res, fmt="%s", delimiter=",", 
                  header='SID,PID,Source,Permission', comments=''), 

rule sample_pop_verdu:
    input:
        loc='regions/location_coords.csv',
        samples='data/verdu.fam'
    output:
        'intermediate/sample_pop_verdu.csv'
    run:
        samp = pd.read_table(input.samples, header=None,  sep=" ")
        samp = samp[[0,1]]
        samp.columns = ['Population', 'SID']
        samp['PID'] = ''
        samp.PID.loc[samp.Population == 'TLI'] = 786   
        samp.PID.loc[samp.Population == 'HAI'] = 785   
        samp.PID.loc[samp.Population == 'NSI'] = 784   
        samp.PID.loc[samp.Population == 'TSI'] = 783   
        samp.PID.loc[samp.Population == 'STS'] = 787   
        samp.PID.loc[samp.Population == 'SPL'] = 788   
        samp['Source'] = 'VERDU'
        samp['Permission'] = 'Private'
        samp[['SID', 'PID', 'Source', 'Permission']].to_csv(output[0], index=False)

rule sample_pop_xing:
    input:
        loc='regions/location_simplified.csv',
        samples='raw/xing_sample_pop.txt'
    output:
        'intermediate/sample_pop_xing.csv'
    run:
        loc = pd.read_csv(input.loc)
        samp = pd.read_table(input.samples, header=None)
        loc = loc[['ID', 'Population', 'Source']]
        xing_loc = loc[loc['Source'] == 'Xing'] 
        samp = samp[[0, 2]]
        samp.columns = ['SampleID', 'Population']  
        m = pd.merge(samp, xing_loc, on=['Population'], how='left')
        m = m[['SampleID', 'ID']]
        m.columns = ['SID', 'PID']
        m['Source'] = 'XING'
        m['Permission'] = 'Public'
        m.to_csv(output[0], index=False)
        
rule sample_pop_lazaridis:
    input:
        loc='regions/location_simplified.csv',
        samples='/home/peterb/old_data/data_clean/raw/EuropeAllData/vdata.ind'
    output:
        'intermediate/sample_pop_lazaridis.csv'
    run:
        loc = pd.read_csv(input.loc)
        samp = pd.read_csv(input.samples, header=None, sep="[ ]+",
            engine='python')
        loc = loc[['ID', 'Group_Population', 'Source']]
        ho_loc = loc[loc['Source'] == 'Human Origins'] 
        samp = samp[[0, 2]]
        samp.columns = ['SampleID', 'Group_Population']  
        m = pd.merge(samp, ho_loc, on=['Group_Population'], how='left')
        m = m[['SampleID', 'ID']]
        m[np.logical_not(np.isnan(m.ID))]
        m.columns = ['SID', 'PID']
        m['Source'] = 'LAZARIDIS'
        m['Permission'] = 'Private_Demo'
        m.to_csv(output[0], index=False)

rule make_curated_table:
    input:
        data='regions/locations_deduplicated.csv',
        dup_dict='duplicate_dict.txt',
        famfiles=['data/Data_for_Ben.fam',
                  'data/hugo.fam',
                  'data/MARITIME_ROUTE.fam',
                  'data/POPRES_Genotypes_QC1_v2.fam',
                  'data/reich2011.fam',
                  'data/vdata.fam',
                  'data/verdu.fam',
                  'data/xing.fam'],
        sample_files=['intermediate/sample_pop_estonians.csv',
                      'intermediate/sample_pop_hugo.csv',
                      'intermediate/sample_pop_lazaridis.csv',
                      'intermediate/sample_pop_paschou.csv',
                      'intermediate/sample_pop_popres.csv',
                      'intermediate/sample_pop_stoneking.csv',
                      'intermediate/sample_pop_verdu.csv',
                      'intermediate/sample_pop_xing.csv'],
        script='curate.py'
    output:
        sample_out='intermediate/sample_pop_all.csv',
        loc_out='intermediate/locations_all.csv'
    script:
        input.script

rule make_individual_list_to_keep:
    input:
        samples='intermediate/sample_pop_all.csv',
    output:
        'plink/indiv_list.txt'
    run:
        f = pd.read_csv(input.samples)
        n = ["%s%s" % (s.Source, s.SID) for s in f.itertuples()]
        with open(output[0], 'w') as fs:
            for i in n:
                fs.write('%s %s\n' % (i, i))
    


names_fam_in=['Data_for_Ben', 'hugo', 'MARITIME_ROUTE', 
    'POPRES_Genotypes_QC1_v2',
    'reich2011',
    'vdata',
    'verdu',
    'xing']
names_fam_out=['Estonians', 'HUGO', 'PASCHOU', 
    'POPRES',
    'Stoneking',
    'LAZARIDIS',
    'VERDU',
    'XING']
name_dict=dict(zip(names_fam_out, names_fam_in))
def rename_fun(wildcards):
    m = name_dict[wildcards.n]
    d  = dict()
    return ['data/%s.bed' %m, 'data/%s.bim' %m, 'data/%s.fam' %m]
    

rule rename_fam_files:
    """rename individuals 1. to ensure identifiers are unique and
        2. to ensure the two ids are the same
    """
    input:
        fam=expand('data2/{n}.fam', n=names_fam_out),
        bed=expand('data2/{n}.bed', n=names_fam_out),
        bim=expand('data2/{n}.bim', n=names_fam_out)

rule rename_single_fam_file:
    input:
        rename_fun,
        indiv_list='plink/indiv_list.txt'
    output:
        tmpbim=temp('tmp/tmp{n}.bim'),
        tmpfam=temp('tmp/tmp{n}.fam'),
        invalid='data2/{n}.invalid',
        fam='data2/{n}.fam',
        bed='data2/{n}.bed',
        bim='data2/{n}.bim'
    run:
        s='awk \'BEGIN{{OFS=" "}}; '
        s+= ' {{ print "{wildcards.n}"$2, "{wildcards.n}"$2, $3, $4, $5, $6 }} \''
        s+= ' {input[2]} > {output.tmpfam}'
        shell(s)
        #shell('ln -sf  `pwd`/{input[1]} `pwd`/{output.bed}')

        data = pd.read_table(input[1], sep="[ \t\s]+", header=None,
            engine='python')
        data[0] = data[0].astype(int)
        data = unify_bim_id(data)
        data.to_csv(output.tmpbim, index=False,
            header=None, sep=" ")

        print(data.columns)

        valid_alleles =  np.logical_and(data.ALLELE1.isin(list('ACGTacgt')),
                                        data.ALLELE1.isin(list('ACGTacgt')))
        ambiguous1 = np.logical_and(data.ALLELE1.isin(list('Aa')),
                                   data.ALLELE2.isin(list('Tt')))
        ambiguous2 = np.logical_and(data.ALLELE1.isin(list('Tt')),
                                   data.ALLELE2.isin(list('Aa')))
        ambiguous3 = np.logical_and(data.ALLELE1.isin(list('Cc')),
                                   data.ALLELE2.isin(list('Gg')))
        ambiguous4 = np.logical_and(data.ALLELE1.isin(list('Gg')),
                                   data.ALLELE2.isin(list('Cc')))
        ambiguous12 = np.logical_or(ambiguous1, ambiguous2)
        ambiguous34 = np.logical_or(ambiguous3, ambiguous4)
        ambiguous = np.logical_or(ambiguous12, ambiguous34)
        valid = np.logical_and(valid_alleles,
                               np.logical_not(ambiguous))
        invalid_ids = data.ID[np.logical_not(valid)]
        invalid_ids.to_csv(output.invalid, index=False,
            header=None, sep=" ")

        shell('cp {output.tmpbim} tmp/test')

        outname = os.path.splitext(output.fam)[0]
        s = 'plink --bed {input[0]} --bim {output.tmpbim} --fam {output.tmpfam} '
        s += '--keep {input.indiv_list} --make-bed --out %s ' % outname
        s += '--autosome --exclude {output.invalid}'
        shell(s)
        

rule merge_everything:
    input:
        expand('data2/{name}.{ext}', 
            name=names_fam_out,
            ext=['bed', 'bim', 'fam']),
        invalid=expand('data2/{name}.{ext}', 
            name=names_fam_out,
            ext=['invalid'])
    output:
        bed='merged/merge_master.bed',
        bim='merged/merge_master.bim',
        fam='merged/merge_master.fam',
        bed2='merged/nohugo_master.bed',
        bim2='merged/nohugo_master.bim',
        fam2='merged/nohugo_master.fam',
        triallelic='merged/triallelic.txt',
        invalid='merged/invalid.txt',
        triallelic2='merged/triallelic2.txt',
        invalid2='merged/invalid2.txt'
    run:
        out = os.path.splitext(output.bed)[0]
        tmp = out + '_tmp%s'

        shell('rm -f {output.triallelic} && touch {output.triallelic}')

        plink_merge('data2/Estonians', 'data2/LAZARIDIS', tmp % 1)
        plink_merge(tmp % 1, 'data2/HUGO', tmp % 2)
        plink_merge(tmp % 2, 'data2/POPRES', tmp % 3)
        plink_merge(tmp % 3, 'data2/PASCHOU', tmp % 4)
        plink_merge(tmp % 4, 'data2/Stoneking', tmp % 5)
        plink_merge(tmp % 5, 'data2/VERDU', tmp % 6)
        plink_merge(tmp % 6, 'data2/XING', tmp % 7)

        shell('cp {output.triallelic} {output.invalid}')
        for invalid in input.invalid:
            shell('cat {invalid} >> {output.invalid}')

        s = 'plink --bfile %s --make-bed --out %s --exclude %s'
        s = s % (tmp % 7, out, output.invalid)
        shell(s)

        ###############################################
        # same w/o hugo
        out = os.path.splitext(output.bed2)[0]
        tmp = out + '_tmp%s'

        shell('rm -f {output.triallelic2} && touch {output.triallelic2}')

        plink_merge('data2/Estonians', 'data2/LAZARIDIS', tmp % 1)
        plink_merge(tmp % 1, 'data2/POPRES', tmp % 3)
        plink_merge(tmp % 3, 'data2/PASCHOU', tmp % 4)
        plink_merge(tmp % 4, 'data2/Stoneking', tmp % 5)
        plink_merge(tmp % 5, 'data2/VERDU', tmp % 6)
        plink_merge(tmp % 6, 'data2/XING', tmp % 7)

        shell('cp {output.triallelic2} {output.invalid2}')
        for invalid in input.invalid:
            shell('cat {invalid} >> {output.invalid2}')

        s = 'plink --bfile %s --make-bed --out %s --exclude %s'
        s = s % (tmp % 7, out, output.invalid2)
        shell(s)

       

rule prune_mat:
    input:
        bed='{name}_master.bed',
        bim='{name}_master.bim',
        fam='{name}_master.fam',
    output:
        pruned_in='{name}_master.prune.in'
    run:
        name = wildcards.name + '_master'
        s = 'plink --bfile %s --indep-pairwise 1000 1000 .1 --out %s'
        shell(s % (name, name))

rule make_rel_mat:
    input:
        bed='{name}_master.bed',
        bim='{name}_master.bim',
        fam='{name}_master.fam',
        pruned_in='{name}_master.prune.in'
    output:
        grm='{name}_master.grm.bin',
    run:
        name = wildcards.name + '_master'
        s = 'plink --bfile {name} --make-grm-bin'
        s += ' --extract {input.pruned_in} ' 
        s += '--out {name}'
        shell(s)

rule rel_filter:
    input:
        bed='{name}_master.bed',
        bim='{name}_master.bim',
        fam='{name}_master.fam',
        pruned_in='{name}_master.prune.in',
        grm='{name}_master.grm.bin',
    output:
        bed='{name}_prune{cutoff}.bed',
        bim='{name}_prune{cutoff}.bim',
        fam='{name}_prune{cutoff}.fam',
    run:
        name = base(input.bed)
        outname = base(output.bed)
        grm = base(base(input.grm))
        cutoff = float(wildcards.cutoff)
        s='plink  --rel-cutoff {cutoff} '
        s += ' --grm-bin {grm} '
        s += ' --out {outname} '
        shell(s)
        s='plink  --bfile {name} --keep {outname}.grm.id '
        s += ' --out {outname} --make-bed '
        shell(s)
        
rule merge_tib_meta:
    input:
        ail="pgs/gvar.indiv_label",
        aip="pgs/gvar.indiv_prov",
        apd="pgs/gvar.pop_display",
        apg="pgs/gvar.pop_geo",
        bil="tib/tibetan.indiv_label",
        bip="tib/tibetan.indiv_prov",
        bpd="tib/tibetan.pop_display",
        bpg="tib/tibetan.pop_geo",
        fam="merged/tib_prune0.6.fam"
    output:
        oil="pgs/gvar2.indiv_label",
        oip="pgs/gvar2.indiv_prov",
        opd="pgs/gvar2.pop_display",
        opg="pgs/gvar2.pop_geo",
    run:
        import pandas as pd
        fam = pd.read_table(input.fam, sep=" ", names=range(6))
        fam = pd.DataFrame({'sampleId':fam[0]})

        p1, p2 = pd.read_csv(input.ail), pd.read_csv(input.bil)
        iil = p1.append(p2)
        assert len(p1.popId.unique()) + len(p2.popId.unique()) ==\
            len(iil.popId.unique())
        iil = iil.merge(fam)
        assert sum(iil.sampleId.duplicated()) == 0
        iil.to_csv(output.oil, index=False, columns=('sampleId', 'popId'))

        p1, p2 = pd.read_csv(input.aip), pd.read_csv(input.bip)
        iip = p1.append(p2)
        iip = iip.merge(fam)
        assert set(iil.sampleId) == set(iip.sampleId)
        assert sum(iip.sampleId.duplicated()) == 0
        iip.to_csv(output.oip, index=False)

        p1, p2 = pd.read_csv(input.apd), pd.read_csv(input.bpd)
        ipd = p1.append(p2)
        assert sum(ipd.popId.duplicated()) == 0
        ipd.to_csv(output.opd, index=False)

        p1, p2 = pd.read_csv(input.apg), pd.read_csv(input.bpg)
        ipg = p1.append(p2)
        assert sum(ipg.popId.duplicated()) == 0
        assert set(ipg.popId) == set(ipd.popId)
        ipg.to_csv(output.opg, index=False)

rule merge_qatar_to_tib_meta:
    input:
        ail="pgs/gvar2.indiv_label",
        aip="pgs/gvar2.indiv_prov",
        apd="pgs/gvar2.pop_display",
        apg="pgs/gvar2.pop_geo",
        bil="qatari/qatari.indiv_label",
        bip="qatari/qatari.indiv_prov",
        bpd="qatari/qatari.pop_display",
        bpg="qatari/qatari.pop_geo",
        fam="merged/tq_prune0.6.fam",
        fix=["pgs/gvar3.names", "pgs/update_pos.csv", "pgs/merge.csv",
            "fix.R"
        ]
    output:
        oil="pgs/gvar3.indiv_label",
        oip="pgs/gvar3.indiv_prov",
        opd="pgs/gvar3.pop_display",
        opg="pgs/gvar3.pop_geo",
    run:
        import pandas as pd
        fam = pd.read_table(input.fam, sep=" ", names=range(6))
        fam = pd.DataFrame({'sampleId':fam[0]})

        p1, p2 = pd.read_csv(input.ail), pd.read_csv(input.bil)
        print(p2.head())
        iil = p1.append(p2)
        assert len(p1.popId.unique()) + len(p2.popId.unique()) ==\
            len(iil.popId.unique())
        iil = iil.merge(fam)
        assert sum(iil.sampleId.duplicated()) == 0
        iil.to_csv(output.oil, index=False, columns=('sampleId', 'popId'))

        p1, p2 = pd.read_csv(input.aip), pd.read_csv(input.bip)
        iip = p1.append(p2)
        iip = iip.merge(fam)
        s1, s2 =set(iil.sampleId), set(iip.sampleId)
        assert set(iil.sampleId) == set(iip.sampleId)
        assert sum(iip.sampleId.duplicated()) == 0
        iip.to_csv(output.oip, index=False)

        p1, p2 = pd.read_csv(input.apd), pd.read_csv(input.bpd)
        ipd = p1.append(p2)
        assert sum(ipd.popId.duplicated()) == 0
        ipd.to_csv(output.opd, index=False)

        p1, p2 = pd.read_csv(input.apg), pd.read_csv(input.bpg)
        ipg = p1.append(p2)
        assert sum(ipg.popId.duplicated()) == 0
        assert set(ipg.popId) == set(ipd.popId)
        ipg.to_csv(output.opg, index=False)

        R("source('fix.R')")
	


rule all:
    input:
        expand('merged/{name}_prune{cf}.bed',
            name=['tib', 'merge', 'nohugo'],
            cf=['0.75', '0.6', '0.25'])


#merging in tib
rule filter_tib:
    input:
        "tib/HGDP_Tibetan_Merged_160509.pop_geo",
        "tib/HGDP_Tibetan_Merged_160509.indiv_label",
        "tib/HGDP_Tibetan_Merged_160509.indiv_prov",
        "tib/HGDP_Tibetan_Merged_160509.pop_display",
        loc='intermediate/locations_all.csv',
    output:
        "tib/tibetan.csv",
        "tib/tib.plink",
        "tib/tibetan.indiv_label",
        "tib/tibetan.indiv_prov",
        "tib/tibetan.pop_geo",
        "tib/tibetan.pop_display",
    script: "tib/loc.R"


rule get_hgdp_dudes:
    input:
        bed="tib/HGDP_Tibetan_Merged_160509.bed",
        bim="tib/HGDP_Tibetan_Merged_160509.bim",
        fam="tib/HGDP_Tibetan_Merged_160509.fam",
        samples="tib/tib.plink",
    output:
        bed="tib/eatw_tibetans.bed",
        bim="tib/eatw_tibetans.bim",
        fam="tib/eatw_tibetans.fam",
        tmp=temp("tib/other.fam")
    run:
        inname=base(input.bed)
        outname=base(output.bed)
        s = "plink --bfile {inname} --make-bed --out {outname}"
        s += " --keep {input.samples}"
        shell(s)
        with open(output.fam, 'r') as f:
            with open(output.tmp, "w") as o:
                for line in f:
                    l = line.split()
                    l[0] = l[1]
                    o.write(" ".join(l) + "\n")
        shell("mv {output.tmp} {output.fam}")

        with open(output.bim, 'r') as f:
            with open(output.tmp, "w") as o:
                for line in f:
                    l = line.split()
                    l[1] = "%s_%s" % (l[0],l[3])
                    o.write(" ".join(l) + "\n")
        shell("mv {output.tmp} {output.bim}")

rule add_tibetans:
    input:
        bed='merged/merge_master.bed',
        bim='merged/merge_master.bim',
        fam='merged/merge_master.fam',
        tibbed='tib/eatw_tibetans.bed',
        tibbim='tib/eatw_tibetans.bim',
        tibfam='tib/eatw_tibetans.fam',
    output:
        bed='merged/tib_master.bed',
        bim='merged/tib_master.bim',
        fam='merged/tib_master.fam',
    run:
        in1 = os.path.splitext(input.bed)[0]
        in2 = os.path.splitext(input.tibbed)[0]
        out = os.path.splitext(output.bed)[0]
        plink_merge(in1, in2, out)

rule flip_qatari:
    input:
        bed="qatari/NWAfrica_HM3_Qat.bed",
        bim="qatari/NWAfrica_HM3_Qat.bim",
        fam="qatari/NWAfrica_HM3_Qat.fam",
    output:
        bed="qatari/hg37.bed",
        bim="qatari/hg37.bim",
        fam="qatari/hg37.fam",
    script: 
        "scripts/flip_henn2012.py"


rule clean_qatari:
    input:
        bed="qatari/hg37.bed",
        bim="qatari/hg37.bim",
        fam="qatari/hg37.fam",
        samples="qatari/qatari.indiv_prov"
    output:
        bed="qatari/qatari.bed",
        bim="qatari/qatari.bim",
        fam="qatari/qatari.fam",
        tmp=temp("qatari/other.fam"),
        tmp_id=temp("qatari/ids.txt"),
    run:
        s = """cut -d, -f1 {input.samples} | sed 's/"//g' > {output.tmp} """
        shell(s)
        s = 'grep -F -f {output.tmp} {input.fam} | cut -f-2 -d" "' 
        s += " > {output.tmp_id}"
        shell(s)

        inname=base(input.bed)
        outname=base(output.bed)
        s = "plink --bfile {inname} --make-bed --out {outname}"
        s += " --keep {output.tmp_id}"
        shell(s)
        with open(output.fam, 'r') as f:
            with open(output.tmp, "w") as o:
                for line in f:
                    l = line.split()
                    l[1] = "QATAR" + l[1]
                    l[0] = l[1]
                    o.write(" ".join(l) + "\n")
        shell("mv {output.tmp} {output.fam}")

        with open(output.bim, 'r') as f:
            with open(output.tmp, "w") as o:
                for line in f:
                    l = line.split()
                    l[1] = "%s_%s" % (l[0],l[3])
                    o.write(" ".join(l) + "\n")
        shell("mv {output.tmp} {output.bim}")

rule add_qatari:
    input:
        bed='merged/tib_master.bed',
        bim='merged/tib_master.bim',
        fam='merged/tib_master.fam',
        tibbed='qatari/qatari.bed',
        tibbim='qatari/qatari.bim',
        tibfam='qatari/qatari.fam',
    output:
        bed='merged/tq_master.bed',
        bim='merged/tq_master.bim',
        fam='merged/tq_master.fam',
    run:
        in1 = os.path.splitext(input.bed)[0]
        in2 = os.path.splitext(input.tibbed)[0]
        out = os.path.splitext(output.bed)[0]
        plink_merge(in1, in2, out)
