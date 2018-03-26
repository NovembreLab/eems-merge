import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from collections import Counter

CUTOFF = .3 #the cutoff for which samples to be assumed to be close/ same
            #location

fname = 'regions/location_coords.csv'

lc = pd.read_csv(fname)
lc.loc[pd.isnull(lc.Group_Population), 'Group_Population'] = ''
lc.not_most_accurate.loc[np.isnan(lc.not_most_accurate)] = 0
print(lc.shape)


exclude_jews = True
if exclude_jews:
    jews = [l.ID for l in lc.itertuples()  
            if 'jew' in str(l.Group_Population).lower() or
              'jew' in str(l.Population).lower()]                                                        

    print("JEWS: ", len(jews))

    lc  =  lc[np.logical_not(lc.ID.isin(jews))]
    print(lc.shape)



ident_cols = ['Population', 'Group_Population', 'Region'] 
ident_cols += ['State / Province / City', 'Country of origin / collected in']
ident_cols += ['Sampling Location (Cambridgesamples)']
simple_duplicate_data = lc[ident_cols]


dup_rows = np.where(simple_duplicate_data.duplicated(keep=False))[0]

G = lc.iloc[dup_rows].groupby('Population')

#G now is an object that contains info on all populations with multiple entries
# - a large fraction of those is from HUGO, which I deal with seperately.
# - several are countries, where we wand to merge them, and thake the mean

coord_means = G[['Latitude', 'Longitude']].mean()
first_ids = G.ID.first()


duplicate_dict = dict()
G_dict = dict(list(G))
for pop, pop_data in dict(list(G)).items():
    is_first = pop_data.ID.isin(first_ids)
    first = pop_data[is_first].ID
    is_not_first = np.logical_not(is_first)
    not_first = pop_data[is_not_first].ID
    for nf in not_first:
        duplicate_dict[nf] = first

lc = lc[np.logical_not(lc.ID.isin(duplicate_dict.keys()))]
deduplicated1 = pd.concat([first_ids, coord_means], axis=1)
for row_id, row in deduplicated1.iterrows():
    lc.loc[lc.ID ==  row.ID, 'Latitude'] = row.Latitude
    lc.loc[lc.ID ==  row.ID, 'Longitude'] = row.Longitude
print(lc.shape)


coords = lc[['Latitude', 'Longitude']]
dmat = squareform(pdist(coords, 'euclidean')) 
np.fill_diagonal(dmat, np.nan)
close_pairs = np.array(np.where(dmat<CUTOFF)) 
close_pairs =np.array([cp for cp in close_pairs.transpose()
                        if cp[0] > cp[1]])


close_pairs_ID = [lc.iloc[cp].ID for cp in close_pairs]   
cp_counter = Counter(np.array(close_pairs_ID).flatten())           


dup_coords = lc[coords.duplicated(keep=False)].groupby(['Latitude',
                                                            'Longitude'])
first_ids = dup_coords.ID.first()
duplicate_dict2 = dict()
G_dict = dict(list(dup_coords))
for pop, pop_data in G_dict.items():
    is_first = pop_data.ID.isin(first_ids)
    first = pop_data[is_first].ID
    is_not_first = np.logical_not(is_first)
    not_first = pop_data[is_not_first].ID
    for nf in not_first:
        duplicate_dict2[nf] = first
lc = lc[np.logical_not(lc.ID.isin(duplicate_dict2.keys()))]

print(lc.shape)

pgp = ['Population', 'Group_Population'] 




lc['Unique'] = lc.Population.copy()
lc.Unique.loc[lc.Unique.duplicated(keep=False)] = 'WTFNONE'
non_unique = lc.Unique == 'WTFNONE'
print(sum(lc.Unique == 'WTFNONE'))
lc.Unique.loc[non_unique] = lc.Group_Population[non_unique]
lc.Unique.loc[lc.Unique.duplicated(keep=False)] = 'WTFNONE'
print(sum(lc.Unique == 'WTFNONE'))


duplicate_dict3 = dict()
#fix altaians
lc.Unique.loc[lc.ID == 155] = 'Altaians_General'
lc.not_most_accurate.loc[lc.ID == 155] = 1
lc.Unique.loc[lc.ID == 156] = 'Altaians_Turochak'
lc.Unique.loc[lc.ID == 157] = 'Altaians_Telengit'

#fix Azeris
lc.Unique.loc[lc.ID == 131] = 'Azeris_Azerbaijan'
lc.Unique.loc[lc.ID == 54] = 'Azeris_Iran'

#fix Basque
duplicate_dict3[812] = 494
lc = lc[lc.ID != 812]

#fix Croats
lc.Unique.loc[lc.ID == 71] = 'Croats_Bosnia'
lc.Unique.loc[lc.ID == 69] = 'Croats_Croatia'


#fix Ethiopians
lc.Unique.loc[lc.ID == 1] = 'Ethopians_Amara'
lc.Unique.loc[lc.ID == 2] = 'Ethopians_Oromo'
lc.Unique.loc[lc.ID == 3] = 'Ethopians_Tigrani'


#fix Evenkis
lc.Unique.loc[lc.ID == 249] = 'Evenkis_General'
lc.not_most_accurate.loc[lc.ID == 249] = 1
lc.Unique.loc[lc.ID == 250] = 'Evenkis_Kuyumba'
lc.Unique.loc[lc.ID == 251] = 'Evenkis_Irokan'
lc.Unique.loc[lc.ID == 252] = 'Evenkis_Strelka'
lc.Unique.loc[lc.ID == 254] = 'Evenkis_Chirinda'
lc.Unique.loc[lc.ID == 255] = 'Evenkis_Tutonchany'
lc.Unique.loc[lc.ID == 256] = 'Evenkis_Tura'
lc.Unique.loc[lc.ID == 257] = 'Evenkis_Ekonda'
lc.Unique.loc[lc.ID == 258] = 'Evenkis_Kislokan'
lc.Unique.loc[lc.ID == 259] = 'Evenkis_Erbogachen'
lc.Unique.loc[lc.ID == 260] = 'Evenkis_Nidym'
#lc.Unique.loc[lc.ID == 261] = 'Evenkis_Kuyumba'
duplicate_dict3[261] = 250
lc = lc[lc.ID != 261]
lc.Unique.loc[lc.ID == 262] = 'Evenkis_Surinda'
lc.Unique.loc[lc.ID == 263] = 'Evenkis_Poligus'
lc.Unique.loc[lc.ID == 264] = 'Evenkis_Nakanno'
lc.Unique.loc[lc.ID == 265] = 'Evenkis_Baykit'

#fix Evens
lc.Unique.loc[lc.ID == 267] = 'Evens_General'
lc.not_most_accurate.loc[lc.ID == 267] = 1
lc.Unique.loc[lc.ID == 268] = 'Evens_Gadlya'
lc.Unique.loc[lc.ID == 269] = 'Evens_Gizhiga'
lc.Unique.loc[lc.ID == 270] = 'Evens_Evensk'
lc.Unique.loc[lc.ID == 271] = 'Evens_Paren'
lc.Unique.loc[lc.ID == 272] = 'Evens_Arman'
lc.Unique.loc[lc.ID == 273] = 'Evens_Topolovka'
lc.Unique.loc[lc.ID == 275] = 'Evens_Garmanda'

#fix Greeks
lc.Unique.loc[lc.ID == 63] = 'Greeks_General'
lc.not_most_accurate.loc[lc.ID == 63] = 1
lc.Unique.loc[lc.ID == 64] = 'Greeks_Thessaly'

#fix Japanese
duplicate_dict3[617] = 445
lc = lc[lc.ID != 617]

#fix Kazakhs
lc.Unique.loc[lc.ID == 152] = 'Kazakhs_CentralWest'
lc.Unique.loc[lc.ID == 153] = 'Kazakhs_General'

#fix Koryaks
lc.Unique.loc[lc.ID == 277] = 'Koryaks_General'
lc.not_most_accurate.loc[lc.ID == 277] = 1
lc.Unique.loc[lc.ID == 278] = 'Koryaks_Paren'
lc.Unique.loc[lc.ID == 279] = 'Koryaks_Evensk'
lc.Unique.loc[lc.ID == 280] = 'Koryaks_Topolovka'
lc.Unique.loc[lc.ID == 281] = 'Koryaks_Gizhiga'
lc.Unique.loc[lc.ID == 282] = 'Koryaks_NorthEast'
lc.not_most_accurate.loc[lc.ID == 282] = 1

#fix Kurds
lc.Unique.loc[lc.ID == 615] = 'Kurds_Iraq'
lc.Unique.loc[lc.ID == 52] = 'Kurds_Kazakhstan'
lc.not_most_accurate.loc[lc.ID == 52] = 1

#fix Kyrg
lc.Unique.loc[lc.ID == 146] = 'Kyrgyzians_TienShan'
lc.Unique.loc[lc.ID == 148] = 'Kyrgyzians_Murghab'
lc.Unique.loc[lc.ID == 149] = 'Kyrgyzians_Alichur'
lc.Unique.loc[lc.ID == 150] = 'Kyrgyzians_General'
lc.not_most_accurate.loc[lc.ID == 150] = 1

#fix Luhya
duplicate_dict3[621] = 339
lc = lc[lc.ID != 621]

#fix Mongolians
lc.Unique.loc[lc.ID == 220] = 'Mongolians_Halha'
lc.Unique.loc[lc.ID == 221] = 'Mongolians_General'

#fix Pamiris
lc.Unique.loc[lc.ID == 136] = 'Pamiris_Ishkashim'
lc.Unique.loc[lc.ID == 137] = 'Pamiris_Rushan'
lc.Unique.loc[lc.ID == 138] = 'Pamiris_Vanch'
lc.Unique.loc[lc.ID == 139] = 'Pamiris_Shugnan'

#fix Russians
lc.Unique.loc[lc.ID == 94] = 'Russians_NorthRussia'
lc.not_most_accurate.loc[lc.ID == 94] = 1
lc.Unique.loc[lc.ID == 93] = 'Russians_Arhangelsk'
lc.Unique.loc[lc.ID == 87] = 'Russians_Voronez'
lc.Unique.loc[lc.ID == 88] = 'Russians_Kursk'
lc.Unique.loc[lc.ID == 89] = 'Russians_Orjol'
lc.Unique.loc[lc.ID == 91] = 'Russians_Central'
lc.not_most_accurate.loc[lc.ID == 91] = 1
lc.Unique.loc[lc.ID == 90] = 'Russians_Smolensk'
lc.Unique.loc[lc.ID == 92] = 'Russians_Tver'
lc.Unique.loc[lc.ID == 95] = 'Russians_Kostroma'
lc.Unique.loc[lc.ID == 160] = 'Russians_Turochak'

#fix  Serbians
lc.Unique.loc[lc.ID == 73] = 'Serbians_Bosnia'
lc.Unique.loc[lc.ID == 72] = 'Serbians_Serbia'

#fix Shors
lc.Unique.loc[lc.ID == 230] = 'Shors_Gorge_Shoria'
lc.Unique.loc[lc.ID == 234] = 'Shors_Tashtagol'
lc.Unique.loc[lc.ID == 235] = 'Shors_Orton'

#fix Italians
lc.Unique.loc[lc.ID == 57] = 'Sicilians_West'
lc.Unique.loc[lc.ID == 58] = 'Sicilians_Central'
lc.Unique.loc[lc.ID == 59] = 'Sicilians_East'
lc.Unique.loc[lc.ID == 60] = 'Sicilians_South'

#fix Tajiks
lc.Unique.loc[lc.ID == 135] = 'Tajiks_Plain'
lc.not_most_accurate.loc[lc.ID == 135] = 1
lc.Unique.loc[lc.ID == 140] = 'Tajiks_General'
lc.not_most_accurate.loc[lc.ID == 140] = 1

#fix Tlingit
lc.Unique.loc[lc.ID == 786] = 'Tlingit_USA'
lc.not_most_accurate.loc[lc.ID == 786] = 1

#fix Tuvinians'
lc.Unique.loc[lc.ID == 223] = 'Tuvinians_General'
lc.Unique.loc[lc.ID == 224] = 'Tuvinians_Upsunur'

#fix Ukraine
lc.Unique.loc[lc.ID == 81] = 'Ukranians_Svetlovodsk'
lc.Unique.loc[lc.ID == 82] = 'Ukranians_Kharkov'
lc.Unique.loc[lc.ID == 83] = 'Ukranians_Poltava'
lc.Unique.loc[lc.ID == 84] = 'Ukranians_Belgorod'
lc.Unique.loc[lc.ID == 85] = 'Ukranians_Lviv'

#fix Yakuts
lc.Unique.loc[lc.ID == 238] = 'Yakuts_General'
lc.not_most_accurate.loc[lc.ID == 238] = 1
lc.Unique.loc[lc.ID == 239] = 'Yakuts_Essey'
lc.Unique.loc[lc.ID == 240] = 'Yakuts_Seymchan'
lc.Unique.loc[lc.ID == 241] = 'Yakuts_Vilyuy'
lc.Unique.loc[lc.ID == 242] = 'Yakuts_Yakutsk'
lc.Unique.loc[lc.ID == 243] = 'Yakuts_Ekonda'
lc.Unique.loc[lc.ID == 244] = 'Yakuts_Chirinda'
lc.Unique.loc[lc.ID == 245] = 'Yakuts_Olenek'
lc.Unique.loc[lc.ID == 246] = 'Yakuts_Ust-Alan'
lc.Unique.loc[lc.ID == 247] = 'Yakuts_Nyurba'
lc.Unique.loc[lc.ID == 248] = 'Yakuts_Nakkano'


#this datastructure is to ensure stuff is unique
G = lc[lc.Unique == 'WTFNONE'].groupby(pgp)
DG = dict(list(G))


#data structure to make sure populations make sense
country = 'Country of origin / collected in'                              
countries = np.unique(lc[country][lc[country].duplicated(keep=False)])    


#fix Bulgaria
duplicate_dict3[594] = 74

#fix Cambodia
duplicate_dict3[619] = 429

#fix China
duplicate_dict3[432] = 608
duplicate_dict3[717] = 608
lc.not_most_accurate.loc[lc.ID == 609] = 1

#fix Croatia
duplicate_dict3[585] = 69

#fix Egypt
duplicate_dict3[7] = 331
duplicate_dict3[797] = 331

#fix Ethiopia
duplicate_dict3[3] = 335

#fix France
duplicate_dict3[570] = 495

#fix Greece
duplicate_dict3[500] = 586

#fix Hungary
duplicate_dict3[577] = 75

#fix Indonesia
duplicate_dict3[721] = 758

#fix Italy
duplicate_dict3[808] = 58
duplicate_dict3[811] = 517
lc.not_most_accurate.loc[lc.ID == 565] = 1

#fix japan
duplicate_dict3[445] = 649

#fix korea
duplicate_dict3[678] = 446

#fix kosovo
duplicate_dict3[595] = 66

#fix Mexico
lc.loc[lc.ID==394, 'Longitude'] *= -1

#fix Morocco
duplicate_dict3[347] = 6

#fix Slovenia
duplicate_dict3[596] = 631

#fix SA
duplicate_dict3[628] = 361

#fix Spain
lc.not_most_accurate.loc[lc.ID == 55] = 1

#fix Sweden
duplicate_dict3[99] = 583

#fix Taiwan
duplicate_dict3[759] = 447
duplicate_dict3[751] = 447
duplicate_dict3[752] = 448
duplicate_dict3[760] = 448

#fix Turkey
duplicate_dict3[51] = 552
duplicate_dict3[591] = 552

#fix US
duplicate_dict3[639] = 607

#fix Russia
duplicate_dict3[160] = 156
duplicate_dict3[408] = 285
lc.not_most_accurate.loc[lc.ID == 587] = 1

#fix Czech
duplicate_dict3[487] = 582

#fix Tunisia
duplicate_dict[793] = 370

#fix UK
duplicate_dict[588] = 490


lc = lc[np.logical_not(lc.ID.isin(duplicate_dict3.keys()))]
lc = lc[lc.not_most_accurate < 0.1]
lc = lc[np.logical_not(np.isnan(lc.Latitude))]
lc = lc[np.logical_not(np.isnan(lc.Longitude))]

lc.to_csv('regions/locations_deduplicated.csv')

duplicate_dict.update(duplicate_dict2)
duplicate_dict.update(duplicate_dict3)
np.savetxt('duplicate_dict.txt', 
           np.array([(k,v) for k,v in  duplicate_dict.items()]), fmt="%d")

