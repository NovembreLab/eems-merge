#!/bin/bash
SCRIPT=scripts/
RAW=raw/
TMP=tmp/


function asc_set {
    python $SCRIPT/convert2ped.py $RAW/EuropeAllData/panel$1/vdata
    sed -i 's/^90/24/' $RAW/EuropeAllData/panel$1/vdata.bim
    python $SCRIPT/get_pos_bim.py $RAW/EuropeAllData/panel$1/vdata $TMP/asc${1}_lazaridis0
    plink --bfile $TMP/asc${1}_lazaridis0 --make-bed --out $TMP/asc${1}_lazaridis1 --exclude $RAW/multiallelic.txt --autosome
    gawk -i inplace 'BEGIN{OFS="\t"}; { print $2, $2, $3, $4, $5, $6 }' $TMP/asc${1}_lazaridis1.fam
    plink --bfile $TMP/asc${1}_lazaridis1 --make-bed --out $TMP/asc${1}_lazaridis2 --keep input_files/valid_names.txt --remove tmp/pca-outliers.txt
}

asc_set $1
