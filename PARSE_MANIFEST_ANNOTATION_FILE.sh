# parsing the manifest annotation file
cd /gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/data/annotation/manifest
wget https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip
unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip 
tail -n+8 infinium-methylationepic-v-1-0-b5-manifest-file.csv > manifest.csv
tail -n+8 MethylationEPIC_v-1-0_B4.csv > manifest_4.csv
