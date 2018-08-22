#!/bin/bash

START=12
END=22
base_url_head='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/multiz100way/maf/'
file_extension='.maf.gz' 
local_base_path='/home/rstudio/data/download/alignments/multiz100way'
s3_base_path='s3://acc-regions-2018/download/alignments/multiz100way'

cd $local_base_path

for i in $(seq $START $END); do 

file_name='chr'$i;
file_url=$base_url_head$file_name$file_extension;
echo $file_url 
wget $file_url
aws s3 cp $local_base_path/$file_name$file_extension $s3_base_path/$file_name$file_extension 
rm $local_base_path/$file_name$file_extension

done

# file_name='chrM'
# file_url=$base_url_head$file_name$file_extension;
# echo $file_url 
# wget $file_url
# aws s3 cp $local_base_path/$file_name$file_extension $s3_base_path/$file_name$file_extension 
# rm $local_base_path/$file_name$file_extension
# 
# file_name='chrX'
# file_url=$base_url_head$file_name$file_extension;
# echo $file_url 
# wget $file_url
# aws s3 cp $local_base_path/$file_name$file_extension $s3_base_path/$file_name$file_extension 
# rm $local_base_path/$file_name$file_extension
# 
# file_name='chrY'
# file_url=$base_url_head$file_name$file_extension;
# echo $file_url 
# wget $file_url
# aws s3 cp $local_base_path/$file_name$file_extension $s3_base_path/$file_name$file_extension 
# rm $local_base_path/$file_name$file_extension

#http://hgdownload.cse.ucsc.edu/goldenPath/hg38/multiz100way/maf/chr1.maf.gz