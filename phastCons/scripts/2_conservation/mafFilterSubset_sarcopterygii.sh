
echo "hola mundo"


#ok: 1 2 3 4 5 6 7
#ok: 8 9 10 11 12
#ok: 13 14 15 16 17 18 19 20 21 22

source_remote_base_path='s3://acc-regions-2018/download/alignments/multiz100way/'
target_remote_base_path='s3://acc-regions-2018/preparation/sarcopterygii/'
local_base_path='/home/rstudio/disco_tmp/'
bin_base_path='/home/rstudio/soft'

begin=1
end=8

for i in $(seq $begin $end)
do

  chr="chr"$i
  echo $chr
  
  sudo aws s3 cp $source_remote_base_path$chr.maf.gz $local_base_path
  
  sudo gunzip $local_base_path$chr.maf.gz $local_base_path$chr.maf
  
  sudo $bin_base_path/mafSpeciesSubset $local_base_path$chr.maf ./species_sarcopterygii.lst $local_base_path$chr'_sarcopterygii.maf'
  
  sudo gzip $local_base_path$chr'_sarcopterygii.maf' $local_base_path$chr'_sarcopterygii.maf.gz'
  
  aws s3 cp $local_base_path$chr'_sarcopterygii.maf.gz' $target_remote_base_path
  
  sudo rm $local_base_path$chr'.maf'
  
  #sudo rm $local_base_path$chr'_sarcopterygii.maf'
  
  sudo rm $local_base_path$chr'_sarcopterygii.maf.gz'

done