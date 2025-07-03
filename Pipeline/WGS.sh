rawdata_dir=$1
minlen=$2
pair=$3
select_file=$4
type=$5


export JAVA_HOME=/home1/LiuJX/00_bin/Java/jdk-22.0.2
export PATH=$PATH:$JAVA_HOME/bin
Bowtie2=/home1/LiuJX/00_bin/share_version/miniconda3/bin/bowtie2
TRIMMO_ADAPTOR_FILE_PE='/home1/LiuJX/00_bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa'
TRIMMO_ADAPTOR_FILE_SE='/home1/LiuJX/00_bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa'
Human_database=/home1/Laisenying/Tools/data/human/human2/hg38
BasePath=/home1/LiuJX/HGMTA_2023/01_datasets/$rawdata_dir
rawdata_dir=$BasePath/01_rawdata

cd $BasePath

echo 'step1: Bowtie2:rmhuman'
if [[ $pair == 'PE' ]];then
 mkdir -p ${BasePath}/01_rmhuman/
 for i in $(cat ${BasePath}/run.txt)
 do
   echo  $i
   $Bowtie2 -x $Human_database  \
     -p 32 \
     --very-sensitive \
     -1 ${rawdata_dir}/${i}_1.fastq \
     -2 ${rawdata_dir}/${i}_2.fastq \
     --un-conc-gz ${BasePath}/01_rmhuman/${i}_clean%.fq.gz \
     --al-conc-gz ${BasePath}/01_rmhuman/${i}_contam%.fq.gz 1>/dev/null 2>&1
 done
else
 mkdir -p ${BasePath}/01_rmhuman/
 for i in $(cat ${BasePath}/run.txt)
 do
   echo  $i
   $Bowtie2 -x $Human_database -p 32 --very-sensitive -U ${rawdata_dir}/${i}.fastq --un-gz ${BasePath}/01_rmhuman/${i}_clean.fq.gz --al-gz ${BasePath}/01_rmhuman/${i}_contam.fq.gz --quiet
 done
fi

echo 'step2: trimmomatic'
if [[ $type == 'rmhuman' ]];then
 mkdir -p ${BasePath}/02_trim

 if [[ $pair == 'PE' ]];then
  for sample in $(cat ${BasePath}/$select_file)
   do
    java -jar /home1/LiuJX/00_bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 \
${BasePath}/01_rmhuman/${sample}_clean1.fq.gz \
${BasePath}/01_rmhuman/${sample}_clean2.fq.gz \
${BasePath}/02_trim/${sample}_trim1.paired.fq ${BasePath}/02_trim/${sample}_trim1.unpaired.fq \
${BasePath}/02_trim/${sample}_trim2.paired.fq ${BasePath}/02_trim/${sample}_trim2.unpaired.fq \
ILLUMINACLIP:${TRIMMO_ADAPTOR_FILE_PE}:2:40:15 \
LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:${minlen}  -phred33 > ${BasePath}/02_trim/${sample}_trim.log
   done
 else
   for sample in $(cat ${BasePath}/$select_file)
   do
    java -jar /home1/LiuJX/00_bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 32 ${BasePath}/01_rmhuman/${sample}_clean.fq.gz  ${BasePath}/02_trim/${sample}_trim.paired.fq  ILLUMINACLIP:${TRIMMO_ADAPTOR_FILE_SE}:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:${minlen} -phred33 > ${BasePath}/02_trim/${sample}_trim.log
  done
 fi
else
 echo 'have trim files! '
fi

echo 'step3: metaphlan4'
mkdir -p ${BasePath}/03_metaphlan
export PATH=/home1/XUCT/miniconda3/envs/mpa/bin:$PATH

if [[ $pair == 'PE' ]];then
 cat ${BasePath}/$select_file|xargs -t -I {} metaphlan ${BasePath}/02_trim/{}_trim1.paired.fq,${BasePath}/02_trim/{}_trim2.paired.fq --bowtie2out ${BasePath}/03_metaphlan/{}.bowtie2.bz2 --bowtie2_exe $Bowtie2 --nproc 12 --input_type fastq -o ${BasePath}/03_metaphlan/{}_prorun.txt
else
 echo 'SE'
 cat ${BasePath}/$select_file|xargs -t -I {} metaphlan ${BasePath}/02_trim/{}_trim.paired.fq --bowtie2out ${BasePath}/03_metaphlan/{}.bowtie2.bz2 --bowtie2_exe  $Bowtie2 --nproc 12 --input_type fastq -o ${BasePath}/03_metaphlan/{}_prorun.txt
fi

#export PATH=/home1/LiuJX/miniconda3/envs/mpa/bin:$PATH
merge_metaphlan_tables.py ${BasePath}/03_metaphlan/*_prorun.txt > ${BasePath}/03_merged_Bacteria_SGB.txt


echo 'step4: kraken2+bracken'
DBNAME=/home1/LiuJX/00_bin/db/kraken
export PATH=/home1/LiuJX/miniconda3/bin:$PATH
mkdir -p ${BasePath}/04_fungus
mkdir -p ${BasePath}/04_fungus/new

if [[ $pair == 'PE' ]];then
 for i in $(cat ${BasePath}/run.txt)
 do
  echo  $i
  /home1/LiuJX/00_bin/kraken2/kraken2 --db $DBNAME --threads 12  --report $BasePath/04_fungus/${i}.report --output $BasePath/04_fungus/${i}.output --paired ${BasePath}/02_trim/${i}_trim1.paired.fq ${BasePath}/02_trim/${i}_trim2.paired.fq
  /home1/LiuJX/00_bin/Bracken/bracken -d $DBNAME -i $BasePath/04_fungus/${i}.report -o $BasePath/04_fungus/${i}.S.bracken -w $BasePath/04_fungus/${i}.S.bracken.report -r ${PE_len} -l S
  python /home1/LiuJX/00_bin/KrakenTools-master/kreport2mpa.py -r $BasePath/04_fungus/${i}.S.bracken.report -o $BasePath/04_fungus/new/${i}.report --display-header --percentages
 done
else
for i in $(cat ${BasePath}/run.txt)
do
  echo  $i
  /home1/LiuJX/00_bin/kraken2/kraken2 --db $DBNAME --threads 12  --report $BasePath/04_fungus/${i}.report --output $BasePath/04_fungus/${i}.output ${BasePath}/02_trim/${i}_trim.paired.fq
  /home1/LiuJX/00_bin/Bracken/bracken -d $DBNAME -i $BasePath/04_fungus/${i}.report -o $BasePath/04_fungus/${i}.S.bracken -w $BasePath/04_fungus/${i}.S.bracken.report -r ${PE_len} -l S
  python /share/home1/LiuJX/00_bin/KrakenTools-master/kreport2mpa.py -r $BasePath/04_fungus/${i}.S.bracken.report -o $BasePath/04_fungus/new/${i}.report --display-header --percentages
done
fi

python /home1/LiuJX/00_bin/KrakenTools-master/combine_mpa.py -i $BasePath/04_fungus/new/*.report -o $BasePath/04_merged_fungi.txt