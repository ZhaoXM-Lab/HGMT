export PATH=/home1/LiuJX/miniconda3/envs/qiime2-2021.4/bin:$PATH

#1.import dada
#Pair-end
qiime tools import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path manifest.txt \
 --input-format PairedEndFastqManifestPhred33V2 \
 --output-path 01_paired-end-demux.qza

#Single-end
#qiime tools import \
#  --type 'SampleData[SequencesWithQuality]' \
#  --input-path manifest.txt \
#  --output-path 01_paired-end-demux.qza \
#  --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
--i-data 01_paired-end-demux.qza \
--o-visualization 02_demux.qzv

#2.dada2
#Note: the following parameters are adjusted according to the quality report
#4 parameters: p-trim-left-f;p-trim-left-r;p-trunc-len-f;p-trunc-len-r

#Pair-end
qiime dada2 denoise-paired \
	--i-demultiplexed-seqs 01_paired-end-demux.qza \
	--p-trim-left-f 24 \
	--p-trim-left-r 10 \
	--p-trunc-len-f 320 \
	--p-trunc-len-r 245 \
	--p-n-threads 0 \
	--o-table 03_table.qza \
	--o-representative-sequences 03_rep-seqs.qza \
	--o-denoising-stats 03_denoising-stats.qza

#Single-end
#qiime dada2 denoise-single \
#	--i-demultiplexed-seqs 01_paired-end-demux.qza \
#	--p-trim-left 0 \
#	--p-trunc-len 150 \
#	--o-table 03_table.qza \
#	--o-representative-sequences 03_rep-seqs.qza \
#	--o-denoising-stats 03_denoising-stats.qza

#3.taxonomic assignment
#greengene2
#V4: GreenGene/2022.10.backbone.v4.nb.qza
#full length: GreenGene/2022.10.backbone.full-length.nb.qza
qiime feature-classifier classify-sklearn \
	--i-classifier GreenGene/2022.10.backbone.v4.nb.qza \
	--i-reads 03_rep-seqs.qza \
	--o-classification 05_taxonomy.qza

#collapse
#Genus level
qiime taxa collapse \
	--i-table 03_table.qza \
	--i-taxonomy 05_taxonomy.qza \
	--p-level 6 \
	--o-collapsed-table 08_table_level6.qza

qiime tools export \
	--input-path 08_table_level6.qza \
	--output-path 08_table_level6

cd 08_table_level6
biom convert \
	-i feature-table.biom \
	-o genus.tsv \
	--to-tsv

