# Introduction

The code is used for the work of “HGMT: a database of human gut microbiota for tumors and immunotherapy response”. 



# Pipeline

* WGS.sh, taxonomic assignment for WGS data (Bowtie2, Trimmomatic, MetaPhlAn4 and Kraken2)

* 16S.sh,  taxonomic assignment for 16S data (QIIME2, DADA2, Greengenes2 database)

* ITS.sh,  taxonomic assignment for ITS data (QIIME2, DADA2, UNITE database v10.0)

  

# Analysis

* Data_statistics.R, data summary of HGMT regarding the distribution of tumor types, assay types, meta variables, etc.
* Pan-tumor_DA_taxa.R, the differential abundance markers of tumor diagnosis and immunotherapy response prediction in pan-tumor
* Tumor_rank.R, tumor rank based on GMTD and GMTP score 



# Contact information

If you have any questions about this code, please contact the corresponding author (xmzhao@fudan.edu.cn) or the first author (21110850029@m.fudan.edu.cn).

