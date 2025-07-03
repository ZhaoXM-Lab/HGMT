#Pan-tumor_DA_taxa.R
rm(list=ls())

library(dplyr)
library(ComplexHeatmap)
library(dendextend)
library(RColorBrewer)

#一、Tumor diagnosis----
#1.Data organization----
DA=read.csv('02_DA_results.csv',header = T)
table(DA$kingdom)
table(DA$level)

#marker taxa
temp_DA=DA[,c('taxid','level')]%>%unique.data.frame()
table(temp_DA$level)

temp_DA=DA[,c('taxa','level')]%>%unique.data.frame()
table(temp_DA$level)

#genus:ncbi_id,species:name
#total:1106
#genus：369
#species 737

DA=subset(DA,kingdom=='Bacteria'&Assay_type!='ITS'&level=='genus')
unique(DA$scientific_name)%>%length()#796
unique(DA$taxa)%>%length()#1756

kk=read.csv('./01_disease_info.csv',header = T)
DA=merge(DA,kk[,c(1:2)],by.x='Case',by.y='Phenotype_name',sort=F)
unique(DA$Short_name)#20
DA$Group=dplyr::case_when(DA$Case%in%c('Colorectal Neoplasms','Rectal Neoplasms',
                                       'Colonic Polyps','Adenomatous Polyps',
                                       'Stomach Neoplasms','Carcinoma, Pancreatic Ductal','Pancreatic Neoplasms',
                                       'Squamous Cell Carcinoma of Head and Neck','Esophageal Squamous Cell Carcinoma')~'GI',
                          TRUE~'non-GI')


#Bacteria: 16S/WGS
DA_1=subset(DA,display==T)#1923 genus
mat=reshape2::dcast(DA_1,scientific_short_name~DA.id, value.var = 'W_statistic',fill = 0,fun.aggregate = median)
rownames(mat)=mat$scientific_short_name
mat$scientific_short_name=NULL#319*73 

df=DA_1%>%group_by(scientific_short_name)%>%summarise(Nr.tumors=length(unique(Case)))
sum(df$Nr.tumors>=4)#84 genus
sum(df$Nr.tumors>=5)#58 genus

mat=mat[df$scientific_short_name[df$Nr.tumors>=4],]
mat=mat[,colSums(mat==0)!=nrow(mat)]
#84 genus*73 DA

#2.Visualization----
max(mat)#13.14894
min(mat)#-11.47035

#（1）W statistic
col_fun = circlize::colorRamp2(c(-8,-6,-4,-2,0,2,4,6,8),
                               c("#2E6E12", "#57C84D", "#83D475", "#C5E8B7",
                                 'white',
                                 "#F6BDC0", "#F1959B", "#EA4C46", "#781426"))

#（2）project
top_anno=DA_1[,c('DA.id','Assay_type','Short_name','Group')]%>%unique.data.frame()
stopifnot(top_anno$DA.id==colnames(mat))
mat=mat[,top_anno$DA.id]
top_anno$DA.id=NULL
colnames(top_anno)[2]='Tumor'

m1=c('16S'="#66C2A5",'WGS'="#FC8D62",'ITS'="#8DA0CB")
m2 = c("#D92B03","#C94E65","#6E86A5","#A5405E","#F4C288",
       "#088C00","#025939","#214EA7","#8B511F","#DB6C76",
       "#03C088","#7552A7","#DCA0DD","#F2E851","#F1B543",
       "#BF7533","#A38277","#592E13","#F2CDCF","#A77D9A")
names(m2)=unique(top_anno$Tumor)
m3=c("GI" ='#ECC97F', "non-GI" = "#8FC9E2")

top_annos = HeatmapAnnotation(df = top_anno, which = "column",
                              col=list(
                                Assay_type=m1,
                                Tumor=m2,
                                Group=m3
                              ))

#（3）taxa
right_anno=DA_1%>%group_by(scientific_short_name)%>%
  reframe(Nr.tumors=length(unique(Case)),
          disease_abbrev=paste0(unique(Short_name),collapse = ';'),
          Assay_type=paste0(unique(Assay_type),collapse = ';'),
          Phylum=sapply(strsplit(scientific_name,';'),function(x){gsub('p__','',x[2])}))%>%unique.data.frame()

table(right_anno$Phylum)
#Bacillota  Firmicutes
#Actinomycetota  Actinobacteria
#Bacteroidota  Bacteroidetes
#Pseudomonadota  Proteobacteria

right_anno$Phylum[!right_anno$Phylum%in%c('Bacillota','Actinomycetota','Bacteroidota','Pseudomonadota')]='Others'
table(right_anno$Phylum)

summary(right_anno$Nr.tumors)
right_anno=as.data.frame(right_anno)
rownames(right_anno)=gsub('g__','',right_anno$scientific_short_name)

rownames(mat)=gsub('g__','',rownames(mat))
right_annos=HeatmapAnnotation(df=right_anno[rownames(mat),c('Phylum','Nr.tumors')], which = "row",
                              col = list(Nr.tumors=circlize::colorRamp2(c(0,2,4,6,8,10), c("#FFFFFFFF","#ECF8F8FF","#C7EAEBFF","#9FDBDDFF","#73CDD0FF","#37BEC3FF")),
                                         Phylum=c('Bacillota'='#9BC985','Pseudomonadota'='#F7D58B',
                                                  'Bacteroidota'='#B595BF', 'Actinomycetota'='#797BB7',
                                                  'Others'='gray')))

#（4）taxa
dend = hclust(dist(mat))
dend = color_branches(dend, k = 5)

pdf('F3_GMTD_DA.pdf',width = max(ncol(mat)*0.2,7)-1,height = nrow(mat)*0.15+2)
Heatmap(mat,col=col_fun, width = unit(ncol(mat)*0.2, "cm"), height = unit(nrow(mat)*0.3, "cm"),
        column_title = 'Differential abundance taxa (tumor vs. health)',
        name = "W statistic",
        row_split = 5,cluster_rows = dend,cluster_columns=T,
        top_annotation=top_annos,right_annotation = right_annos,show_column_names = F
)
dev.off()
col_list=list(m1,m2,m3)
save(col_list,file='F3_GMTD_col_list.RData')

#3.Subset----
clusters <- cutree(dend, k = 5)
cluster_members <- list()
for (i in 1:5) {
  cluster_members[[i]] <- names(clusters[clusters == i])
}
save(cluster_members,file='F3_GMTD_cluster.RData')

write.csv(right_anno,'F3_GMTD_taxa_anno.csv',row.names = F)
write.csv(right_anno[c(cluster_members[[3]],cluster_members[[5]]),],'F3_GMTD_taxa_anno_healthy.csv',row.names = F)

#二、Tumor prognosis----
rm(list=ls())
#1.Data organization----
DA=read.csv('F3_DA_results_ICI.csv',header = T)
table(DA$kingdom)
table(DA$level)
DA=subset(DA,kingdom=='Bacteria'&Assay_type!='ITS'&level=='genus')

DA_1=subset(DA,display==T)#168 genus
unique(DA_1$Tumor)#4
mat=reshape2::dcast(DA_1,scientific_short_name~DA.id, value.var = 'W_statistic',fill = 0,fun.aggregate = median)
rownames(mat)=mat$scientific_short_name
mat$scientific_short_name=NULL#94*16 

df=DA_1%>%group_by(scientific_short_name)%>%summarise(Nr.tumors=length(unique(Tumor)))
sum(df$Nr.tumors>=2)#27 genus

mat=mat[df$scientific_short_name[df$Nr.tumors>=2],]
mat=mat[,colSums(mat==0)!=nrow(mat)]
#27*16 DA

#2.Visualization----
max(mat)#4.406114
min(mat)#-3.881063

#（1）W statistic
col_fun = circlize::colorRamp2(c(-8,-6,-4,-2,0,2,4,6,8),
                               c("#2E6E12", "#57C84D", "#83D475", "#C5E8B7",
                                 'white',
                                 "#F6BDC0", "#F1959B", "#EA4C46", "#781426"))

#（2）project
top_anno=DA_1[,c('DA.id','Assay_type','Tumor','Group')]%>%unique.data.frame()
stopifnot(top_anno$DA.id==colnames(mat))
mat=mat[,top_anno$DA.id]
top_anno$DA.id=NULL

load('F3_GMTD_col_list.RData')
m2=col_list[[2]]
rgb(red=241, green=101, blue=034, maxColorValue = 255)
rgb(red=244, green=148, blue=149, maxColorValue = 255)
m2=c(m2,'GIN'="#F16522",'RCC'="#F49495")
top_annos = HeatmapAnnotation(df = top_anno, which = "column",
                              col=list(Assay_type=col_list[[1]],Tumor=m2,Group=col_list[[3]]))

#（3）taxa
right_anno=DA_1%>%group_by(scientific_short_name)%>%
  reframe(Nr.tumors=length(unique(Tumor)),
          disease_abbrev=paste0(unique(Tumor),collapse = ';'),
          Assay_type=paste0(unique(Assay_type),collapse = ';'),
          Phylum=sapply(strsplit(scientific_name,';'),function(x){gsub('p__','',x[2])}))%>%unique.data.frame()

right_anno$Phylum[!right_anno$Phylum%in%c('Bacillota','Actinomycetota','Bacteroidota','Pseudomonadota')]='Others'
table(right_anno$Phylum)

summary(right_anno$Nr.tumors)
right_anno=as.data.frame(right_anno)
rownames(right_anno)=gsub('g__','',right_anno$scientific_short_name)

rownames(mat)=gsub('g__','',rownames(mat))
right_annos=HeatmapAnnotation(df=right_anno[rownames(mat),c('Phylum','Nr.tumors')], which = "row",
                              col = list(Nr.tumors=circlize::colorRamp2(c(0,2,4,6,8,10), c("#FFFFFFFF","#ECF8F8FF","#C7EAEBFF","#9FDBDDFF","#73CDD0FF","#37BEC3FF")),
                                         Phylum=c('Bacillota'='#9BC985','Pseudomonadota'='#F7D58B',
                                                  'Bacteroidota'='#B595BF', 'Actinomycetota'='#797BB7',
                                                  'Others'='gray')))

pdf('F3_GMTP_DA.pdf',width = max(ncol(mat)*0.2,7),height = nrow(mat)*0.15+2)
Heatmap(mat,col=col_fun, width = unit(ncol(mat)*0.2, "cm"), height = unit(nrow(mat)*0.3, "cm"),
        column_title = 'Differential abundance taxa (Responder vs. Non-responder)',
        name = "W statistic",row_split = 2,
        cluster_columns=T,cluster_rows = T,
        top_annotation=top_annos,right_annotation = right_annos,show_column_names = F
)
dev.off()


#3.Subset----
write.csv(right_anno,'F3_GMTP_taxa_anno.csv',row.names = F)
write.csv(right_anno[c('Methylobacterium',
                       'Faecalibacterium',
                       'Oscillibacter',
                       'Faecalibacillus',
                       'Fusicatenibacter',
                       'Roseburia'),],'F3_GMTP_taxa_anno_Responder.csv',row.names = F)

#三、Comparison between diagnosis markers and immunotherapy response prediction markers----
rm(list=ls())
Diagnosis=read.csv('F3_GMTD_taxa_anno.csv')
Prognosis=read.csv('F3_GMTP_taxa_anno.csv')

taxid=read.csv('06_taxa_final.csv',header = T)#06_taxa_taxid_update.R
unique(taxid$taxid[taxid$level=='genus'])%>%length()#1995,N.genus

library(VennDiagram)
a=venn.diagram(list(Diagnosis=Diagnosis$scientific_short_name,Prognosis=Prognosis$scientific_short_name),
               fill =  c("Diagnosis" = "#E8D283", "Prognosis" = "#9193B4"),
               cat.col = c("Diagnosis" = "#E8D283", "Prognosis" = "#9193B4"),
               lty = "dashed",
               filename=NULL)
grid.draw(a)
grid.newpage()
p1=phyper(89-1, 319, 1995-319, 94, lower.tail = F)#3.66715e-69

Diagnosis_sub=read.csv('F3_GMTD_taxa_anno_healthy.csv')
Prognosis_sub=read.csv('F3_GMTP_taxa_anno_responder.csv')
b=venn.diagram(list(Diagnosis_sub=Diagnosis_sub$scientific_short_name,Prognosis_sub=Prognosis_sub$scientific_short_name),
               fill =  c("Diagnosis_sub" = "#E8D283", "Prognosis_sub" = "#9193B4"),
               cat.col = c("Diagnosis_sub" = "#E8D283", "Prognosis_sub" = "#9193B4"),
               lty = "dashed",
               filename=NULL)
grid.draw(b)

intersect(Diagnosis$scientific_short_name,Prognosis$scientific_short_name)
# "g__Faecalibacterium" "g__Roseburia"        "g__Fusicatenibacter"
p2=phyper(3-1, 13, 1995-13, 6, lower.tail = F)#4.280125e-06
