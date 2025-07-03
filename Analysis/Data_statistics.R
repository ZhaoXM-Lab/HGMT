#Data_statistics.R

rm(list=ls())

library(dplyr)
library(ggplot2)

#1.project/run distributions----
sample_meta=read.csv('04_sample_meta_all_QC_v2.csv',header=T)#18630
table(sample_meta$QC_state)#17323+68 processed runs
table(sample_meta$Phenotype_name)
table(sample_meta$Assay_type,useNA = 'always')
# 16S     ITS   WGS  
# 11356    69  7205  

xx=read.csv('01_disease_info.csv',header=T)
mat=sample_meta[,c("Project_ID",'Assay_type','Phenotype_name','Phenotype_ID')]%>%unique.data.frame()
mat=merge(mat,xx,by=c('Phenotype_name','Phenotype_ID'))
mat$count=1

#1.1 Disease types----
df2=mat[,c('Short_name','Assay_type','count')]%>%group_by(Short_name)%>%summarise(Count=sum(count))%>%ungroup()
df2 = df2[order(df2$Count, decreasing = TRUE),]
df2$Short_name
sum(df2$Count>2)#14-1 cancer
sum(df2$Count>1)#20-1
df2$Short_name=factor(df2$Short_name, levels=df2$Short_name)

P1 =  df2 %>% 
  ggplot(aes(x=factor(Short_name,levels=rev(Short_name)), y=Count)) +
  geom_bar(stat="identity",width=0.8,color="black",fill="#BD9AAD") + 
  geom_text(aes(label=Count),vjust = 0.5,hjust=0,size=4)+ 
  theme_classic() +labs(title='Number of projects',x='Phenotype')+
  theme(text=element_text(face = 'bold',size=10),
        axis.text.x=element_text(face = 'bold',size=10),
        axis.text.y=element_text(face = 'bold',size=10),
        plot.title=element_text(hjust=0.5)) + 
  scale_y_continuous(trans = "log1p")+coord_flip()
P1
ggsave(P1,filename = 'F2_project_distribution.pdf',height = 7,width = 5)

#1.2 Assay type----
mat=sample_meta[,c("Project_ID",'Assay_type')]%>%unique.data.frame()#56
df=table(mat$Assay_type)%>%as.data.frame()
# 16S ITS WGS 
# 74   3  45 
df = df[order(df$Freq, decreasing = TRUE),]
df$Var1=factor(df$Var1,levels=df$Var1)
label_value <- paste('(', round(df$Freq/sum(df$Freq) * 100, 2), '%)', sep = '')
df$label<- paste(df$Var1, label_value, sep = '\n')

P2=df%>%mutate(label=factor(label,levels=df$label))%>%
  ggplot(aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4, color = "black") + 
  labs(x = '', y = '', title = 'Assay type')+
  coord_polar("y", start = 0)+
  theme_void()+
  theme(text= element_text(face = 'bold',size = 10),
        plot.title=element_text(hjust=0.5),
        legend.position = 'none',
        axis.line = element_blank(),       
        axis.ticks = element_blank(),      
        axis.text = element_blank())+
  scale_fill_brewer(palette="Set2")
P2
ggsave(P2,filename = 'F2_project_assay_distribution.pdf',height = 4,width = 4)

#2.meta----
cbPalette <- c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A","#E31A1C","#FB9A99","#FF7F00","#FDBF6F","#6A3D9A","#CAB2D6","#FFFF99","#666666")

df=sample_meta[,c('Age','Sex','BMI','Country')]
kk=1-colSums(is.na(df))/nrow(df)
df=data.frame(names=names(kk),percentages=round(kk*100,digits = 2))

P3 = df %>%
  arrange(desc(percentages)) %>%
  mutate(names=factor(names,levels = names))%>%
  ggplot(aes(x=names, y=percentages)) +
  geom_bar(stat="identity", aes(fill=names),color="black",fill="#BD9AAD",width = 0.8) +
  geom_text(aes(label=percentages),size=4,vjust = 0.5,hjust=0)+ 
  theme_classic() + scale_y_continuous(limits=c(0,100),breaks = seq(0,100,20))+
  theme(text=element_text(face = 'bold',size=10),
        axis.text.x=element_text(hjust=1,face = 'bold',size=10),
        axis.text.y=element_text(face = 'bold',size=10),
        plot.title=element_text(hjust=0.5),legend.position = 'none') + coord_flip()+
  labs(title="",x='meta information',y='% samples') 
P3

#The frequency of co-occurrence
df=sample_meta[,c('Age','Sex','BMI','Country')]
kk=rowSums(is.na(df))
kk=table(kk)
df=data.frame(names=names(kk),Freq=as.numeric(kk),label=c('None','Only one','With two','With three','All four'))

df = df[order(df$Freq, decreasing = TRUE),]
label_value <- paste('(', round(df$Freq/sum(df$Freq) * 100, 2), '%)', sep = '')
df$label<- paste(df$label, label_value, sep = ' ')

P4=df%>%mutate(label=factor(label,levels=df$label))%>%
  ggplot(aes(x = "", y = Freq, fill = label)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  labs(x = '', y = '', title = '',fill='Nr. of available meta variables: \n Age,Sex,BMI,Country')+
  coord_polar("y", start = 0)+
  scale_fill_manual(values=cbPalette)+
  theme_void()+
  theme(text= element_text(face = 'bold',size = 10),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", color = "black"),
        plot.title=element_text(hjust=0.5),
        axis.line = element_blank(),       
        axis.ticks = element_blank(),    
        axis.text = element_blank())
P=cowplot::plot_grid(P3,P4,nrow=1)
P
ggsave(P,filename = 'F2_meta_percentages.pdf',height = 5,width = 10)

#3.runs:ICI----
mat=subset(sample_meta,grepl('^response:',Sample_description))#2207
mat=merge(mat,xx,by=c('Phenotype_name','Phenotype_ID'))
table(mat$Short_name,mat$Assay_type)
table(mat$Assay_type)
# 16S  WGS 
# 584 1623 

mat$count=1
df=mat[,c('Short_name','Assay_type','count')]%>%group_by(Short_name,Assay_type)%>%summarise(Count=sum(count))%>%ungroup()
unique(df$Short_name)
df$Assay_type=factor(df$Assay_type,levels=c('16S','WGS'))

df2=mat[,c('Short_name','Assay_type','count')]%>%group_by(Short_name)%>%summarise(Count=sum(count))%>%ungroup()
df2 = df2[order(df2$Count, decreasing = TRUE),]

library(ggbreak)
P6 =  df %>% 
  ggplot(aes(x=factor(Short_name,levels=df2$Short_name), y=Count,group=Assay_type)) +
  geom_bar(stat="identity",position='stack', aes(fill=Assay_type),color="black") +
  geom_text(aes(label=Count),position=position_stack(),size=4)+
  theme_classic() +labs(title='Number of runs (with ICI therapy)',x='Phenotype')+
  theme(text=element_text(face = 'bold',size=10),
        axis.text.x=element_text(angle=90, hjust=1,vjust=0.5,face = 'bold',size=10),
        axis.text.y=element_text(face = 'bold',size=10),
        plot.title=element_text(hjust=0.5)) + scale_fill_brewer(palette="Set2")
P6
ggsave(P6,filename = 'F2_ICI_runs_distribution.pdf',height = 5,width = 9)

library(RColorBrewer)
library(stringr)
brewer.pal(3, "Set2")
col_assay=c('16S'="#66C2A5",'WGS'="#FC8D62",'ITS'="#8DA0CB")
P6_2 = df %>% subset(Short_name %in% df2$Short_name[-(1:7)]) %>% 
  mutate(Short_name=factor(Short_name, levels=df2$Short_name[-(1:7)])) %>%
  ggplot(aes(x=Short_name, y=Count,group=Assay_type)) +
  geom_bar(stat="identity",position='stack', aes(fill=Assay_type),color="black") +
  geom_text(aes(label=Count),position=position_stack(),size=5)+
  theme_classic() +labs(title='',x='',y='')+
  theme(text=element_text(face = 'bold',size=7),
        axis.text.x=element_text(angle=90, hjust=1,vjust=0.5,face = 'bold',size=7),
        axis.text.y=element_text(face = 'bold',size=7),) + 
  theme(panel.border = element_blank(), axis.line = element_line(),
        panel.background = element_rect(fill = 'grey90'),
        legend.position = 'none')+
  scale_x_discrete(labels = function(x) str_sub(x, 1, 1))+
  scale_fill_manual(values=col_assay)
P6_2
ggsave(P6_2,filename = 'F2_ICI_runs_2.pdf',height = 2.5,width = 3)
