#Tumor_rank.R
rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggpubr)

#1.GMTD----
mat=read.csv('10_GMTD_2.csv',header = T)#20
mat$Short_name=factor(mat$Short_name,levels = mat$Short_name)

# 1.1（Score_1 and Short_name）
p1=ggdotchart(mat, x = "Short_name", y = "Score_1",color = "#FC4E07",
              add = "segments",                            
              ggtheme = theme_pubr())+
  scale_x_discrete(limits = rev(mat$Short_name)) +  
  labs(x = '', y = 'GMTD score')+coord_flip()

# 1.2（nproject and Short_name）
p2 <- ggplot(mat, aes(x = factor(Short_name,levels = rev(Short_name)), y = nsample, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge",width=0.8, color = "black") +
  geom_text(aes(label=label2),vjust = 0.5,hjust=0,size=4)+
  labs(x = '', y = 'No.samples(batches)') +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +scale_y_continuous(trans = "log1p")+
  scale_fill_manual(values = c("GI" ='#ECC97F', "non-GI" = "#8FC9E2")) + 
  coord_flip() 

P=cowplot::plot_grid(p2,p1,align = "vh",nrow=1)
ggsave(P,filename = 'F4_GMTD.pdf',height = 6,width = 10)

#2.GMTP----
mat=read.csv('10_GMTP_2.csv',header = T)#4
mat$Short_name=factor(mat$Short_name,levels = mat$Short_name)

p1=ggdotchart(mat, x = "Short_name", y = "Score_1",color = "#FC4E07",
              add = "segments",                            
              ggtheme = theme_pubr())+
  scale_x_discrete(limits = rev(mat$Short_name))+
  scale_y_continuous(breaks = seq(0, 15, by = 5)) + 
  labs(x = '', y = 'GMTP score')+coord_flip()

p2 <- ggplot(mat, aes(x = factor(Short_name,levels = rev(Short_name)), y = nsample, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge",width=0.8, color = "black") +
  geom_text(aes(label=label2),vjust = 0.5,hjust=0,size=4)+
  labs(x = '', y = 'No.samples(batches)') +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +scale_y_continuous(trans = "log1p")+
  scale_fill_manual(values = c("GI" ='#ECC97F', "non-GI" = "#8FC9E2")) + 
  coord_flip() 

P=cowplot::plot_grid(p2,p1,align = "vh",nrow=1)
ggsave(P,filename = 'F4_GMTP.pdf',height = 4,width = 10)
