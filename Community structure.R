rm(list=ls()) #Delete all variables in the current environment

library(tidyverse)#Data collation and data conversion package, using some functions that are easier to use and understand
library(ggprism)#The ggprism package provides various themes, color palettes, and other useful functions to customize ggplots, and provides them with
library(vegan)
otu <- read.csv('D:/Liu Yi-Nan/Desktop/Liu Yinan data/otu_mat_filtered.csv',row.names = 1)
otu=t(otu)
head(otu, n = 3)
tax <- read.csv('D:/Liu Yi-Nan/Desktop/Liu Yinan data/tax_mat_filtered_del.csv',row.names = 1)
head(tax, n = 3)
metadata<- read.csv('D:/Liu Yi-Nan/Desktop/Liu Yinan data/sample_data_update.csv',row.names = 1)
head(metadata, n = 3)
dat <- merge(x=otu,y=tax,by='row.names') 
head(dat, n = 3)
dat =dplyr::rename(dat,OTUID = Row.names)
head(dat, n = 3)


####Group summary by genus level (change the species level to be displayed according to your own needs)
aa<-aggregate(dat[,1:ncol(otu)+1],by=list(dat$Genus),FUN=sum)
#Divide the dataframe according to the variable behind by (automatically as a factor),
#and bring it into the Fun function behind, and return different head(aa)
row.names(aa)=aa$Group.1   
head(aa)
aa<-dplyr::select(aa,-Group.1)
head(aa, n = 3)
#Sort data based on row sum results
order<-sort(rowSums(aa),index.return=TRUE,decreasing=T)   #Sort from largest to smallest and return
#Sort the table according to the column summation result, the two return value lists of order, one is the arranged sequence,
#and the other is its serial number position in the original list
cc<-aa[order$ix,]
head(cc, n = 3)
write.csv(cc,'D:/Liu Yi-Nan/Desktop/Liu Yinan data/community_str_genus.csv',quote=FALSE)
write.table(cc,'D:/Liu Yi-Nan/Desktop/Liu Yinan data/community_str_genus.txt')






install.packages("ggalluvial")
library(tidyverse)
library(ggpubr)
library(ggalluvial)

aa=read.csv('D:/code/Liu Yinan data/community_str_genus.csv',row.names=1)
metadata=read.csv('D:/code/Liu Yinan data/sample_data_update.csv',row.names=1)

Kind_Position=vector()
Full_char=vector()

for(i in 1:nrow(metadata)){
  Kind_Position[i]=paste(metadata$Kind[i],metadata$Position[i])
  Full_char[i]=paste(metadata$Kind[i],metadata$Position[i],metadata$Stage[i])
}

metadata_plus=mutate(metadata,Kind_Position,Full_char)


n_select=which(metadata_plus$Kind_Position=="HDPE adhering"|metadata_plus$Kind_Position=="HDPE original")
aa_select=aa[,n_select]
order<-sort(rowSums(aa_select),index.return=TRUE,decreasing=T)   ##Sort from big to small and return
#Sort the table according to the column summation results, the two return value lists of order, one is the arranged sequence,
#and the other is its serial number position in the original list
cc<-aa_select[order$ix,]
dd<-rbind(colSums(cc[9:nrow(cc),]),cc[8:1,])#Get the dimensions nrow and length of the dataframe
head(dd, n = 3)
rownames(dd)[1]<-"Others"
head(dd, n = 3)
bb<-merge(t(dd),dplyr::select(metadata_plus,Stage),##Then merge with metadata
          by.x = "row.names",by.y ="row.names")#Select selects several columns, 
#its first parameter is the selection range (that dastaframe), and the latter is the column to be selected
#by.x by.y are the interface variables when the two dataframes are aligned and spliced respectively
head(bb, n = 3)##wide data variable length data
bb_value=select(bb,3:ncol(bb)-1)
bb_value_nor=bb_value/rowSums(bb_value)
bb_value_nor$Stage=bb$Stage

n=nrow(bb_value_nor)
bb_value_nor$Index=as.character(c(1:n))
ccol=colnames(bb_value_nor)[1:9]
bb_value_long<-gather(bb_value_nor,Genus,Abundance,-Stage,-Index)#The syntax is that the variables in the brackets include, the original data frame
##Add column, keep column
Yanse=c('#2F4F4F','#FF6C65','#00B800','#D39100','#ACE278','#00B8DF','#619BFA','#DB72F6','#FF61BF')
##sorted stacked column chart
ggplot(bb_value_long,aes(x = Index,y = Abundance,fill = Genus)) + #If the grouping information is not changed to a factor, the order of the legend of the character type is alphabeta,
  #which does not match the graphics
  geom_bar(colour='white',alpha=0.85,position="stack",stat = "identity",width = 0.1,size=0.1)+
  facet_grid(~Stage)+
  scale_fill_manual(values = Yanse)+
  labs(x='Incubation stage',y='Abundance')+
  scale_x_discrete(labels = c("Source","Early","Middle","Late"),breaks=c('Source','Early stage','Middle stage','Late stage'),limits=c('Source','Early stage','Middle stage','Late stage'))+
  guides(fill=guide_legend(reverse = F))+
  ggprism::theme_prism(border = TRUE)+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.05))+
  ggtitle('HDPE adhering')
p1
ggsave(file = "HDPE_Stage.jpg", plot =p1, dpi = 300,width = 6, height = 3.75 )






