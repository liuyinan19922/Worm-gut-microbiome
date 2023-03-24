library(phyloseq)
library(ggplot2)
library(ggprism)
library(vegan)
library(gridExtra)#gridExtra: Miscellaneous Functions for "Grid" Graphics
library(picante)#integrating Phylogenies and Ecology

###########  Analytical data preprocessing and import   ###########

setwd("D:/Liu Yi-Nan/Desktop/Re__NMDS_results")#As the name suggests, it is used to set the current working directory. 
#Note that the calling format is: setwd("target path"), remember to add double quotes
sample <- read.delim('metadata.txt', row.names=1)
sample$Time=factor(sample$Time)
Tax <- read.delim('Classification.txt', row.names=1) 
OTU_mat <- as.matrix(ASV)
tax_mat <- as.matrix(Tax)
dim(OTU_mat)
OTU_f = otu_table(OTU_mat, taxa_are_rows = FALSE)
TAX = tax_table(tax_mat)
samples = sample_data(sample)
NAC_phy_master <- phyloseq(OTU_f, TAX, samples) #Make the phyloseq object
set.seed(100)#Set the random number seed to ensure that the results can be repeated each time
rar <- rarefy_even_depth(NAC_phy_master, sample.size = min(sample_sums(NAC_phy_master)),
                         rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)#Rarefy to even depth
write.csv(otu_table(rar), file = "rarefied_otu.csv")
rarecurve(otu_table(rar), step=50, cex=0.5)
Fac1 <- sample
head(Fac1)
OTU <- otu_table(pe_rar)#Storing the rarefied abundacne matrix
OTUt <- data.frame(OTU) # Column name should be genera and samples should be row name
head(OTUt)
dim(OTUt)
x<- row.names(OTUt)
OTUt$Sample <- x # Adding sample names as a column to merge for subsetting based on kind of plastic later
y <- row.names(Fac1)
Fac1$Sample <- y
FOR_SUB <- merge(OTUt, Fac1, by = "Sample", all = TRUE) #merge metadata and count #如果有同名的by#如果没有同名by.x by.y
head(FOR_SUB)#It is a combination of an otu form plus a taxon form
dim(FOR_SUB)
dim(Fac1)
dim(OTUt)
phy <- read.tree("tree_asv")
class(phy)
phy
plot(phy, cex = 0.5)
phy.dist <- cophenetic(phy)


######### Select grouping for analysis  pe #########
Prefix='HDPE'
beta.reps = 999# number of randomizations
plastic <- subset(FOR_SUB, Kind == Prefix) #Subset by kind of different plastics

head(FOR_SUB)#It is a combination of an otu form plus a taxon form
dim(FOR_SUB)#Split the composite table
dim(plastic)
rn_plastic <- plastic$Sample
row.names(plastic) <- rn_plastic # Add row names back
rn_plastic_pic <- subset(plastic, select = -c(Sample, Kind,Position, Time,Stage)) #Remove the metadata
plastic_meta <- subset(plastic, select = c(Sample, Kind,Position,Time,Stage)) #Store metadata for later use
keep <- colSums(rn_plastic_pic) >=100#Filter ASV based on count

rn_plastic_pic_f <- rn_plastic_pic[,keep]
dim(rn_plastic_pic_f)
apply(rn_plastic_pic_f, 1, sum) #Check sum of samples
#The second parameter of #apply specifies which kind of margin it is, 1 corresponds to each column, and 2 corresponds to each row. 
#That is the so-called margin #Other several apply family functions#lapply and sapply operate similarly, 
#but are suitable for vector and list, so there is no need to choose margin, only input data
#And the function operation performed (this function operation uses each element as an input variable, rather than the previous input of a whole line)
#tapply is to select a variable in the dataframe and perform a grouping operation on a variable,
# Then use the first parameter as the input value to perform the operation of the last function
plastic_work <- decostand(rn_plastic_pic_f, method = "total") #Relative abundance calculation 
##Standardization Methods for Community Ecology total is normalized by dividing by row sum 
#total: divide by margin total (default MARGIN = 1).
apply(plastic_work, 1, sum)
head(plastic_work)

match.phylo.otu.plastic = match.phylo.data(phy,data.frame(t(plastic_work))) #compare taxa present in phylogenies with community or trait data sets, 
#pruning and sorting the two kinds of data to match one another for subsequent analysis.
str(match.phylo.otu.plastic) #check classification of viriables
beta.mntd.weighted.plastic = as.matrix(comdistnt(t(match.phylo.otu.plastic$data),
                                                 cophenetic(match.phylo.otu.plastic$phy),abundance.weighted=TRUE,exclude.conspecifics = FALSE))
dim(beta.mntd.weighted.plastic)
beta.mntd.weighted.plastic[1:5,1:5] #write.csv(beta.mntd.weighted,'betaMNTD_weighted.csv',quote=F);

identical(colnames(match.phylo.otu.plastic$data),colnames(beta.mntd.weighted.plastic)) # just a check, should be TRUE
identical(colnames(match.phylo.otu.plastic$data),rownames(beta.mntd.weighted.plastic)) #Judging whether two objects are equal in R language 
#There are also functions with similar functions such as all.equal(x,y)
rand.weighted.bMNTD.comp.plastic = array(c(-999),dim=c(ncol(match.phylo.otu.plastic$data),
                                                       ncol(match.phylo.otu.plastic$data),beta.reps));#array R language creates an array, 
#and the input of the function is data and dimensions
dim(rand.weighted.bMNTD.comp.plastic);#The semicolon in the #r language indicates the end of the statement. 
#If it has already been divided into lines, the semicolon is optional

for (rep in 1:beta.reps) {
  rand.weighted.bMNTD.comp.plastic[,,rep] = as.matrix(comdistnt(t(match.phylo.otu.plastic$data),
                                                                taxaShuffle(cophenetic(match.phylo.otu.plastic$phy)),abundance.weighted=TRUE,exclude.conspecifics = FALSE));
  #taxaShuffle  Matrix with taxa names shuffled 
  print(c(date(),rep));
  #print prints its argument and returns it invisibly
}


weighted.bNTI.plastic = matrix(c(NA),nrow=ncol(match.phylo.otu.plastic$data),
                               ncol=ncol(match.phylo.otu.plastic$data))
dim(weighted.bNTI.plastic)

for (columns in 1:(ncol(match.phylo.otu.plastic$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu.plastic$data)) {
    
    rand.vals.plastic = rand.weighted.bMNTD.comp.plastic[rows,columns,]
    weighted.bNTI.plastic[rows,columns] = (beta.mntd.weighted.plastic[rows,columns] - mean(rand.vals.plastic)) / sd(rand.vals.plastic)
    rm("rand.vals.plastic")#This line deletes the variable from memory
    
  }
}

rownames(weighted.bNTI.plastic) = colnames(match.phylo.otu.plastic$data)
colnames(weighted.bNTI.plastic) = colnames(match.phylo.otu.plastic$data)
weighted.bNTI.plastic
#write.csv(weighted.bNTI,"weighted_bNTI.csv",quote=F);

hist(weighted.bNTI.plastic)
weighted.bNTI_df.plastic <- as.data.frame(weighted.bNTI.plastic)
weighted.bNTI_df.plastic$Sample <- row.names(weighted.bNTI_df.plastic) 
head(plastic_meta)
dim(plastic_meta)

weighted.bNTI_df_merged.plastic <- merge(weighted.bNTI_df.plastic, plastic_meta, by = "Sample")
dim(weighted.bNTI_df_merged.plastic)

write.csv(weighted.bNTI_df_merged.plastic,"weighted_bNTI_merged_pe_RAREFIED_WEIGHTED_try.csv",quote=F)

###########Import data for plotting##############

table.pe=read.csv('D:/Liu Yi-Nan/Desktop/Re__NMDS_results/weighted_bNTI_merged_pe_RAREFIED_WEIGHTED_try.csv',row.names = 1)
metadata=read.delim('D:/Liu Yi-Nan/Desktop/Re__NMDS_results/metadata.txt',row.names=1)
Table.pe=merge(table.pe[,1:52],metadata,by.x="Sample",by.y="row.names")
n.early=which(Table.pe$Stage=='Early stage'&Table.pe$Position=='adhering')
n.middle=which(Table.pe$Stage=='Middle stage'&Table.pe$Position=='adhering')
n.late=which(Table.pe$Stage=='Late stage'&Table.pe$Position=='adhering')
n.source=which(Table.pe$Stage=='Source')



n.pair.source.early=expand.grid(n.source,n.early)
bNTI.source.early=vector()
for (i in 1:nrow(n.pair.source.early))
{     if(n.pair.source.early[i,1]<n.pair.source.early[i,2])
{
  temp=n.pair.source.early[i,1]
  n.pair.source.early[i,1]=n.pair.source.early[i,2]
  n.pair.source.early[i,2]=temp
  bNTI.source.early[i]=Table.pe[n.pair.source.early[i,1],n.pair.source.early[i,2]+1]
}
  else{
    bNTI.source.early[i]=Table.pe[n.pair.source.early[i,1],n.pair.source.early[i,2]+1]
  }
}
order1=Table.pe$Order[n.pair.source.early[,1]]
order2=Table.pe$Order[n.pair.source.early[,2]]
bNTI.source.early=data.frame(value=bNTI.source.early,name=rep('Early stage',length(bNTI.source.early)),order1,order2)
bNTI.source.early=subset(bNTI.source.early,(order1==order2))
summary(bNTI.source.early)




n.pair.source.middle=expand.grid(n.source,n.middle)
bNTI.source.middle=vector()
for (i in 1:nrow(n.pair.source.middle))
{     if(n.pair.source.middle[i,1]<n.pair.source.middle[i,2])
{
  temp=n.pair.source.middle[i,1]
  n.pair.source.middle[i,1]=n.pair.source.middle[i,2]
  n.pair.source.middle[i,2]=temp
  bNTI.source.middle[i]=Table.pe[n.pair.source.middle[i,1],n.pair.source.middle[i,2]+1]
}
  else{
    bNTI.source.middle[i]=Table.pe[n.pair.source.middle[i,1],n.pair.source.middle[i,2]+1]
  }
}  
order1=Table.pe$Order[n.pair.source.middle[,1]]
order2=Table.pe$Order[n.pair.source.middle[,2]]
bNTI.source.middle=data.frame(value=bNTI.source.middle,name=rep('Middle stage',length(bNTI.source.middle)),order1,order2)
bNTI.source.middle=subset(bNTI.source.middle,(order1==order2))
summary(bNTI.source.middle)


n.pair.source.late=expand.grid(n.source,n.late)
bNTI.source.late=vector()
for (i in 1:nrow(n.pair.source.late))
{     if(n.pair.source.late[i,1]<n.pair.source.late[i,2])
{
  temp=n.pair.source.late[i,1]
  n.pair.source.late[i,1]=n.pair.source.late[i,2]
  n.pair.source.late[i,2]=temp
  bNTI.source.late[i]=Table.pe[n.pair.source.late[i,1],n.pair.source.late[i,2]+1]
}
  else{
    bNTI.source.late[i]=Table.pe[n.pair.source.late[i,1],n.pair.source.late[i,2]+1]
  }
}  
order1=Table.pe$Order[n.pair.source.late[,1]]
order2=Table.pe$Order[n.pair.source.late[,2]]
bNTI.source.late=data.frame(value=bNTI.source.late,name=rep('Late stage',length(bNTI.source.late)),order1,order2)
bNTI.source.late=subset(bNTI.source.late,(order1==order2))
summary(bNTI.source.late)


colorassorted=c('#8997ab',"#ff6c67","#00bd0f","#4c9eff")#
bNTI=rbind(bNTI.source.early,bNTI.source.middle,bNTI.source.late)

pebNTI=ggplot(bNTI, aes(x = name, y = value)) + 
  geom_boxplot(width=0.6,fatten = NULL,lwd=1.5,color=colorassorted[2])+
  geom_hline(yintercept = 2,lwd=1.5,linetype="dashed",color='grey')+
  geom_jitter(width = 0.2)+
  stat_summary(fun.y = "mean", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width =0.6, size = 1, linetype = "solid",color=colorassorted[2])+
  labs(x='Incubation stage',y=expression(paste(beta*'NTI')))+
  ggtitle("bNTI Pe")+
  scale_x_discrete(limits = c("Early stage", "Middle stage", "Late stage"),breaks = c("Early stage", "Middle stage", "Late stage") ,labels =c("Early", "Middle", "Late"))+
  scale_y_continuous(expand = c(0,0),limits = c(-2,10))+
  theme_prism(border = TRUE)+
  guides(fill = "none")
pebNTI
ggsave(file = "D:/Liu Yi-Nan/Desktop/pebNTI.jpg", plot =pebNTI, dpi = 300,width = 3.4, height = 3.75 )

#####Data saving and statistical testing####
data=data.frame(bNTI$value,bNTI$name)
sd=tapply(data[,1],data[,2],sd)
mean=tapply(data[,1],data[,2],mean)

kresult=kruskal.test(data[,1]~data[,2], data =data)
kresult
wresult=pairwise.wilcox.test(data[,1],data[,2], p.adjust.method = "BH")
wresult
prefix='pe'
prefix01="_bNTI_data.csv"
prefix03='_bNTI_test_wilcxon.csv'
prefix1=paste0(prefix,prefix01)
prefix3=paste0(prefix,prefix03)
write.csv(data.frame(mean,sd),prefix1)
write.csv(wresult$p.value,prefix3)

