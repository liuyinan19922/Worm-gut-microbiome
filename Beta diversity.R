source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
BiocManager::install("phyloseq")
BiocManager::install("vegan")
install.packages("ggalt")
library(phyloseq)
library(ggplot2)
library(vegan)
library(ggprism)
library(ggalt)
setwd("D:/Liu Yi-Nan/Desktop/summary/code/Liu Yinan data")
OTU <- read.csv('otu_mat_filtered.csv', row.names=1)# read.delim does not require all columns to be equal. 
#When read.delim is an empty string, you can specify the corresponding characters 
#or numbers to fill in for reading
sample <- read.csv('sample_data_update.csv', row.names=1)
Tax <- read.csv('tax_mat_filtered.csv', row.names=1) 

sample$'Time'=factor(sample$'Time')
OTU_mat <- as.matrix(OTU)
tax_mat <- as.matrix(Tax)

OTU_f = otu_table(OTU_mat, taxa_are_rows = FALSE)
TAX = tax_table(tax_mat)
samples = sample_data(sample)

head(OTU_f)
head(TAX)
head(samples)


NAC_phy_master <- phyloseq(OTU_f, TAX, samples)

rar <- rarefy_even_depth(NAC_phy_master, sample.size = min(sample_sums(NAC_phy_master)),
                         rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)


#HDPE community nmds
coloruse=c("#9A1207",'#FF6C67',"#FBB0AA")
NAC_phy_pe <- subset_samples(rar, Kind=='HDPE' )#subset data according to analysis
my.NMDS.bray_pe <- ordinate(NAC_phy_pe, "NMDS", "bray", trymax = 100)
stress=paste('pe',"Stress = ",as.character(round(my.NMDS.bray_pe[["stress"]],3)))
df_pe=plot_ordination(NAC_phy_pe, my.NMDS.bray_pe,type = "samples" ,justDF=TRUE)
po_pe=ggplot(df_pe,aes(x=NMDS1,y=NMDS2))+
  ggtitle(stress)+
  geom_point(aes(shape=Stage,color=Position),size =2, stroke = 1.5)+
  stat_ellipse(type = "t",level=0.68,aes(color=Position),lwd=0.8)+
  scale_shape_manual(values=c(16, 15, 17,18),breaks=c('Source','Early stage','Middle stage','Late stage'),labels=c('Source','Early','Middle','Late'))+
  scale_linetype_manual(values=c("dashed","twodash","solid"))+
  scale_color_manual(values=coloruse,breaks=c('original','adhering','planktonic'),labels=c('Source','Biofilm','Planktonic'))+
  theme_prism(border = TRUE)+
  theme(legend.text = element_text(size = 14))

po_pe
ggsave(file = "D:/Liu Yi-Nan/Desktop/beta_pe.jpg", plot =po_pe, dpi = 300,width = 5.625, height = 3.75 )

#############permanova############

group=data.frame(rownames(NAC_phy_ps@otu_table))
group$class=NAC_phy_ps@sam_data[["Position"]]
adonis_result=adonis(NAC_phy_ps@otu_table~class,group,distance='bray',psrmutations=999 )
adonis_result#
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
otu.pairwise.adonis=pairwise.adonis(x=NAC_phy_ps@otu_table,factors=group$class,sim.function="vegdist",
sim.method='bray',p.adjust.m = 'BH',reduce=NULL,perm=999)
otu.pairwise.adonis
prefix='ps'
#im.function	Function used to calculate the similarity matrix, one of 'daisy' or 'vegdist' default is 'vegdist'. Ignored if x is a distance matrix.
#reduce	String. Restrict comparison to pairs including these factors. If more than one factor, separate by pipss like reduce = 'setosa|versicolor'
write.csv(adonis_result[["aov.tab"]],paste0(prefix,'.permanova.csv'))
write.csv(otu.pairwise.adonis,paste0(prefix,'.Pairwisepermanova.csv'))



