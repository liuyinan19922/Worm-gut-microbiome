#community rarify

library(phyloseq)
library(ggplot2)
library(ggprism)
library(vegan)
library(gridExtra)#gridExtra: Miscellaneous Functions for "Grid" Graphics
library(picante)#integrating Phylogenies and Ecology
library(tidyverse)
setwd("D:/Liu Yi-Nan/Desktop/Re__NMDS_results")

ASV =read.delim('ASV_Mat.txt', row.names=1)
sample =read.delim('metadata.txt', row.names=1)
sample$Time=factor(sample$Time)
Tax =read.delim('Classification.txt', row.names=1) 
OTU_mat <- as.matrix(ASV)
tax_mat <- as.matrix(Tax)
OTU_f = otu_table(OTU_mat, taxa_are_rows = FALSE)
TAX = tax_table(tax_mat)
samples = sample_data(sample)


NAC_phy_master <- phyloseq(OTU_f, TAX, samples)
set.seed(100)
rar <- rarefy_even_depth(NAC_phy_master, sample.size = min(sample_sums(NAC_phy_master)),
                         rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)


NAC_phy_pe <- subset_samples(rar, Kind=="HDPE")


#calculate alpha diversity index
plotrich.pe <- plot_richness(NAC_phy_pe, measures=c("Observed","Chao1", "Shannon",'Simpson')) 
data.pe=plotrich.pe$data

#drawing
shannon_pe=subset(data.pe,(variable=="Shannon") )
ColorContainGrey=c('#8997ab',"#ff6c67","#00bd0f","#4c9eff")




p_pe=ggplot(shannon_pe, aes(x =Position, y =value,color=Position)) + geom_boxplot(width=0.7,fatten = 1,lwd=1.5)+ geom_jitter(width = 0.2)+
  labs(x='',y="Shannon index")+
  ggtitle("Shannon pe")+
  scale_color_manual(values=c("#FF6C67",'#9A1207',"#FBB0AA"))+
  scale_x_discrete(limits = c('original','planktonic', "adhering"),breaks = c('original','planktonic', "adhering") ,labels = c("Source",'Planktonic', "Biofilm"))+
  scale_y_continuous(expand = c(0,0),limits = c(2.8,4.5))+
  theme_prism(border = TRUE)+
  theme(axis.text.x=element_text(size=14,face = "bold",angle=15),legend.position="none")+
  guides(fill = "none")
p_pe

ggsave(file = "D:/Liu Yi-Nan/Desktop/alpha_pe.jpg", plot =p_pe, dpi = 300,width = 3.4, height = 3.75 )


#####Data saving and statistical testing####
data=shannon_pe[c('value','Position')]
sd=tapply(data[,1],data[,2],sd)
mean=tapply(data[,1],data[,2],mean)


kresult=kruskal.test(data[,1]~data[,2], data =data)
kresult
wresult=pairwise.wilcox.test(data[,1],data[,2], p.adjust.method = "BH")
wresult
prefix='pe'
prefix01="_alpha_data.csv"
prefix03='_alpha_test_wilcxon.csv'
prefix1=paste0(prefix,prefix01)
prefix3=paste0(prefix,prefix03)
write.csv(data.frame(mean,sd),prefix1)
write.csv(wresult$p.value,prefix3)




