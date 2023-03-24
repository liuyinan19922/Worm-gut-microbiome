library(phyloseq)
library(ggplot2)
library(vegan)
library(gridExtra)
library(picante)
library(ggprism)
library(beepr)

setwd("D:/Liu Yi-Nan/Desktop/summary/code/debug")
ASV <- read.delim('ASV_Mat.txt', row.names=1)
sample <- read.delim('sample_data.txt', row.names=1)
Tax <- read.delim('Classification.txt', row.names=1) 

OTU_mat <- as.matrix(ASV)
tax_mat <- as.matrix(Tax)
dim(OTU_mat)

OTU_f = otu_table(OTU_mat, taxa_are_rows = FALSE)
TAX = tax_table(tax_mat)
samples = sample_data(sample)
Dim(samples)
head(OTU_f)
head(TAX)
head(samples)

NAC_phy_master <- phyloseq(OTU_f, TAX, samples) #Make the phyloseq object

set.seed(100)

#Rarefy to even depth
ps_rar <- rarefy_even_depth(NAC_phy_master, sample.size = min(sample_sums(NAC_phy_master)),
                            rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

write.csv(otu_table(ps_rar), file = "rarefied_otu.csv")
rarecurve(otu_table(ps_rar), step=50, cex=0.5)




Fac1 <- sample
head(Fac1)
OTU <- otu_table(ps_rar)#Storing the rarefied abundacne matrix
OTUt <- data.frame(OTU) # Column name should be genera and samples should be row name
head(OTUt)
dim(OTUt)
x<- row.names(OTUt)
OTUt$Sample <- x # Adding sample names as a column to merge for subsetting based on kind of plastic later

y <- row.names(Fac1)
Fac1$Sample <- y
FOR_SUB <- merge(OTUt, Fac1, by = "Sample", all = TRUE) #merge metadata and count
head(FOR_SUB)

dim(FOR_SUB)
Dim(Fac1)
Dim(OTUt)

 #Subset by kind of different plastics
pe <- subset(FOR_SUB, Kind == "HDPE")
dim(pe)
phy <- read.tree("tree_asv")

plot(phy, cex = 0.5)
phy.dist <- cophenetic(phy)# change tree to distance


rn_pe <- pe$Sample
row.names(pe) <- rn_pe # Add row names back
rn_pe_pic <- subset(pe, select = -c(Sample, Kind,Position, Time,Stage,Order)) #Remove the metadata
pe_meta <- subset(pe, select = c(Sample, Kind,Position, Time,Stage)) #Store metadata for later use
dim(rn_pe_pic)
#Filter ASV based on count
keep <- colSums(rn_pe_pic) >= 10 
rn_pe_pic_f <- rn_pe_pic[,keep]
dim(rn_pe_pic_f)
pe_work <- decostand(rn_pe_pic_f, method = "total") #Relative abundance calculation.
apply(pe_work, 1, sum)
head(pe_work)


pe.sesmntd <- ses.mntd(pe_work, phy.dist, null.model = "taxa.labels", abundance.weighted = TRUE, 
                       runs = 999)
head(pe.sesmntd)
pe.sesmntd$NTI <- -1*pe.sesmntd$mntd.obs.z
write.csv(pe.sesmntd, "pe_NTI_weighted")


    
    pe.sesmntd=read.csv('pe_NTI_weighted',row.names = 1)
  pNTIpe=ggplot(pe.sesmntd, aes(x =pe_meta$Position,y=pe.sesmntd$NTI, color = pe_meta$Position))+
    geom_boxplot(width=0.7,fatten = 1,lwd=1.5)+
    geom_jitter(width = 0.2)+
    labs(title="Pe_NTI", y = "NTI",x='')+
    scale_color_manual(values=c("#FF6C67",'#9A1207',"#FBB0AA"))+
    scale_x_discrete( limits = c('original','planktonic', "adhering"),
                      labels = c("Source",'Planktonic', "Biofilm"))+
  theme_prism(border = TRUE)+
  theme(axis.text.x=element_text(size=14,face = "bold",angle=15),legend.position="none")+
  guides(fill = "none")+
  ylim(3.2,5.2)
  pNTIpe
  ggsave(file = "D:/Liu Yi-Nan/Desktop/pNTIpe.jpg", plot =pNTIpe, dpi = 300,width = 3.2, height = 3.75 )
  
  #####Data saving and statistical testing####
  data=data.frame(pe.sesmntd$NTI,pe_meta$Position)
  sd=tapply(data[,1],data[,2],sd)
  mean=tapply(data[,1],data[,2],mean)
  
  kresult=kruskal.test(data[,1]~data[,2], data =data)
  kresult
  wresult=pairwise.wilcox.test(data[,1],data[,2], p.adjust.method = "BH")
  wresult
  prefix='pe'
  prefix01="_NTI_data.csv"
  prefix03='_NTI_test_wilcxon.csv'
  prefix1=paste0(prefix,prefix01)
  prefix3=paste0(prefix,prefix03)
  write.csv(data.frame(mean,sd),prefix1)
  write.csv(wresult$p.value,prefix3)
  
beep(sound = "mario")
