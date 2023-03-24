library(phyloseq)
library(vegan)
library(ape)

setwd("D:/Liu Yi-Nan/Desktop/Re__NMDS_results")#It was originally used to set the current working directory. Note that the call format is: setwd("target path"), remember to add double quotes
ASV <- read.delim('ASV_Mat.txt', row.names=1)
sample <- read.delim('metadata.txt', row.names=1)
sample$Time=factor(sample$Time)
Tax <- read.delim('Classification.txt', row.names=1) 
tax_mat <- as.matrix(Tax)
dim(OTU_mat)
OTU_f = otu_table(OTU_mat, taxa_are_rows = FALSE)
TAX = tax_table(tax_mat)
samples = sample_data(sample)

########rarefy data#######

NAC_phy_master <- phyloseq(OTU_f, TAX, samples) #Make the phyloseq object
set.seed(100)#Set the random number seed to ensure that the results can be repeated each time
pe_rar <- rarefy_even_depth(NAC_phy_master, sample.size = min(sample_sums(NAC_phy_master)),
                            rngseed = FALSE, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)#Rarefy to even depth
write.csv(otu_table(pe_rar), file = "rarefied_otu.csv")
rarecurve(otu_table(pe_rar), step=50, cex=0.5)

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
dim(Fac1)
FOR_SUB <- merge(OTUt, Fac1, by = "Sample", all = TRUE) #merge metadata and count #如果有同名的by #如果没有同名by.x by.y

phy <- read.tree("tree_asv")
class(phy)
plot(phy, cex = 0.5)
phy.dist <- cophenetic(phy)


NAC_phy_pe <- subset_samples(pe_rar, Kind=="HDPE")
head(FOR_SUB)#It is a combination of an otu form plus a taxon form
dim(FOR_SUB)#Split the composite table
pe <- subset(FOR_SUB, Kind == "HDPE") #Subset by kind of different plastics #select 是用来选择各个列名的，里面还可以嵌套多种选择的函数例如startwith endwith #subset(母集合，筛选条件)
dim(pe)
head(pe)#Before that, you need to introduce data to the level of pe
rn_pe <- pe$Sample
row.names(pe) <- rn_pe # Add row names back
rn_pe_pic <- subset(pe, select = -c(Sample, Kind,Position, Time,Stage,Order)) #Remove the metadata
pe_meta <- subset(pe, select = c(Sample, Kind,Position, Time,Stage,Order)) #Store metadata for later use
dim(rn_pe_pic)
keep <- colSums(rn_pe_pic) >= 50
rn_pe_pic_f <- rn_pe_pic[,keep]
dim(rn_pe_pic_f)
apply(rn_pe_pic_f, 1, sum)
pe_work <- decostand(rn_pe_pic_f, method = "total") #Relative abundance calculation
apply(pe_work, 1, sum)


#######iCAMP######

#install.packages("Bioconductor")
BiocManager::install("iCAMP")
BiocManager::install("NST")
BiocManager::install("Rtools")
BiocManager::install("MicEco")
BiocManager::install("devtools")
install.packages("githubinstall")
install_github("Russel88/MicEco")
install.packages("devtools")
library(iCAMP)
library(ape)
library(NST)
#library(Rtools)
library("devtools")
library(MicEco)
save.wd="D:/Liu Yi-Nan/Desktop/icamp/iCAMP_ASV_Final_1" #Directory where you save ur output files

#In continuationwith previous script. To be executed one at a time

comm <- pe_work
tree <- phy
treat <- pe_meta
env <- data.frame(row.names= pe_meta$Sample, "Position" = as.numeric(factor(pe_meta$Position)),"Stage"=as.numeric(factor(pe_meta$Stage)))
# match sample IDs in OTU table and treatment information table
sampid.check=match.name(rn.list=list(comm=comm,treat=treat))
# match OTU IDs in OTU table and tree file
#colnames(comm)[colnames(comm) == "ZOR0006"] <- "'Firmicutes bacterium ZOR0006'"
spid.check=match.name(cn.list=list(comm=comm),tree.list=list(tree=tree))
# 6 # calculate pairwise phylogenetic distance matrix.
prefix="pe20_2"  # prefix of the output file names. usually use a project ID.
rand.time=100  # randomization time, 1000 is usually enough. For example test, you may set as 100 or less to save time.
nworker=1 # nworker is thread number for parallel computing, which depends on the CPU core number of your computer.
memory.G=50 # to set the memory size as you need (but should be less than the available space in your hard disk), so that calculation of large tree will not be limited by physical memory. unit is Gb.

setwd(save.wd)
if(!file.exists("pd.desc")) 
{
  #Calculates between-species phylogenetic distance matrix from a tree, using bigmemory to deal with too large dataset.
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
  #The main parameters
  #comm: data, behavior samples, listed as species (OTU or ASV)
  #tree: Phylogenetic tree
  #pd.desc: Phylogenetic distance matrix file
  #pd.spname: phylogenetic distance matrix taxa id
  #pd.wd: file save path
  #rand: number of randomizations, default 1000
  #ds: Phylogenetic signal threshold, default 0.2
  #pd.cut: The position where the phylogenetic tree is truncated to get strict Bins. Default NA
  #phylo.rand.scale: Phylogenetic null model randomization method, default within.bin
  #taxa.rand.scale: taxonomy zero model randomization method, default across.all
  #phylo.metric: Zero model calculation method, default MPD.
  #sig.index: Null model significance test method
  #bin.size.limit: The size of the smallest bin (nmin). default 24
  #nworker: parallel operation, default 4 threads
  #Other parameters are not very important, and basically do not need to be moved.
  
  # icamp.bins: output after sorting the results;
  # icamp.boot: Perform bootstrapping analysis on the results;
  # qp.bin.js: Calculate the community construction of each bin (that is, the pairwise comparison between samples/communities), and then calculate the relative importance of community construction;
  # taxa.binphy.big: binning of the phylogenetic tree;
  # ps.bin: Phylogenetic signals inside bin;
  # change.sigindex: Quickly switch between different indicators for zero model significance test;
  # null.norm: Calculate the normality of the null model;
  # pdist.big: Phylogenetic tree interspecies phylogenetic distance matrix, using big memory to handle overly large datasets;
  # Calculate MNTD (mntdn), MPD (mpdn), NRI (NRI.p), NTI (NTI.p), RC (RC.bin.bigc, RC.pc), betaNTI (bNTIn.p, bNTI.bin.big ), betaMNTD(bmntd), betaMPD(bmpd), betaNRI(bNRI.bin.big,bNRIn.p)
  #
  
  }else{
  # if you already calculated the phylogenetic distance matrix in a previous run
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}

setwd(save.wd)
niche.dif=iCAMP::dniche(env = env,comm = comm,method = "niche.value",
                        nworker = nworker,out.dist=FALSE,bigmemo=FALSE,
                        nd.wd=save.wd)
#Calculate niche difference between species based on each environmental variable, 
#directly output the matrix or save the result matrix as big.matrix.

#Input must be numeric, changed the previous data "location, period" to numeric



# 8 # within-bin phylogenetic signal assessment.
# For real data, you may try several different settings of binning, and choose the one leading to the best within-bin phylogenetic signal.
# env is required for this step.
# 8.1 # try phylogenetic binning using current setttings.
ds = 0.2 # setting can be changed to explore the best choice
bin.size.limit = 20 # setting can be changed to explore the best choice. # here set as 5 just for the small example dataset. For real data, usually try 12 to 48.

phylobin=taxa.binphy.big(tree = tree, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = nworker)
#If the distance between other species and him is less than the phylogenetic threshold (ds)
#Minimum value for #bin size to include species
#The initial analysis is set to 24

#Phylogenetic binning for iCAMP analysis. To handle large phylogenetic tree, phylogenetic
#distance matrix should be calculated and saved using the package 'bigmemory' in advance.
#pd.spname：character vector, taxa id in the same rank as the big matrix of phylogenetic distances.


# 8.2 # test within-bin phylogenetic signal.
sp.bin=phylobin$sp.bin[,3,drop=FALSE]
#The function is to keep the attributes of the fetched data still as a data frame
#one-column matrix or data.frame, indicating the bin ID for each species (OTU or ASV), 

sp.ra=colMeans(comm/rowSums(comm))
#The calculation should be the ratio of the total asv
#one-column matrix or data.frame, or a vector with name for each element, 
#indicating mean relative abundance of each species

abcut=3 # you may remove some species, if they are too rare to perform reliable correlation test.
commc=comm[,colSums(comm)>=abcut,drop=FALSE]
dim(commc)
spname.use=colnames(commc)

#Use Mantel test to evaluate phylogenetic signal within each bin
#i.e. correlation between phylogenetic distance and niche difference.

binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "pearson",r.cut = 0.1, p.cut = 0.05, min.spn = 5)
#spname.use:to specify which species will be used for phylogenetic signal test. 
#           Default is NULL, means to use all species.
#r.cut/p.cut		
#the cutoff of correlaiton coefficient or p to identify significant correlation.

if(file.exists(paste0(prefix,".PhyloSignalSummary.csv"))){appendy=TRUE;col.namesy=FALSE}else{appendy=FALSE;col.namesy=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0(prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)
if(file.exists(paste0(prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE}else{appendy2=FALSE;col.namesy2=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),file = paste0(prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)
# since this example small data is randomly generated, the correlation should be very weak.
# usually, you are looking for a binning setting lead to higher RAsig.abj (relative abundance of bins with significant phylogenetic signal) and relative high meanR (mean correlation coefficient across bins).
# see help document of the function "ps.bin" for the meaning of output.

# 9 # iCAMP analysis
# 9.1 # without omitting small bins.
# commonly use # set sig.index as Confidence instead of SES.RC (betaNRI/NTI + RCbray)
bin.size.limit = 20 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
sig.index="Confidence" # see other options in help document of icamp.big.
icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
#Infer community assembly mechanism by phylogenetic-bin-based null model analysis
#
#sig.index:character, the index for null model significance test.



# there are quite a few parameters in this function, please check the help document of "icamp.big".
# output files:
# Test.iCAMP.detail.rda: the object "icres" saved in R data format. it is a list object. The first element bNRIiRCa is the result of relative importance of each assembly process in each pairwise comparison (each turnover). The second element "detail" including binning information (named taxabin), phylogenetic and taxonomic metrics results in each bin (named like bNRIi, RCa, etc.), relative abundance of each bin (bin.weight), relative importance of each process in each turnover between communities (processes), input settings (setting), and input community data matrix (comm). See help document of the function icamp.big for more details.


##########################################################################################
treat[4]
icamp.detail <- icres$detail
treat_1 <- treat[, -(1:2)]
# 10 # iCAMP bin level statistics
icbin=icamp.bins(icamp.detail,treat = treat_1,
                 clas=NULL,silent=TRUE, boot = TRUE,
                 rand.time = rand.time,between.group = TRUE)
#This function is to calculate various statistic index to assess relative 
#importance of each process in each bin and each turnover, and bin's contribution to each process.

save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) # just to archive the result. rda file is automatically compressed, and easy to load into R.
write.csv(icbin$Pt,file = paste0(prefix,".ProcessImportance_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptk,file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptuv,file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),row.names = FALSE)
write.csv(icbin$BPtk,file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),row.names = FALSE)

# output files:
# iCAMP.Summary.rda: the object "icbin" saved in R data format. see help document of the function icamp.bins for description of each element in the object.
# 1 Relative importance of each process in governing the turnovers in a group of samples.
# 2 Relative importance of each process in governing the turnovers of each bin among a group of samples.
# 3 Relative importance of each process in governing the turnovers between each pair of communities (samples).
# 4 Bin contribution to each process, measuring the contribution of each bin to the relative importance of each process in the assembly of a group of communities.

# 11 # Bootstrapping test
# please specify which column in the treatment information table.
i=3
treat.use=treat[,i,drop=FALSE]
icamp.result=icres$CbMPDiCBraya
icboot=iCAMP::icamp.boot(icamp.result = icamp.result,treat = treat.use,rand.time = rand.time,
                         compare = TRUE,silent = FALSE,between.group = TRUE,ST.estimation = TRUE)
save(icboot,file=paste0(prefix,".iCAMP.Boot_TIME_position",colnames(treat)[i],".rda"))
write.csv(icboot$summary,file = paste0(prefix,".iCAMP.BootSummary_TIME_position.",colnames(treat)[i],".csv"),row.names = FALSE)
write.csv(icboot$compare,file = paste0(prefix,".iCAMP.Compare_TIME_position.",colnames(treat)[i],".csv"),row.names = FALSE)

# output files:
# Test.iCAMP.Boot.Management.rda: the object "icboot" saved in R data format. see help document of the function icamp.boot for description of each element in the object.
# Test.BootSummary.Management.csv: a table to summarize bootstrapping results. see help document of the function icamp.boot for description of the output element "summary".
# Test.Compare.Management.csv: a table to summarize comparison index, effect size, and significance between each two groups. see help document of the function icamp.boot for description of the output element "compare".



# Details found at this link http://www.sthda.com/english/wiki/one-way-anova-test-in-r
# For our results can be used like following example
res.aov <- aov(pe.pd$PD ~ pe_meta$Position, data = pe_work)
summary(res.aov)
TukeyHSD(res.aov)
plot(res.aov, 1)
plot(res.aov, 2)
# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )


kruskal.test(pe.pd$PD ~ pe_meta$Position, data = pe_work)
pairwise.wilcox.test(pe.pd$PD, pe_meta$Position,
                     p.adjust.method = "BH")
