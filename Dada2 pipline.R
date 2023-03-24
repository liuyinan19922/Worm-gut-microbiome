#####DADA2 Pipeline######

path="d:/220412seq"
list.files(path)

install.packages("Rcpp")
library(dada2)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# #First perform quality filtering and cutting on the sequence. We use the filterAndTrim 
# command to trim the forward sequence and the reverse sequence to 240 and 160 bp respectively
# (truncLen=c(240,160)), the maximum number of N is set to 0 (maxN=0), and the maxEE parameter 
# indicates that it is allowed to expect in the reads The maximum value of the error 
# (see "A filtering strategy better than average quality": https://academic.oup.com/bioinformatics/article/31/21/3476/194979 )
# ; truncQ is the threshold for filtering the highest quality , the default is 2, 
# will remove sequences lower than this quality; rm.phix defaults to TRUE, 
# will remove the sequence compared to the reference PhiX genome, 
# which is commonly used in amplicons; compress, the default output is the result of the compressed format; 
# multithread is The default is single-threaded, you can change it to T or an integer to increase the speed,
# but make sure your computer has enough memory and computing resources.
# #If we need to shorten the filtering time, set the maxEE parameter to be smaller, 
# and want to keep most of the reads, then we need to set the maxEE parameter to be larger, 
# especially for reverse sequencing data (eg maxEE=c(2,5) ). In addition, when setting the trunclen parameters,
# it must be noted that the sequence after trimming can retain enough length (minimum 20bp) for splicing
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,190),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#The core algorithm of #DADA2 is also parameter learning, and the amount of calculation is very considerable.
# At present, the size of the data we measured is ten times or even a hundred times that of this example file.
# Faced with such a huge amount of data and the computing resources that need to be consumed, 
# the model’s The display is not suitable for our actual large amount of data, and the degree of fitting 
# can be adjusted by increasing the nbase parameter to reduce the amount of calculation.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)#from 1.16 We are now ready to apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)#Sequence splicing is modified with reference to 1.16
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)#Generate ASV table
dim(seqtab)

table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 361:418]#Remove the length that I think is impossible, the calculated length is 372:408


seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)#Chimera removal
dim(seqtab2.nochim)
sum(seqtab2.nochim)/sum(seqtab2)

getN <- function(x) sum(getUniques(x))#Statistical analysis steps above
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Sequence Species Annotation
taxa <- assignTaxonomy(seqtab2.nochim, paste0(path, "/silva_nr_v132_train_set.fa.gz"), multithread=TRUE)
taxa <- addSpecies(taxa, paste0(path, "/silva_species_assignment_v132.fa.gz"))
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


#The dicipher method cannot be used with the above
BiocManager::install("DECIPHER", version = "3.14")
library(DECIPHER)
packageVersion("DECIPHER")
##


#calculate the phelogenic tree
library("phangorn")
phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)








setwd(path)
seqtable.taxa.plus <- cbind('#seq'=rownames(taxa), t(seqtab2.nochim), taxa)
# ASV table write
write.table(seqtab2.nochim, "dada2_counts.txt", sep="\t", quote=F, row.names = T)
# ASV table export of annotated files
write.table(seqtable.taxa.plus , "dada2_counts.taxon.species.txt", sep="\t", quote=F, row.names = F)
# track save
write.table(track , "dada2_track.txt", sep="\t", quote=F, row.names = F)
# save RDdata
save(list = ls(all=TRUE), file = "dada2.RData")

library(DECIPHER); packageVersion("DECIPHER")


dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET

ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

#Evaluate accuracy
unqs.mock <- seqtab2.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

#Bonus: Handoff to phyloseq

library(phyloseq); packageVersion("phyloseq")#用于把语句写在一行里面时的断句，硬要说含义的话，那就是代表一条语句的终止，告诉R可以执行了
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

samples.out <- rownames(seqtab2.nochim)

rownames(df) <- samples.out


ps <- phyloseq(otu_table(seqtab2.nochim, taxa_are_rows=FALSE), 
               sample_data(df), 
               tax_table(taxa))
