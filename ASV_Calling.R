library(devtools)
library(dada2)
library(phyloseq)
library(stringr)

fnFs <- sort(list.files(path="../Amplicon Analysis/demuxed/",pattern="_R1.fastq.gz", full.names = TRUE,recursive = T))
fnRs <- sort(list.files(path="../Amplicon Analysis/demuxed/",pattern="_R2.fastq.gz", full.names = TRUE,recursive = T))
sample.names <- word(fnFs,1,sep="_R1")
sample.names <- word(sample.names,-1,sep="/")

path<-"../"
filtFs <- file.path(path, "filtered", paste0(sample.names, "_1.fastq.gz")) #Save sample names and paths and F orientation as a separate object
filtRs <- file.path(path, "filtered", paste0(sample.names, "_2.fastq.gz")) #Save sample names and paths and R orientation as a separate object
names(filtFs) <- sample.names #assign actual sample names to file path - F orientation
names(filtRs) <- sample.names #assign actual sample names to file path - R orientation

#Filter and Trim files according parameters set out below
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=FALSE) # On Windows set multithread=FALSE
#head(out) #Make sure its worked

errF <- learnErrors(filtFs, multithread=6) #Train error rate for F orientation
errR <- learnErrors(filtRs, multithread=6) #Train error rate for R orientation

plotErrors(errF, nominalQ=TRUE) #Show error rates of base to base substitutions for F orientation
plotErrors(errF, nominalQ=TRUE) #Show error rates of base to base substitutions for R orientation

exists <- file.exists(filtFs) #Save sample names of samples that passed filtering and trimming
deRepFs <- derepFastq(filtFs[exists]) #Dereplicate F samples for further processing
deRepRs <- derepFastq(filtRs[exists]) #Dereplicate R samples for further processing


#Carry out sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) #Maybe look at pooling in the future, see how it affects ASV creation?
dadaRs <- dada(filtRs, err=errR, multithread=TRUE) #Maybe look at pooling in the future, see how it affects ASV creation?
dadaFs[[1]] #Make sure it worked


mergers <- mergePairs(dadaFs, deRepFs, dadaRs, deRepRs, verbose=TRUE) #Merge F and R samples keeping in mind error rates
# Inspect the merger data.frame from the first sample
head(mergers[[1]]) #Make sure its worked


##Construct a sequence Table
seqtab <- makeSequenceTable(mergers) #Make the by row merged samples into a table.
dim(seqtab) #Check dimensions to see if the expected number of samples have made it into the table with the expected number of reads

# Inspect distribution of sequence lengths

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) #Check dimensions to see if the a good number of samples have made it into the table with a good number of reads

sum(seqtab.nochim)/sum(seqtab) #Check the proportion of reads that made it past the chimera check

saveRDS(seqtab.nochim, file=paste(path, "./seqtab.nochim.rds", sep = "")) #Save table as .rds object so we can use it later in a local version of R to create a tree


##Track reads through the pipeline
getN <- function(x) sum(getUniques(x)) #See how many unique sequences there are? Very tired, will check later.
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) #As it says on the box
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim") #Set column names
rownames(track) <- sample.names #Label rownames for each sample
head(track) #Check to see the proportion of samples that made it through the whole process.

write.csv(track, "Track_reads_through_pipeline.csv") #Save the tracked number of reads through the process

#################
## downloading DECIPHER-formatted SILVA v138 reference
download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")

## loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")

## loading DECIPHER
library(DECIPHER)
packageVersion("DECIPHER") # v2.6.0 when this was initially put together, though might be different in the binder or conda installation, that's ok!

## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=6)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)
