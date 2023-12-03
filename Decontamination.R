library(decontam)
packageVersion("decontam") # 1.1.2 when this was put together

asv_tab<-read.table("./ASVs_counts.tsv", header=T, row.names=1,
                    check.names=F, sep="\t")
asv_fasta<-Biostrings::readDNAStringSet("./ASVs.fa")

asv_seqs <- (as.vector(asv_fasta))
asv_headers <- names(asv_seqs)

for (i in 1:length(asv_headers)) {
  asv_headers[i] <- paste(">", asv_headers[i], sep="")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))

map_file <- "metadata.txt" #Define pathway for mapping file
metdat<-as.data.frame(read.delim(map_file))
names(asv_tab)<-gsub("_1.fastq.gz","",names(asv_tab))
names(asv_tab)[!names(asv_tab) %in% metdat$Sample]
asv_tab<-asv_tab[names(asv_tab) %in% metdat$Sample]

library(phyloseq)

meta_reformed<-(metdat[match(colnames(asv_tab), metdat$Sample),])
rownames(meta_reformed)<-1:nrow(meta_reformed)
meta_reformed$Sample_or_Control<-ifelse(meta_reformed$Sample_or_Control=="Sample",FALSE,TRUE)
contam_df <- isContaminant(t(asv_tab), neg=meta_reformed$Sample_or_Control,threshold = 0.1)

table(contam_df$contaminant) # identified 6 as contaminants
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

# making new taxonomy table
#asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
#write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
#            sep="\t", quote=F, col.names=NA)
