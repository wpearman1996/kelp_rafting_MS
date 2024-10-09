

library(phyloseq)
inplace_samples <- sample_names(ASV_InPlace)
raft_samples <- sample_names(ASV_Munida_Kelp)

# Create an empty vector to store the means
means <- numeric(length(raft_samples))
for (i in 1:length(raft_samples)) {
  # Get the ASVs for the current sample in Raft
  raft_asvs <- rownames(otu_table(ASV_Munida_Kelp)[,i][otu_table(ASV_Munida_Kelp)[,i]>0])
  
  # Calculate the shared ASVs between the current sample in Raft and all samples in InPlace
  shared_asvs <- sapply(inplace_samples, function(x) length(intersect(raft_asvs, rownames(otu_table(ASV_InPlace)[,x][otu_table(ASV_InPlace)[,x]>0]))))
  # Calculate the mean number of shared ASVs
  means[i] <- mean(shared_asvs)
}



means_shared_raft<-data.frame(MeanShared=means,Name=raft_samples,SpecSample=ASV_Munida_Kelp@sam_data$Specific.Sample )
means_shared_raft$RaftTime<-RaftTimes_mins$RaftTime[match(means_shared_raft$SpecSample,RaftTimes_mins$Sample)]
means_shared_raft$SD<-rich_preds$MeanSD[match(means_shared_raft$SpecSample,rich_preds$Spec_Sample)]
means_shared_raft$SST<-rich_preds$AverageSST[match(means_shared_raft$SpecSample,rich_preds$Spec_Sample)]
shared_asvs<-ggplot(means_shared_raft[!is.na(means_shared_raft$RaftTime),]) +
  aes(x = RaftTime/24, y = MeanShared) +
  geom_point(shape = "circle", size = 4L, colour = "#112446") +
  geom_smooth(span = 1L) +
  theme_minimal() + xlab("Days Spent Rafting") + ylab ("Mean number of shared ASVs with non-rafts")
