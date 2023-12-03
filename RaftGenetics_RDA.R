
library(tidyverse)

keep<-meta_reformed %>% 
  filter(Location %in% c("Munida")) 
keep<-keep[!grepl("c|C",keep$Sample),]
keep<-keep[!keep$Species %in% c("Seawater"),]
keep<-keep[!keep$Sample %in% c("S9B","S9B_1"),]
keep<-keep %>% group_by(Specific.Sample) %>% slice_sample(n=1)
keep<-keep[!keep$Species == "Control",]                  
#keep$Sample<-keep$Specific.Sample
keep_inplace<-meta_reformed %>% 
  filter(!Location %in% c("Munida")) %>%
  filter(Species == "Antarctica") 
keep_inplace$Specific.Sample[keep_inplace$Location=="StewartIsland"]<-substring(keep_inplace$Specific.Sample[keep_inplace$Location=="StewartIsland"],1,4)
keep_inplace$Group<-paste(keep_inplace$SampleType,keep_inplace$Specific.Sample)

keep_inplace <- keep_inplace %>% group_by(Group) %>% slice_sample(n=1)

ASV_physeq <- phyloseq(count_tab_phy, sample_data(meta_reformed),fastas)
ASV_physeq<-rarefy_even_depth(ASV_physeq,4000)
ASV_physeq<- prune_samples(
ASV_physeq@sam_data$Sample[!grepl(c("Wai|Wha|WRD|RRG|Kai|CC|AK|LB"),ASV_physeq@sam_data$Sample)],
ASV_physeq)
keep_details<-rbind(keep,keep_inplace)
keep_details<-keep_details[keep_details$Sample %in% ASV_physeq@sam_data$Sample,]


ASV_physeq<-prune_samples(keep_details$Sample,ASV_physeq)
ASV_physeq@sam_data$Sample


ASV_Munida_Total<-prune_samples(keep_details$Sample[keep_details$Location=="Munida"],ASV_physeq)

ASV_InPlace<-prune_samples(keep_details$Sample[!keep_details$Location=="Munida"],ASV_physeq)
ASV_Munida_Total_temp<-prune_samples(ASV_Munida_Total@sam_data$Sample[
  ASV_Munida_Total@sam_data$SampleType=="Blade"
],ASV_Munida_Total)

my_subset <- ASV_Munida_Total@sam_data$Sample[!grepl("[Cc]$", ASV_Munida_Total@sam_data$Sample)]
my_subset<-my_subset[! my_subset %in% c("MunF1c_1")]
ASV_Munida_Total_temp<-prune_samples(my_subset,ASV_Munida_Total)

ASV_Merge_Gen<-merge_phyloseq(ASV_InPlace,ASV_Munida_Total_temp)
sample_names(ASV_Merge_Gen)<-ASV_Merge_Gen@sam_data$Specific.Sample


###################
library(vcfR)
library(adegenet)
library(poppr)
library(dartR)
library(StAMPP)

eq_dat<-read.vcfR("../for_analysis_NZ.recode.vcf")
sampes_gen<-colnames(eq_dat@gt)
sampes_gen[sampes_gen=="MunF2b"]<-"MunF2"
sampes_gen[sampes_gen=="MunF8_reextract"]<-"MunF8"
sampes_gen[sampes_gen=="K5_original"]<-"K5"
sampes_gen[sampes_gen=="MunF3_c"]<-"MunF3"
sampes_gen[sampes_gen=="PPT1_backup"]<-"PPT1"
sampes_gen[sampes_gen=="PPT3_c"]<-"PPT3"
sampes_gen[sampes_gen=="PPT4_c"]<-"PPT4"
sampes_gen[sampes_gen=="PPT6_backup"]<-"PPT6"
sampes_gen[sampes_gen=="PPT7_c"]<-"PPT7"
sampes_gen[sampes_gen=="PPT8_c"]<-"PPT8"
sampes_gen[sampes_gen=="CHR1_c"]<-"CHR1"
sampes_gen[sampes_gen=="SI1"]<-"SIA1a"
sampes_gen[sampes_gen=="SI4"]<-"SIA4a"
sampes_gen[sampes_gen=="SIA5b"]<-"SIA5b"
sampes_gen[sampes_gen=="SIA6"]<-"SIA6a"
colnames(eq_dat@gt)<-sampes_gen
library(SNPfiltR)
vcfR<-hard_filter(vcfR=eq_dat, depth = 5, gq = 30)
mergedGP = merge_samples(ASV_Merge_Gen, "Specific.Sample")
sample_names(mergedGP)[!sample_names(mergedGP) %in% colnames(vcfR@gt)]
sample_names(mergedGP)[!sample_names(mergedGP) %in% colnames(vcfR@gt)]
vcfR_gen<-vcfR
vcfR_gen@gt<-cbind(vcfR@gt[,1], vcfR@gt[,colnames(vcfR@gt) %in% sample_names(mergedGP)])
colnames(vcfR_gen@gt)<-c("FORMAT",colnames(vcfR_gen@gt)[2:59])
vcfR_gen<-missing_by_snp(vcfR_gen, cutoff = .7)
vcfR_gen<-missing_by_sample(vcfR=vcfR_gen,cutoff = 0.7)
vcfR_gen<-min_mac(vcfR_gen, min.mac = 2)

snp_dat<-vcfR2genind(vcfR_gen)
snp_dat_clean<-missingno(snp_dat,type="mean")

library(poppr)
library(pegas)
library(hierfstat)

genspeclist<-data.frame(ASV_merged@sam_data$GeneticSpecies,ASV_merged@sam_data$Specific.Sample)


snp_dat_clean$pop<-as.factor(ifelse(grepl("WPT|CHR|PPT|SIA",rownames(snp_dat_clean$tab)),
                        "NonRaft","Raft"))


Gen_raft<-popsub(snp_dat_clean,"Raft")
Gen_nonraft<-popsub(snp_dat_clean,"NonRaft")

raft_pca<-prcomp(Gen_raft)
#raft_pca$x<-raft_pca$x[-25,]
raft_summ<-summary(raft_pca)
raft_summ<-as.data.frame(t(raft_summ$importance));colnames(raft_summ)<-c("SD","PropVar","CumPropVar")
raft_pcadat<-as.data.frame(raft_pca$x[,1:2])

nonraft_pca<-prcomp(Gen_nonraft)
nonraft_summ<-summary(nonraft_pca)
nonraft_summ<-as.data.frame(t(nonraft_summ$importance));colnames(nonraft_summ)<-c("SD","PropVar","CumPropVar")
nonraft_pcadat<-as.data.frame(nonraft_pca$x[,1:2])

ASV_NonRaft<-prune_samples(rownames(Gen_nonraft@tab),mergedGP)
ASV_NonRaft_com<-microbiome::transform(ASV_NonRaft,transform = "compositional")
ASV_NonRaft_com_hell<-microbiome::transform(ASV_NonRaft_com,"hellinger")

ASV_Raft<-prune_samples(rownames(Gen_raft@tab),mergedGP)
ASV_Raft_com<-microbiome::transform(ASV_Raft,transform = "compositional")
ASV_Raft_com_hell<-microbiome::transform(ASV_Raft_com,"hellinger")

sum(nonraft_summ$PropVar[1])
sum(raft_summ$PropVar[1:7])
raft_vec<-as.factor(genspeclist$ASV_merged.sam_data.GeneticSpecies[match(rownames(raft_pca$x),genspeclist$ASV_merged.sam_data.Specific.Sample)])
raft_rda<-rda((ASV_Raft_com@otu_table@.Data)~  raft_pca$x[,1:7]+raft_vec)
vegan::RsquareAdj(raft_rda)

anova(raft_rda,permutations = 1000)
nonraft_vec<-as.factor(genspeclist$ASV_merged.sam_data.GeneticSpecies[match(rownames(nonraft_pca$x),genspeclist$ASV_merged.sam_data.Specific.Sample)])
nonraft_rda<-rda((ASV_NonRaft_com@otu_table@.Data)~nonraft_pca$x[,1],nonraft_vec)
vegan::RsquareAdj(nonraft_rda)
anova(nonraft_rda,permutations = 1000)
