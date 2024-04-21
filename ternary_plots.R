
ASV_InPlace<-readr::read_rds("../../Munida_Analysis/ASV_InPlace.RDS")
ASV_Munida_Kelp<-readr::read_rds("../../Munida_Analysis/ASV_Munida_Kelp.RDS")
library(phyloseq)

ASV_InPlace_temp<-prune_taxa(taxa_names(ASV_InPlace)[taxa_sums(ASV_InPlace)>0],ASV_InPlace)
ASV_Munida_Kelp_temp<-prune_taxa(taxa_names(ASV_Munida_Kelp)[taxa_sums(ASV_Munida_Kelp)>0],ASV_Munida_Kelp)
ASV_InPlace_temp<-prune_taxa(taxa_names(ASV_InPlace_temp)[taxa_sums(ASV_InPlace_temp)>0],ASV_InPlace_temp)
ASV_InPlace_temp<-prune_samples(ASV_InPlace_temp@sam_data$Sample[!ASV_InPlace_temp@sam_data$Location %in% c("Wainui","Ward","Kaikoura","CapeCampbell","LeBons","Rarangi")],
                           ASV_InPlace_temp)
ASV_InPlace_temp<-prune_taxa(taxa_names(ASV_InPlace_temp)[taxa_sums(ASV_InPlace_temp)>0],ASV_InPlace_temp)




core_inplace<-create_core_comm(ASV_InPlace_temp,4000)
#core_inplace[[1]]
core_plots<-cowplot::plot_grid(core_inplace[[1]],core_inplace[[5]],nrow=1)
ggsave(file="../../coreplots_supp.pdf",core_plots,width=6,height=3)
core_raft<-create_core_comm(ASV_Munida_Kelp_temp,4000)

library(phyloseq)
ASV_Inplacecore<-prune_taxa(core_inplace[[2]]$otu[1:14],ASV_InPlace_temp)

#core_raft<-create_core_comm(ASV_Munida_Kelp,4000)

BC_ranked<-core_inplace[[4]]
elbow <- which.max(BC_ranked$fo_diffs)

lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])))

graph<-ggplot(BC_ranked[1:200,], aes(x=factor(BC_ranked$rank[1:200], levels=BC_ranked$rank[1:200]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
  geom_vline(xintercept=lastCall, lty=3, col='blue', cex=.5) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbow+10, y=.15, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
  annotate(geom="text", x=lastCall-4, y=.08, label=paste("Last 2% increase (",lastCall,')',sep=''), color="blue")


inplace_core<-core_inplace[[2]]$otu[1:14]
munida_core<-core_raft[[2]]$otu[1:18]


InPlace_Comp <- microbiome::transform(ASV_InPlace, "compositional")
Munida_Comp <- microbiome::transform(ASV_Munida_Kelp, "compositional")



#munida_core<-core_members(Munida_Comp, detection = 0.0, prevalence = 80/100,include.lowest = F)
InPlace_Comp<-prune_taxa(rownames(InPlace_Comp@otu_table@.Data)[!taxa_sums(InPlace_Comp)==0],InPlace_Comp)
rmeans_inplace<-data.frame(RMean=rowMeans(InPlace_Comp@otu_table@.Data))
rmeans_inplace$Abundance_Class<-ifelse(rownames(rmeans_inplace) %in% inplace_core,"Core",
                                       ifelse(rmeans_inplace$RMean > 0.001,"Abundant","Rare"))
InPlace_Comp_DF<-as.data.frame(InPlace_Comp@otu_table@.Data)
InPlace_Comp_DF$Group<-rmeans_inplace$Abundance_Class
InPlaceTernarySums<-as.data.frame(t(rowsum((InPlace_Comp@otu_table@.Data), rmeans_inplace$Abundance_Class)))
InPlaceTernarySums$SamplePop<-InPlace_Comp@sam_data$Location

Munida_Comp<-prune_taxa(rownames(Munida_Comp@otu_table@.Data)[!taxa_sums(Munida_Comp)==0],Munida_Comp)
rmeans_munida<-data.frame(RMean=rowMeans(Munida_Comp@otu_table@.Data))
rmeans_munida$Abundance_Class<-ifelse(rownames(rmeans_munida) %in% inplace_core,"Core",
                                      ifelse(rmeans_munida$RMean > 0.001,"Abundant","Rare"))
munida_Comp_DF<-as.data.frame(Munida_Comp@otu_table@.Data)
munida_Comp_DF$Group<-rmeans_munida$Abundance_Class
munidaTernarySums<-as.data.frame(t(rowsum((Munida_Comp@otu_table@.Data), rmeans_munida$Abundance_Class)))
munidaTernarySums$Sample<-rownames(munidaTernarySums)
munidaTernarySums$Spec_Sample<-Munida_Comp@sam_data$Specific.Sample
#munidaTernarySums$RaftTime<-RaftTimes_mins$RaftTime[match(munidaTernarySums$Spec_Sample,RaftTimes_mins$Sample)]
#munidaTernarySums$RaftTime[is.na(munidaTernarySums$RaftTime)]<-0
#plot(munidaTernarySums$Core~munidaTernarySums$RaftTime)


InPlaceTernarySums$RaftType<-"Non-Raft"
munidaTernarySums$RaftType<-"Raft"

library(RColorBrewer)
colfunc<-colorRampPalette(brewer.pal(8, "Blues"))
yourScale <- c(colfunc(30))

Ternary_Sums<-rbind(munidaTernarySums[c(1,2,3,6)],
                    InPlaceTernarySums[c(1,2,3,5)])

#Ternary_Sums$MeanSD<-rich_preds$MeanSD[match(rownames(Ternary_Sums),rownames(rich_preds))]
#Ternary_Sums$RaftTime<-rich_preds$RaftTime[match(rownames(Ternary_Sums),rownames(rich_preds))]

colors<-c("#5DBB63","#3281a8")
library(ggtern)
#set.seed(1)
library(ggalt)
TernaryPlot.pdf<-ggtern(data=Ternary_Sums,mapping=aes(Abundant, Core, Rare,col=RaftType)) + 
  labs(title = "",
       color ="Type of Sample") +
  atomic_percent()   + #make ternary scales on atomic %
  geom_point(size=4,shape=16)  + scale_colour_manual(values=colors) +
  theme_rgbw() 
TernaryPlot.pdf
mean_terns<-as.data.frame(rbind(
colMeans(Ternary_Sums[Ternary_Sums$RaftType=="Raft",1:3]),
colMeans(Ternary_Sums[Ternary_Sums$RaftType=="Non-Raft",1:3])))
mean_terns$RaftType<-c("Raft","Non-Raft")
tern_sideplot<-ggtern(data=mean_terns,mapping=aes(Abundant, Core, Rare,col=RaftType)) + 
  labs(title = "",
       color ="Type of Sample") +
  atomic_percent()   + #make ternary scales on atomic %
  geom_point(size=4,shape=16)  + scale_colour_manual(values=colors) +
  theme_rgbw() 
ggsave(TernaryPlot.pdf,filename="TernaryPlot.pdf",height=5,width=5)
ggsave(tern_sideplot,filename="tern_sideplot.pdf",height=5,width=5)

tern_sideplot


ggsave(TernaryPlot.pdf,filename="TernaryPlot.pdf",height=5,width=5)

mean(Ternary_Sums$Abundant[Ternary_Sums$RaftType=="Raft"])
mean(Ternary_Sums$Core[Ternary_Sums$RaftType=="Raft"])
mean(Ternary_Sums$Rare[Ternary_Sums$RaftType=="Raft"])

mean(Ternary_Sums$Abundant[Ternary_Sums$RaftType=="Non-Raft"])
mean(Ternary_Sums$Core[Ternary_Sums$RaftType=="Non-Raft"])
mean(Ternary_Sums$Rare[Ternary_Sums$RaftType=="Non-Raft"])
