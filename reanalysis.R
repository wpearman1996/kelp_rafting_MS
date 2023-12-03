library(ggpubr)
######################## 
library(tidyverse)
library(microbiome)
library(phyloseq)
count_tab<-read.table("../../../Microbiome_DADA2_Results/ASVs_counts-no-contam.tsv", header=T, row.names=1,
                      check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("../../../Microbiome_DADA2_Results/ASVs_taxonomy-no-contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
tax_tab_lotus<-as.matrix(read.delim("../../../Microbiome_DADA2_Results/ASVs-no-contam.fa.hier"))
tax_tab_lotus<-gsub("[?]",NA,tax_tab_lotus)
rownames(tax_tab_lotus)<-tax_tab_lotus[,1]
tax_tab_lotus<-tax_tab_lotus[,2:8]
colnames(tax_tab_lotus)<-colnames(tax_tab)
tax_tab<-tax_tab_lotus
tax_tab<-as.data.frame(tax_tab_lotus)
tax_tab$domain<-paste0("k__",tax_tab$domain)
tax_tab$phylum<-paste0("p__",tax_tab$phylum)
tax_tab$class<-paste0("c__",tax_tab$class)
tax_tab$order<-paste0("o__",tax_tab$order)
tax_tab$family<-paste0("f__",tax_tab$family)
tax_tab$genus<-paste0("g__",tax_tab$genus)
tax_tab$species<-paste0("s__",tax_tab$species)
tax_tab_taxa<-paste(tax_tab$domain,tax_tab$phylum,tax_tab$class,tax_tab$order,
                    tax_tab$family,tax_tab$genus,tax_tab$species,sep="; ")
tax_tab_taxa<-data.frame(tax_tab_taxa,rownames(tax_tab))
tax_tab2keep<-tax_tab_taxa$rownames.tax_tab.[!(grepl("Chloroplast|Mitochondria|Eukaryota",(tax_tab_taxa$tax_tab_taxa)))]
tax_tab_lotus_filt<-tax_tab_lotus[rownames(tax_tab_lotus) %in% tax_tab2keep,]


fastas<-Biostrings::readDNAStringSet("../../../Microbiome_DADA2_Results/ASVs-no-contam.fa")
library(ape)

otu_rel<-read.delim("../../../Microbiome_DADA2_Results/99_otu_taxonomy.txt",head=F)
library(castor)
cast_tree<-read.tree("../../../Microbiome_DADA2_Results/ASV_retree2_maxit_0fasttree_gtr.tree")
ASV_list<-cast_tree$tip.label[!grepl("ASV",cast_tree$tip.label)]

library(phyloseq)
metdat<-as.data.frame(read.delim("../../../Microbiome_DADA2_Results/metadata.txt"))
metdat$AgeClass<-ifelse(metdat$Species=="Antarctica",ifelse(grepl("J",metdat$Sample),"Juvenile","Adult"),NA)
meta_reformed<-metdat[metdat$Sample %in% intersect(metdat$Sample,colnames(count_tab)),]
meta_reformed$AgeClass[meta_reformed$Species=="Willana"]<-"Adult"
count_tab<-count_tab[,colnames(count_tab) %in% meta_reformed$Sample]
meta_reformed$SampleType[meta_reformed$Sample=="CCWater4"]<-"Seawater"

rownames(meta_reformed)<-meta_reformed$Sample
count_tab_phy<-otu_table(count_tab,taxa_are_rows = T)

tax_tab_phy_lotus<-tax_table(tax_tab_lotus)
tax_tab_lotus_df<-as.data.frame(tax_tab_lotus)
euk_tax<-rownames(tax_tab_lotus_df[tax_tab_lotus_df$domain=="Eukaryota",])
sum(count_tab[rownames(count_tab) %in% euk_tax,])/
  sum(count_tab)
rownames(count_tab)



ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy_lotus, sample_data(meta_reformed),fastas,cast_tree)

ASV_physeq<-prune_taxa(tax_tab2keep,ASV_physeq)

library(tidyverse)
keep<-meta_reformed %>% 
  filter(Location %in% c("Munida")) 
keep_inplace<-meta_reformed %>% 
  filter(!Location %in% c("Munida")) %>%
  filter(!Species == "ISCA") %>% 
  filter(!Species %in% c("Substrate","ExtractionNegative","Control","Willana"))

keep_inplace$Group<-paste(keep_inplace$SampleType,keep_inplace$Specific.Sample)
keep_inplace <- keep_inplace %>% group_by(Group) %>% slice_sample(n=1)
ASV_physeq<-prune_samples(c(keep$Sample,keep_inplace$Sample),ASV_physeq)
set.seed(21283)
ASV_physeq<-rarefy_even_depth(ASV_physeq,4000)
readr::write_rds(ASV_physeq,"ASV_physeq_Rafting.RDS")
ASV_Munida_Total<-prune_samples(keep$Sample,ASV_physeq)
ASV_Munida_Kelp<-prune_samples(ASV_Munida_Total@sam_data$Sample[
  ASV_Munida_Total@sam_data$Species=="Antarctica"
],ASV_Munida_Total)
ASV_InPlace<-prune_samples(keep_inplace$Sample,ASV_physeq)
merged2keep<-ASV_InPlace@sam_data[!grepl(c("Wai|WRD|Wha"),ASV_InPlace@sam_data$Sample),]
merged2keep<-merged2keep$Sample[!merged2keep$SamplingYear == 2020]
ASV_InPlace<-prune_samples(merged2keep,ASV_InPlace)
ASV_InPlace<-prune_samples(ASV_InPlace@sam_data$Sample[ASV_InPlace@sam_data$Species=="Antarctica"],ASV_InPlace)
ASV_Munida_Kelp@sam_data$TissType<-ifelse(grepl("c",ASV_Munida_Kelp@sam_data$Sample),"Palmate",
                                          ifelse(grepl("C",ASV_Munida_Kelp@sam_data$Sample),"Palmate","Blade"))

########################################
# Get in raft times and standard deviations etc
RaftTimes<-read.csv("../../Munida_Analysis/target_times_v5_complete.txt",head=F)
RaftTimes<-data.table::melt(RaftTimes)
RaftTimes<-RaftTimes[!RaftTimes$value == "NaN",]
RaftTimes$value<-abs(RaftTimes$value)
RaftTimes$Hours<-RaftTimes$value/3600
RaftTimes_mins<-aggregate(Hours ~ V1, data=RaftTimes, function(x){median(x)})
colnames(RaftTimes_mins)<-c("Sample","RaftTime")


sst_rafts<-readr::read_rds("../../SST_Raft_Data/raftsds(2).RDS")

sst_rafts_means<-lapply(sst_rafts,function(x){
  #quantile(x$RaftSD,0.1,na.rm=T)
  median(x$RaftSD,na.rm=T)
})

sst_inplace<-readr::read_rds("./sst_inplace_11days.RDS")
sst_inplace<-lapply(sst_inplace,function(x){
  sd(x,na.rm=T)
})

sst_rafts_means<-do.call("rbind",sst_rafts_means)
sd_inplace<-data.frame(Sample=ASV_InPlace@sam_data$Sample,Location=ASV_InPlace@sam_data$Location)

sst_inplace<-do.call("rbind",sst_inplace)
sd_inplace$SD<-sst_inplace[,1][match(sd_inplace$Location,rownames(sst_inplace))]
sd_inplace$Location<-NULL
SD_Values<-data.frame(Sample=rownames(sst_rafts_means),SD=sst_rafts_means[,1])
SD_Values<-rbind(SD_Values,sd_inplace)
SD_Values$Sample<-gsub("_ostia_interp.zip","",SD_Values$Sample)
sst_rafts_mean<-readr::read_rds("../../SST_Raft_Data/raftmeans.RDS")
sst_rafts_means<-lapply(sst_rafts_mean,function(x){
  #quantile(x$RaftSD,0.1,na.rm=T)
  median(x$RaftSD,na.rm=T)
  })
names(sst_rafts_means)<-gsub("_ostia_interp.zip","",names(sst_rafts_means))
sst_inplace<-readr::read_rds("../../SST_Raft_Data/inplace(1).RDS")
sst_inplace<-lapply(sst_inplace,function(x){
  #quantile(x,0.1,na.rm=T)
  median(x,na.rm=T)
  })
sst_rafts_means<-do.call("rbind",sst_rafts_means)
mean_inplace<-data.frame(Sample=ASV_InPlace@sam_data$Sample,Location=ASV_InPlace@sam_data$Location)
sst_inplace<-do.call("rbind",sst_inplace)
mean_inplace$Mean<-sst_inplace[,1][match(mean_inplace$Location,rownames(sst_inplace))]
mean_inplace$Location<-NULL
Mean_Values<-data.frame(Sample=rownames(sst_rafts_means),Mean=sst_rafts_means[,1])
Mean_Values<-rbind(Mean_Values,mean_inplace)

##############
# Beta diversity analyses

ASV_merged<-merge_phyloseq(ASV_Munida_Kelp,ASV_InPlace)

merged2keep<-ASV_merged@sam_data[!grepl(c("Wai|WRD|Wha"),ASV_merged@sam_data$Sample),]
merged2keep<-merged2keep$Sample[!merged2keep$SamplingYear == 2020]
merged2keep<-c(merged2keep,ASV_physeq@sam_data$Sample[ASV_physeq@sam_data$Species=="Seawater" &
                                                        !ASV_physeq@sam_data$SamplingYear == 2020])
ASV_merged_temp<-prune_samples(merged2keep, ASV_physeq)
set.seed(12352346)
ASV_ordination<-ordinate(ASV_merged_temp,method="NMDS",distance="bray",trymax=200,k=3)

ASV_ord<-as.data.frame(ASV_ordination$points[,1:2])
ASV_ord$Sample<-rownames(ASV_ord)
temp_df<-as.data.frame(sample_data(ASV_merged_temp))
ASV_ord<-cbind(ASV_ord[,1:2],temp_df)#(ASV_ord,temp_df,by="Sample")
ASV_ord$Species[ASV_ord$Sample=="SIA4b"] <- "Antarctica"

ASV_ord$Location<-ifelse(ASV_ord$Location =="Munida","Munida","InPlace")
ASV_ord$SampleGroup<-ifelse(ASV_ord$SampleType=="Seawater","Seawater","Kelp")
ASV_ord$SampleGroup<-paste(ASV_ord$SampleGroup, ASV_ord$Location,sep=" ")

ASV_ord$RaftTime<-RaftTimes_mins$RaftTime[match(ASV_ord$Specific.Sample,RaftTimes_mins$Sample)]
ASV_ord$RaftTime<-ifelse(is.na(ASV_ord$RaftTime),0,ASV_ord$RaftTime)
ASV_ord$Group<-ifelse(ASV_ord$Specific.Sample =="Munida","Raft","NonRaft")
library(ggplot2);library(vegan)
bdiv<-betadiver(t(ASV_merged_temp@otu_table@.Data),triangular=F,method=1)
bdiv_dist<-phyloseq::distance(ASV_merged_temp,method="bray")
bdiv_mun<-(betadisper(bdiv_dist,
                      ASV_ord$SampleGroup,type="centroid",bias.adjust=T))

anova(bdiv_mun)
boxplot(bdiv_mun)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

bdiv_mun<-data.frame(Distance=bdiv_mun$distances)
bdiv_mun$SampleGroup <-ASV_ord$SampleGroup
bdiv_mun<-bdiv_mun[bdiv_mun$SampleGroup %in% c("Kelp InPlace","Kelp Munida"),]
xy.df<-as.data.frame(t(combn(c("Kelp InPlace","Kelp Munida"), 2)))
xy.list <- split(xy.df, seq(nrow(xy.df)))
xy.list <- lapply(xy.list,unlist)
colors<-c("#5DBB63","#3281a8","#9e32a8","#a85f32")


BDivPlot<-ggplot(bdiv_mun) +
  aes(x = SampleGroup, y = Distance,fill=SampleGroup) +
  geom_violin(adjust = 1L, scale = "area") + 
  theme_classic2()+ theme(text = element_text(size=18,face="bold"),legend.position="none") + geom_signif(comparisons = xy.list, step_increase = 0.01,
                                                                                                         map_signif_level=TRUE,
                                                                                                         test="wilcox.test") +ylab("Distance to Centroid") + 
  scale_fill_manual(
    values = c(`Kelp InPlace` = "#5DBB63",
               `Kelp Munida` = "#3281a8")
  ) + stat_summary(fun.data=data_summary)
BDivPlot


ASV_ord$GeneticSpecies<-ifelse(ASV_ord$GeneticSpecies %in% c("Poha","Antarctica"),ASV_ord$GeneticSpecies,
                               ifelse(ASV_ord$SampleType=="Seawater","Seawater","Unknown"))
ASV_ord$GeneticSpecies<-ifelse(ASV_ord$Species=="Seawater","Seawater",ASV_ord$GeneticSpecies)
ASV_ord$SampleGroup<-ifelse(ASV_ord$GeneticSpecies=="Seawater","Seawater","Kelp")
ASV_ord$SampleGroup<-paste(ASV_ord$SampleGroup, ASV_ord$Location,sep=" ")


ord_plot<-ggplot(ASV_ord) +
  aes(x = MDS1, y = MDS2, fill = SampleGroup, colour = SampleGroup) +
  geom_point(aes(shape=GeneticSpecies), size = 3.45) +
  scale_colour_manual(values=colors) + theme_classic2() +
  theme(text = element_text(size=18,face="bold"),legend.position="none" ) + stat_ellipse() #+ theme_classic2()
ord_plot

plots_bdiv<-(cowplot::align_plots(ord_plot,BDivPlot))

ASV_merged_dist<-phyloseq::distance(ASV_merged,method="bray")
ASV_merged_dist<-as.matrix(ASV_merged_dist)
library(reshape)
ASV_merged_dist <- melt(ASV_merged_dist)[melt(upper.tri(ASV_merged_dist))$value,]
raft_dist<-as.matrix(dist(RaftTimes_mins$RaftTime))
ASV_merged_dist$Loc1<-metdat$Location[match(ASV_merged_dist$X1,metdat$Sample)]
ASV_merged_dist$Loc2<-metdat$Location[match(ASV_merged_dist$X2,metdat$Sample)]
ASV_merged_dist$Loc2<-ifelse(ASV_merged_dist$Loc2=="Munida","Munida","Happy")
ASV_merged_dist$Loc1<-ifelse(ASV_merged_dist$Loc1=="Munida","Munida","Happy")


ASV_merged_dist$Pair<-paste0(ASV_merged_dist$Loc1,ASV_merged_dist$Loc2)
ASV_merged_dist$Pair<-ifelse(ASV_merged_dist$Pair=="MunidaMunida",
                             "Munida",ifelse(ASV_merged_dist$Pair=="HappyHappy","Happy","Diff"))
ASV_merged_dist_mun<-ASV_merged_dist[ASV_merged_dist$Pair %in% c("Diff","Happy"),]
ASV_merged_dist_mun$SpecSample<-metdat$Specific.Sample[match(ASV_merged_dist_mun$X1,
                                                             metdat$Sample)]
ASV_merged_dist_mun$RaftDist<-RaftTimes_mins$RaftTime[match(ASV_merged_dist_mun$SpecSample,
                                                            RaftTimes_mins$Sample)]
ASV_merged_dist_mun$RaftDist[is.na(ASV_merged_dist_mun$RaftDist)]<-0


means_dists_raft_meansd<-aggregate(ASV_merged_dist_mun$value, list(ASV_merged_dist_mun$X1), FUN=mean) 
means_dists_raft_meansd$SpecSample<-ASV_merged_dist_mun$SpecSample[match(means_dists_raft_meansd$Group.1,ASV_merged_dist_mun$X1)]
means_dists_raft_meansd$RaftDist<-RaftTimes_mins$RaftTime[match(means_dists_raft_meansd$SpecSample,
                                                                RaftTimes_mins$Sample)]
means_dists_raft_meansd$MunSamp<-ifelse(means_dists_raft_meansd$SpecSample %in% metdat$Specific.Sample[metdat$Location=="Munida"],
                                        "Munida","InPlace")
means_dists_raft_meansd_mun<-means_dists_raft_meansd[means_dists_raft_meansd$MunSamp=="Munida" &
                                                       means_dists_raft_meansd$RaftDist>0,]
#means_dists_raft_meansd_mun<-means_dists_raft_meansd[means_dists_raft_meansd$MunSamp=="Munida" &
#                                                       means_dists_raft_meansd$RaftDist<2000,]
means_dists_raft_meansd<-rbind(means_dists_raft_meansd_mun,means_dists_raft_meansd[means_dists_raft_meansd$MunSamp=="InPlace",])
means_dists_raft_meansd$RaftDist[is.na(means_dists_raft_meansd$RaftDist)]<-0

mean_meansd_samples<-aggregate(SD_Values$SD~SD_Values$Sample,FUN=mean)
means_dists_raft_meansd$MeanSD<-mean_meansd_samples$`SD_Values$SD`[match(means_dists_raft_meansd$SpecSample,mean_meansd_samples$`SD_Values$Sample`)]
colnames(means_dists_raft_meansd)<-c("Sample","AvgDiss","SpecSample","RaftDist","MunSamp","MeanSD")
means_dists_raft_meansd$AverageSST<-Mean_Values$Mean[match(means_dists_raft_meansd$SpecSample,Mean_Values$Sample)]
means_dists_raft_meansd$Inplace<-ifelse(means_dists_raft_meansd$Sample %in% sample_names(ASV_InPlace),"InPlace","Raft")
means_dists_raft_meansd<-means_dists_raft_meansd[complete.cases(means_dists_raft_meansd),]
library(ggplot2)

lm_bdiv<-lm(means_dists_raft_meansd$AvgDiss~means_dists_raft_meansd$Inplace+means_dists_raft_meansd$AverageSST+
              means_dists_raft_meansd$RaftDist+means_dists_raft_meansd$MeanSD)
AIC(lm_bdiv)
lm_bdiv<-lm(means_dists_raft_meansd$AvgDiss~
              means_dists_raft_meansd$RaftDist+means_dists_raft_meansd$MeanSD)
AIC(lm_bdiv)
lm_bdiv<-lm(AvgDiss~ AverageSST+ MeanSD + RaftDist,data=means_dists_raft_meansd)
AIC(lm_bdiv)
lm_bdiv<-lm(AvgDiss~  MeanSD + RaftDist,data=means_dists_raft_meansd)
AIC(lm_bdiv)
lm_bdiv<-lm(AvgDiss~ AverageSST+ MeanSD,data=means_dists_raft_meansd)
AIC(lm_bdiv)
means_dists_raft_meansd$SpecSample
lm_bdiv<-lm(AvgDiss~  MeanSD +  Inplace,data=means_dists_raft_meansd)
AIC(lm_bdiv)
bdiv_lm_flex<-flextable::as_flextable(lm_bdiv)

flextable::save_as_docx(bdiv_lm_flex,path="bdiv_lm_flex.docx")


DHARMa::testDispersion(lm_bdiv)
DHARMa::simulateResiduals(lm_bdiv,n=1000,refit=F,plot=T)

AIC(lm_bdiv)

summary(lm_bdiv)
lm_bdiv<-lm(AvgDiss~Inplace+AverageSST+
              RaftDist+MeanSD,data=means_dists_raft_meansd)

pred<-ggeffects::ggpredict(lm_bdiv,terms = c("MeanSD [0.1:0.8,by=0.01"))
pred<-data.frame(MeanSD=pred$x,AvgDiss=pred$predicted,lowconf=pred$conf.low,highconf=pred$conf.high)


kval_beta=11
gam_beta2<-(mgcv::gam(AvgDiss~ 
                        s(MeanSD,k=kval_beta) +
                        s(AverageSST,k=kval_beta)
                      + Inplace,
                       data=means_dists_raft_meansd, 
                      method = "REML",gamma=1))
AIC(gam_beta2)

summary(gam_beta2)
gam.check(gam_beta2)
bdiv_gam_flex<-flextable::as_flextable(gam_beta2)

flextable::save_as_docx(bdiv_gam_flex,path="bdiv_gam_flex.docx")

means_dists_raft_meansd



####################
library(gdm)
envTab_gdm<-data.frame(MeanSD=means_dists_raft_meansd$MeanSD,
                       SST=means_dists_raft_meansd$AverageSST,
                       RaftDist=means_dists_raft_meansd$RaftDist/24,
                  #     Inplace=means_dists_raft_meansd$Inplace,
                       row.names = means_dists_raft_meansd$Sample)
envTab_gdm$site=rownames(envTab_gdm)
ASV_merged_gdm<-prune_samples(rownames(envTab_gdm),ASV_merged)
ASV_merged_gdm<-phyloseq::distance(ASV_merged_gdm,method="bray")
ASV_merged_gdm<-as.matrix(ASV_merged_gdm)
site <- unique(rownames(ASV_merged_gdm))
gdmDissim <- cbind(site, ASV_merged_gdm)
envTab_gdm$Lat=1;envTab_gdm$Long=1
gdmTab.dis <- formatsitepair(bioData=gdmDissim, 
                               bioFormat=3, #diss matrix 
                               XColumn="Long", 
                               YColumn="Lat", 
                               predData=envTab_gdm, 
                               siteColumn="site")
gdm.1 <- gdm(data=gdmTab.dis, geo=F)
summary(gdm.1)
gdm.varImp(gdmTab.dis,geo=F,nPerm = 1000,cores = 6)
dev.off()
varSet <- vector("list", 3)
names(varSet) <- c("SST", "MeanSD","RaftDist")
varSet$SST="SST";varSet$MeanSD="MeanSD";varSet$RaftDist="RaftDist"

gdm::gdm.partition.deviance(gdmTab.dis,varSet,partSpace=FALSE)
preds <- length(gdm.1$predictors)
predmax <- 0
splineindex <- 1
PSAMPLE=300
preddata <- rep(0,times=PSAMPLE)
    z <- .C( "GetPredictorPlotData",
             pdata = as.double(preddata),
             as.integer(PSAMPLE),
             as.double(gdm.1$coefficients[7:9]),
             as.double(gdm.1$knots[7:9]),
             as.integer(3),
             PACKAGE = "gdm" )
raft_Dist_plotdat<-data.frame(pdata=z$pdata,
                    Data=seq(from=gdm.1$knots[7],to=gdm.1$knots[9],length=300))

preddata <- rep(0,times=PSAMPLE)
z <- .C( "GetPredictorPlotData",
         pdata = as.double(preddata),
         as.integer(PSAMPLE),
         as.double(gdm.1$coefficients[4:6]),
         as.double(gdm.1$knots[4:6]),
         as.integer(3),
         PACKAGE = "gdm" )
raft_sst_plotdat<-data.frame(pdata=z$pdata,
                              Data=seq(from=gdm.1$knots[4],to=gdm.1$knots[6],length=300))
preddata <- rep(0,times=PSAMPLE)
z <- .C( "GetPredictorPlotData",
         pdata = as.double(preddata),
         as.integer(PSAMPLE),
         as.double(gdm.1$coefficients[1:3]),
         as.double(gdm.1$knots[1:3]),
         as.integer(3),
         PACKAGE = "gdm" )
raft_sd_plotdat<-data.frame(pdata=z$pdata,
                             Data=seq(from=gdm.1$knots[1],to=gdm.1$knots[3],length=300))
raft_sd_plotdat$Variable<-"SD";raft_sst_plotdat$Variable<-"SST";raft_Dist_plotdat$Variable="RaftTime"
gdm_plot_data<-rbind(raft_sd_plotdat,raft_sst_plotdat,raft_Dist_plotdat)
esquisse::esquisser(gdm_plot_data)

gdm_plot<-ggplot(gdm_plot_data) +
  aes(x = Data, y = pdata) +
  geom_line(colour = "#112446",size=1) +
  theme_minimal() +
  facet_wrap(vars(Variable), scales = "free") + theme_classic2() + scale_y_continuous(limits=c(0,1))


####################

library(ggplot2)
plot(ggeffects::ggpredict(gam_beta2,terms = c("MeanSD [0.1:0.8 by=0.001]")))
pred<-ggeffects::ggpredict(gam_beta2,terms = c("MeanSD [0.1:0.82 by=0.01]"))
pred<-data.frame(MeanSD=pred$x,AvgDiss=pred$predicted,lowconf=pred$conf.low,highconf=pred$conf.high)
pred$lowconf[pred$lowconf<=0]<-0


raft_bdiv_time<-ggplot(data=means_dists_raft_meansd,aes(x =MeanSD, y = AvgDiss))+
  geom_point(
    shape = "circle filled",
    size = 3.45,
    aes(fill =  RaftDist/24)
  ) +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  theme_minimal() +
  ylab("Average Dissimilarity to Non-Raft") + xlab("SST SD") +
  theme(text = element_text(size=18,face="bold")) +
  geom_ribbon(data = pred, aes(x = MeanSD, ymin = lowconf, ymax = highconf),
              alpha = 0.5, fill = "grey80") +
  geom_line(data = pred, aes(x = MeanSD, y = AvgDiss),
            color = "black", size = 1L, linetype = "solid",inherit.aes = F)



library(mgcv)
library(cowplot)
plots_bdiv<-(cowplot::align_plots(ord_plot,BDivPlot))
bdiv_plot<-plot_grid(plots_bdiv[[1]],plots_bdiv[[2]],
                     ncol=2,nrow=1,rel_widths = c(3,3),rel_heights = c(2,2))
bdiv_plot<-plot_grid(bdiv_plot,gdm_plot,ncol=1,rel_heights = c(3,2))
bdiv_plot

ggsave(filename = "sd_raft_bdiv_time.pdf",plot=bdiv_plot,width = 9,height = 9)

######################
#Alpha diversity plots
norm_rich<-estimate_richness(ASV_InPlace)
norm_rich_even<-evenness(ASV_InPlace)
norm_rich$Evenness<-norm_rich_even$pielou
norm_rich$Spec_Sample<-rownames(norm_rich)
norm_rich_agg<-reshape2::melt(norm_rich)
norm_rich_agg$Group<-paste0(norm_rich_agg$Spec_Sample,"_",norm_rich_agg$variable)
norm_rich_agg$Variable<-stringr::word(norm_rich_agg$Group,2,sep="_") 
norm_rich_agg$Sample<-stringr::word(norm_rich_agg$Group,1,sep="_") 
norm_rich_agg$RaftPeriod<-0
norm_rich_agg<-norm_rich_agg[norm_rich_agg$Variable %in% c("Observed","Chao1","Simpson","Shannon","Evenness"),]

mun_rich<-estimate_richness(ASV_Munida_Kelp)
mun_rich_even<-evenness(ASV_Munida_Kelp)
mun_rich$Evenness<-mun_rich_even$pielou
mun_rich$Spec_Sample<-ASV_Munida_Kelp@sam_data$Specific.Sample
library(vegan)

mun_rich$TissType<-ASV_Munida_Kelp@sam_data$TissType

mun_rich_preds<-mun_rich
mun_rich_preds$MeanSD<-SD_Values$SD[match(mun_rich_preds$Spec_Sample,SD_Values$Sample)]
library(tidyverse)
mun_rich_preds$RaftTime<-RaftTimes_mins$RaftTime[match(mun_rich_preds$Spec_Sample,RaftTimes_mins$Sample)]
mun_rich_preds$Raft<-rownames(mun_rich_preds)


inplace_dat<-sample_data(ASV_InPlace)
dist_nonraft<-data.frame(as.numeric(inplace_dat$Long),as.numeric(inplace_dat$Lat))
library(terra)
library(vegan)
library(microbiome)
norm_rich<-estimate_richness(ASV_InPlace)
norm_rich_even<-evenness(ASV_InPlace)
norm_rich$Evenness<-norm_rich_even$pielou
norm_rich$Spec_Sample<-rownames(norm_rich)

norm_rich$MeanSD<-SD_Values$SD[match(norm_rich$Spec_Sample,SD_Values$Sample)]
norm_rich$RaftTime<-0

mun_rich_preds$TissType<-NULL
mun_rich_preds$Raft<-NULL
rich_preds<-rbind(mun_rich_preds,norm_rich)
rich_preds<-rich_preds[!is.na(rich_preds$MeanSD),]
rich_preds$Inplace<-ifelse(rich_preds$RaftTime==0, "Non-Raft","Raft")
rich_preds$AverageSST<-Mean_Values$Mean[match(rich_preds$Spec_Sample,Mean_Values$Sample)]


kval=11
gam_rich1<-(mgcv::gam(Observed~ s(RaftTime,k=kval)+s(MeanSD,k=kval) + s(AverageSST, k = kval)+ 
                        Inplace,data=rich_preds, family = "nb", method = "REML"
                      ,gamma=0.6)) 
AIC(gam_rich1)
gam_rich2<-(mgcv::gam(Observed~ s(MeanSD,k=kval) + s(AverageSST, k = kval)+ Inplace,data=rich_preds, family = "nb", method = "REML",gamma=0.6)) 
AIC(gam_rich2)
gam_rich3<-(mgcv::gam(Observed~ s(RaftTime,k=kval)+s(AverageSST,k=kval)+ Inplace,data=rich_preds, 
                      family = "nb", method = "REML",gamma=0.6)) 
AIC(gam_rich3)



gam_rich4<-(mgcv::gam(Observed~ s(MeanSD,k=kval) + s(AverageSST,k=kval),data=rich_preds,
                      family = "nb", method = "REML",gamma=0.6)) 
AIC(gam_rich4)
#plot(gam_rich4)

gam.check(gam_rich4)


gam_rich5<-(mgcv::gam(Observed~ s(AverageSST,k=kval) +Inplace,data=rich_preds, family = "nb", method = "REML",gamma=0.6)) 
AIC(gam_rich5)
gam_rich6<-(mgcv::gam(Observed~ s(MeanSD,k=kval)+Inplace,data=rich_preds, family = "nb", method = "REML",gamma=0.6)) 
AIC(gam_rich6)
flextable::as_flextable(gam_rich4)
rich_flex<-flextable::as_flextable(gam_rich4)

kval=11
gam_even1<-(mgcv::gam(Evenness~ s(MeanSD,k=kval) + s(AverageSST, k = kval)+ Inplace,data=rich_preds, family = "betar", method = "REML"
                      ,gamma=0.6)) 
gam_even2<-(mgcv::gam(Evenness~ s(RaftTime,k=kval)+s(MeanSD,k=kval) + s(AverageSST, k = kval),data=rich_preds, family = "betar", method = "REML"
                      ,gamma=0.6)) 
gam_even3<-(mgcv::gam(Evenness~ 
                        s(RaftTime,k=kval)+
                        s(MeanSD,k=kval) + 
                        s(AverageSST, k = kval)
                      ,data=rich_preds, family = "betar", method = "REML"
                      ,gamma=0.6)) 

AIC(gam_even3)
summary(gam_even3)
concurvity(gam_even3,full = F)
concurvity(gam_rich4,full = F)
gam.check(gam_even3)


even_flex<-flextable::as_flextable(gam_even3)
flextable::save_as_docx(rich_flex,path="gamrich4.docx")
flextable::save_as_docx(even_flex,path="gam_even3.docx")



library(ggeffects);library(ggplot2)
pred<-ggeffects::ggpredict(gam_rich4,terms = c("MeanSD [0.13:0.82 by=0.01]"))
pred<-data.frame(MeanSD=pred$x,Observed=pred$predicted,lowconf=pred$conf.low,highconf=pred$conf.high)
pred$lowconf[pred$lowconf<=0]<-0
observed_pred<-ggplot(rich_preds) +
  geom_ribbon(data = pred, aes(x = MeanSD, ymin = lowconf, ymax = highconf),
              alpha = 0.8, fill = "grey80") +
   aes(x = MeanSD, y = Observed,fill=RaftTime/24)  +
  geom_point(
    shape = "circle filled",
    size = 4L,
    colour = "#112446"
  ) +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  theme_minimal() + xlab("SST SD") + ylim(0,400) + ylab("") + xlim(0.1,0.85) + 
 
  geom_line(data = pred, aes(x = MeanSD, y = Observed),
            color = "black", size = 1L, linetype = "solid",inherit.aes = F) + 
  theme_classic2() + theme(panel.grid = element_line(color = "black",linewidth=1.5))


observed_pred

pred<-ggeffects::ggpredict(gam_even3,terms = c("MeanSD [0.13:0.82 by=0.001]"))
pred<-data.frame(MeanSD=pred$x,Evenness=pred$predicted,lowconf=pred$conf.low,highconf=pred$conf.high)

evenness_pred<-ggplot(rich_preds) +
  geom_ribbon(data = pred, aes(x = MeanSD, ymin = lowconf, ymax = highconf),
              alpha = 0.5, fill = "grey80") +
  
  aes(x = MeanSD, y = Evenness,fill=RaftTime/24)+
  geom_point(
    shape = "circle filled",
    size = 4L,
    colour = "#112446"
  ) +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  xlab("SST SD") + ylim(0,1) + ylab("")+
  geom_line(data = pred, aes(x = MeanSD, y = Evenness),
            color = "black", size = 1L, linetype = "solid",inherit.aes = F) +
  theme_classic2() + theme(panel.grid = element_line(color = "black",linewidth=1.5))
evenness_pred



richness_box<-rich_agg %>% filter(Variable %in% c("Observed")) %>%  filter(!is.na(Inplace)) %>%
  ggplot() +  aes(x = Inplace, y = value,fill=Inplace) +  geom_violin(adjust = 1L, scale = "area") +
  ylim(0,400) +stat_summary(fun.data=data_summary)+
  theme_minimal() + xlab("") + ylab("Observed Richness") +
  scale_fill_manual(values = c(`InPlace` = "#5DBB63",`Raft` = "#3281a8")) +
  geom_signif(comparisons = list(c("Raft","InPlace")),
                                                map_signif_level=TRUE,
                                                test="wilcox.test") + theme_classic2() + theme(legend.position = "none")


evenness_box<-rich_agg %>% filter(Variable %in% c("Evenness")) %>%  filter(!is.na(Inplace)) %>%
  ggplot() +  aes(x = Inplace, y = value,fill=Inplace) +  geom_violin(adjust = 1L, scale = "area") +
  ylim(0,1) +stat_summary(fun.data=data_summary)+
  theme_minimal() + ylab("Evenness") + xlab("") +
  scale_fill_manual(values = c(`InPlace` = "#5DBB63",`Raft` = "#3281a8")) +
  geom_signif(comparisons = list(c("Raft","InPlace")),
                                                map_signif_level=TRUE,
                                                test="wilcox.test") + theme_classic2() + theme(legend.position = "none")

library(cowplot)
richness_set<-align_plots(richness_box,observed_pred)
eveness_set<-align_plots(evenness_box,evenness_pred)
richness<-plot_grid(richness_set[[1]],richness_set[[2]],nrow=1,rel_widths=c(2,4),labels=c("A","B"))
evenness<-plot_grid(eveness_set[[1]],eveness_set[[2]],nrow=1,rel_widths=c(2,4),labels=c("C","D"))
alphaplot<-align_plots(richness,evenness)
alphaplot<-plot_grid(alphaplot[[1]],alphaplot[[2]],nrow=2)
ggsave("alphaplot.pdf",alphaplot,width=8,height=8)

##########
#Disease analysis
library(iCAMP)
library(phyloseq)
library(ape)
library(readr)
#save(ASV_Munida_Kelp,file="ASV_Munida_Kelp_27Apr.RData")
#save(ASV_InPlace,file="ASV_InPlace_27Apr.RData")
icamp_broad_munida<-readr::read_rds("../../icamp_broad_munida.RDS")
icamp_broad_inplace<-readr::read_rds("../../icamp_broad_inplace.RDS")
q<-rbind(colSums(icamp_broad_munida$CbMNTDiCBraya[,3:7])/sum(colSums(icamp_broad_munida$CbMNTDiCBraya[,3:7])),
         colSums(icamp_broad_inplace$CbMNTDiCBraya[,3:7])/sum(colSums(icamp_broad_inplace$CbMNTDiCBraya[,3:7])))
rownames(q)<-c("Raft","Non-Raft")
colnames(q)<-c("Heterogeneous Selection", "Homogeneous Selection","Dispersal Limitaation","Homogenizing Dispersal","Drift")

q<-reshape2::melt(q)
colnames(q)<-c("variable","Driver","value")
library(tidyverse)
icamp_props<-ggplot(q) +
  aes(x = variable, y = value, fill = Driver, colour = Driver) +
  geom_col() +
  scale_fill_brewer(palette = "RdYlGn", direction = 1) +
  scale_color_brewer(palette = "RdYlGn", direction = 1) +
  theme_minimal() + ylab("Proportion") + 
  theme(axis.title.y = element_text(size = 13L),
        axis.title.x = element_text(size = 0),
        legend.position="bottom") 


colSums(icamp_broad_munida$CbMNTDiCBraya[,3:7])/sum(colSums(icamp_broad_munida$CbMNTDiCBraya[,3:7]))
colSums(icamp_broad_inplace$CbMNTDiCBraya[,3:7])/sum(colSums(icamp_broad_inplace$CbMNTDiCBraya[,3:7]))

round(100*(colSums(icamp_broad_munida$CbMNTDiCBraya[,3:7])/sum(colSums(icamp_broad_munida$CbMNTDiCBraya[,3:7]))),1)
round(100*(colSums(icamp_broad_inplace$CbMNTDiCBraya[,3:7])/sum(colSums(icamp_broad_inplace$CbMNTDiCBraya[,3:7]))),1)
############
icamp_broad_joint_munida<-readr::read_rds("icamp_broad_joint_munida.RDS")
colSums(icamp_broad_munida$CbMNTDiCBraya[,3:7])/sum(colSums(icamp_broad_munida$CbMNTDiCBraya[,3:7]))
colSums(icamp_broad_joint_munida$CbMNTDiCBraya[,3:7])/sum(colSums(icamp_broad_joint_munida$CbMNTDiCBraya[,3:7]))
colSums(icamp_broad_inplace$CbMNTDiCBraya[,3:7])/sum(colSums(icamp_broad_inplace$CbMNTDiCBraya[,3:7]))
munda_joint_df<-icamp_broad_joint_munida$CbMNTDiCBraya
ASV_Munida_Inplace<-prune_samples(c(sample_names(ASV_InPlace),sample_names(ASV_Munida_Kelp)),ASV_physeq)
ASV_Munida_Inplace<-prune_taxa(c(taxa_names(ASV_InPlace),taxa_names(ASV_Munida_Kelp)),ASV_Munida_Inplace)
ASV_Munida_Inplace<-prune_taxa(taxa_names(ASV_Munida_Inplace)[taxa_sums(ASV_Munida_Inplace)>0],ASV_Munida_Inplace)


munda_joint_df$Pop1<-ASV_Munida_Inplace@sam_data$Location[match(munda_joint_df$sample1,ASV_Munida_Inplace@sam_data$Sample)]
munda_joint_df$Pop2<-ASV_Munida_Inplace@sam_data$Location[match(munda_joint_df$sample2,ASV_Munida_Inplace@sam_data$Sample)]

munda_joint_df_munon<-munda_joint_df[munda_joint_df$Pop1=="Munida",]
munda_joint_df_munon<-munda_joint_df_munon[munda_joint_df_munon$Pop2=="Munida",]

munda_joint_df_nonmunon<-munda_joint_df[!munda_joint_df$Pop1=="Munida",]
munda_joint_df_nonmunon<-munda_joint_df_nonmunon[!munda_joint_df_nonmunon$Pop2=="Munida",]


munda_joint_df<-munda_joint_df[munda_joint_df$Pop1==munda_joint_df$Pop2,]

colSums(munda_joint_df_nonmunon[,3:7])/sum(colSums(munda_joint_df_nonmunon[,3:7]))
colSums(munda_joint_df_munon[,3:7])/sum(colSums(munda_joint_df_munon[,3:7]))
meta.group
icamp_broad_joint_munida$CbMNTDiCBraya
library(iCAMP)
icamp_bin_res<-icamp.bins(icamp_broad_joint_munida$detail)

dom_proc<-icamp_bin_res$Wtuvk

dom_proc_munida<-dom_proc[dom_proc$samp1 %in% sample_names(ASV_Munida_Kelp) &
                            dom_proc$samp2 %in% sample_names(ASV_Munida_Kelp),]

dom_proc_inplace<-dom_proc[dom_proc$samp1 %in% sample_names(ASV_InPlace) &
                             dom_proc$samp2 %in% sample_names(ASV_InPlace),]

library(ggplot2)
library(gridExtra)
factor_levels <- c("DR","HoS","HeS","DL")
prop_table_munida <- apply(dom_proc_munida[,4:ncol(dom_proc_munida)], 2, function(x) {
  prop <- table(factor(x, levels = factor_levels)) / length(x)
  prop
})
# Convert the list to a dataframe
prop_df_munida <- as.data.frame(t(prop_table_munida))
prop_df_munida$Bin<-rownames(prop_df_munida)
library(reshape2)
prop_df_munida<-melt(prop_df_munida)
munida_prop_plot<-ggplot(prop_df_munida) +
  aes(x = Bin, y = value, fill = variable) +
  geom_col() +
  scale_fill_hue(direction = 1) +
  theme_minimal()

prop_table_inplace <- apply(dom_proc_inplace[,4:ncol(dom_proc_inplace)], 2, function(x) {
  prop <- table(factor(x, levels = factor_levels)) / length(x)
  prop
})
# Convert the list to a dataframe
prop_df_inplace <- as.data.frame(t(prop_table_inplace))
prop_df_inplace$Bin<-rownames(prop_df_inplace)
prop_df_inplace<-melt(prop_df_inplace)
inplace_prop_plot<-ggplot(prop_df_inplace) +
  aes(x = Bin, y = value, fill = variable) +
  geom_col() +
  scale_fill_hue(direction = 1) +
  theme_minimal()


prop_df_diffs<-prop_df_munida
prop_df_diffs$InPlaceValues<-prop_df_inplace$value
prop_df_diffs$Diff<-prop_df_diffs$value-prop_df_diffs$InPlaceValues





prop_df_munida$Group<-"Munida"
prop_df_inplace$Group<-"InPlace"
prop_df_bins<-rbind(prop_df_inplace,prop_df_munida)
prop_df_bins$BinNum<-as.numeric(gsub("bin","",prop_df_bins$Bin))

prop_df_bins <- prop_df_bins %>%
  arrange(Group, BinNum)

prop_df_diffs$BinNum<-as.numeric(gsub("bin","",prop_df_diffs$Bin))
prop_df_bins$
prop_df_diffs <- prop_df_diffs %>%
  arrange(Group, BinNum)

prop_df_bins$Group



ggplot(prop_df_diffs) +
  aes(x = BinNum, y = Diff, fill = variable, label = BinNum) +
  geom_col() +
  scale_fill_hue(direction = 1) +
  geom_text(aes(y = -0.6), size = 4) +
  theme_minimal() +
  facet_wrap(vars(Group), ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

prop_df_bins_sub<-prop_df_bins[prop_df_bins$Bin %in% paste0("bin",c(64,70,39,4,8,40,71)),]
esquisse::esquisser(prop_df_bins_sub)


ggplot(prop_df_bins_sub) +
  aes(x = Bin, y = value, fill = variable, label = Bin) +
  geom_col() +
  scale_fill_hue(direction = 1) +
  geom_text(aes(y = -0.02), size = 4) +
  theme_minimal() +
  facet_wrap(vars(Group), ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(prop_df_bins_sub) +
  aes(x = Bin, y = Diff, fill = variable, label = Bin) +
  geom_col() +
  scale_fill_hue(direction = 1) +
  geom_text(aes(y = -0.6), size = 4) +
  theme_minimal() +
  facet_wrap(vars(Group), ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


prop_df_bins<-prop_df_bins[prop_df_bins$variable=="HoS",]

inplace_prop_plot
munida_prop_plot

bindetails<-icamp_broad_joint_munida$detail$taxabin$sp.bin
bindetails<-bindetails[bindetails$bin.id.new %in% c(70,39,4,8,40,71),]
bindetails<-icamp_broad_joint_munida$detail$taxabin$sp.bin
bindetails<-bindetails[bindetails$bin.id.new %in% c(22,6,15,16,41,53,64),]

bindetails$ASV<-rownames(bindetails)
taxtab_joint<-as.data.frame(tax_table(ASV_Munida_Inplace)@.Data)
taxtab_joint$ASV<-rownames(taxtab_joint)
bintaxa<-left_join(bindetails,taxtab_joint,by="ASV")
#############
#Selection on disease microbes
find_rows <- function(df, keywords) {
  rows <- c()
  for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
      if (length(grep(paste(keywords, collapse = "|"), df[i,j])) > 0) {
        rows <- c(rows, i)
        break
      }
    }
  }
  return(unique(rows))
}
rows <- find_rows(taxtab_joint, c("Aquimarina", "Colwellia", "Gaetbulibacter", "Kordia", 
                                  "Pseudoalteromonas", "Vibrio", "Alteromonas", "Pseudomonas", 
                                  "Flavobacterium", "Nautella", "Phaeobacter", "Thalassospira", 
                                  "Pleurocapsa", "Planococcus", "Bacillus", "Halomonas", "Croceitalea", 
                                  "Aeromonas", "Moraxella", "Plectonema"))
df_subset <- taxtab_joint[rows,]



plot_of_int<-prop_df_bins_sub %>%
  filter(Bin %in% "bin64") %>%
  mutate(Group = reorder(Group, value)) %>%
  mutate(variable = factor(variable, levels = c("DL", "DR", "HoS", "HeS"))) %>%
  ggplot() +
  aes(x = Group, y = value, fill = variable) +
  geom_col() +
  scale_fill_manual(
    values = c(DL = "#006AFF",
               DR = "#006AFF",
               HoS = "#CD2317",
               HeS = "#CD2317")
  ) +
  theme_minimal()

#ggsave("plot_of_int.pdf",plot_of_int)



ggplot(prop_df_bins_sub) +
  aes(x = Bin, y = value, fill = variable) +
  geom_col() +
  scale_fill_manual(
    values = c(DR = "#CD2317",
               HoS = "#31B425",
               HeS = "#20AFEC",
               DL = "#006AFF")
  ) +
  theme_minimal()



ASV_Munida_Inplace_core<-prune_taxa(rownames(df_subset),
                                    ASV_Munida_Inplace)
ASV_Munida_Inplace_core<-prune_samples(sample_names(ASV_Munida_Inplace_core)[sample_sums(ASV_Munida_Inplace_core)>0],ASV_Munida_Inplace_core)
#comm=t(ASV_Munida_Inplace_core@otu_table@.Data)
#tree=ASV_Munida_Inplace_core@phy_tree
#meta.group=data.frame(ASV_Munida_Inplace_core@sam_data$Location,row.names=sample_names(ASV_Munida_Inplace_core))
#icamp_broad_joint_munida_core<-icamp.cm(comm=t(ASV_Munida_Inplace_core@otu_table@.Data),
#                                        tree=ASV_Munida_Inplace_core@phy_tree,
#                                        meta.group=data.frame(ASV_Munida_Inplace_core@sam_data$Location,row.names=sample_names(ASV_Munida_Inplace_core)),
#                                        phylo.metric = "bMNTD",
#                                        sig.index = "Confidence",nworker = 8,
#                                        bin.size.limit=24,pd.wd="./munida_icamp_joint_core/",rand = 1000)
icamp_broad_joint_munida_core
munda_core_joint_df<-icamp_broad_joint_munida_core$CbMNTDiCBraya
munda_core_joint_df$Pop1<-ASV_Munida_Inplace@sam_data$Location[match(munda_core_joint_df$sample1,ASV_Munida_Inplace@sam_data$Sample)]
munda_core_joint_df$Pop2<-ASV_Munida_Inplace@sam_data$Location[match(munda_core_joint_df$sample2,ASV_Munida_Inplace@sam_data$Sample)]

munda_core_joint_df_munon<-munda_core_joint_df[munda_core_joint_df$Pop1=="Munida",]
munda_core_joint_df_munon<-munda_core_joint_df_munon[munda_core_joint_df_munon$Pop2=="Munida",]

munda_core_joint_df_nonmunon<-munda_core_joint_df[!munda_core_joint_df$Pop1=="Munida",]
munda_core_joint_df_nonmunon<-munda_core_joint_df_nonmunon[!munda_core_joint_df_nonmunon$Pop2=="Munida",]


munda_core_joint_df<-munda_core_joint_df[munda_core_joint_df$Pop1==munda_core_joint_df$Pop2,]

round(100*colSums(munda_core_joint_df_nonmunon[,3:7])/sum(colSums(munda_core_joint_df_nonmunon[,3:7])),1)
round(100*colSums(munda_core_joint_df_munon[,3:7])/sum(colSums(munda_core_joint_df_munon[,3:7])),1)
#############

cast_tree<-phyloseq::read_tree("../../../Microbiome_DADA2_Results//ASV_retree2_maxit_0fasttree_gtr.tree")
cast_tree<-drop.tip(cast_tree,cast_tree$tip.label[!cast_tree$tip.label %in% taxa_names(ASV_InPlace)])
ASV_InPlace@phy_tree<-phy_tree(cast_tree)

ASV_Munida_Kelp@phy_tree<-phy_tree(cast_tree)
comm=t(ASV_Munida_Kelp@otu_table@.Data)
tree=ASV_Munida_Kelp@phy_tree
meta.group=data.frame(ASV_Munida_Kelp@sam_data$Location,row.names=sample_names(ASV_Munida_Kelp))
#icamp_broad_munida<-icamp.cm(comm=t(ASV_Munida_Kelp@otu_table@.Data),
#                             tree=ASV_Munida_Kelp@phy_tree,
#                             meta.group=data.frame(ASV_Munida_Kelp@sam_data$Location,row.names=sample_names(ASV_Munida_Kelp)),
##                             phylo.metric = "bMNTD",
#                             sig.index = "SES.RC",nworker = 8,
#                             bin.size.limit=48,pd.wd="./munida_icamp/",rand = 1)

comm=t(ASV_Munida_Kelp@otu_table@.Data)
tree=ASV_Munida_Kelp@phy_tree
meta.group=data.frame(ASV_Munida_Kelp@sam_data$Location,row.names=sample_names(ASV_Munida_Kelp))
#icamp_broad_munida<-icamp.cm(comm=t(ASV_Munida_Kelp@otu_table@.Data),
#                             tree=ASV_Munida_Kelp@phy_tree,
#                             meta.group=data.frame(ASV_Munida_Kelp@sam_data$Location,row.names=sample_names(ASV_Munida_Kelp)),
#                             phylo.metric = "bMNTD",
#                             sig.index = "SES.RC",nworker = 8,
#                             bin.size.limit=48,pd.wd="./munida_icamp/",rand = 1)

icamp_broad_munida<-readr::read_rds("icamp_broad_joint_munida.RDS")
icamp_broad_inplace<-readr::read_rds("icamp_broad_inplace.RDS")



colSums(icamp_broad_munida$bNTIiRCa[,3:7])/sum(colSums(icamp_broad_munida$bNTIiRCa[,3:7]))
comm=t(ASV_InPlace@otu_table@.Data)
tree=ASV_InPlace@phy_tree
meta.group=data.frame(ASV_InPlace@sam_data$Location,row.names=sample_names(ASV_InPlace))
icamp_broad_inplace<-icamp.cm(comm=t(ASV_InPlace@otu_table@.Data),
                              tree=ASV_InPlace@phy_tree,
                              meta.group=data.frame(ASV_InPlace@sam_data$Location,row.names=sample_names(ASV_InPlace)),
                              phylo.metric = "bMNTD",
                              sig.index = "SES.RC",nworker = 8,
                              bin.size.limit=48,pd.wd="./inplace_icamp/",rand = 1000)

colSums(icamp_broad_inplace$bNTIiRCa[,3:7])/sum(colSums(icamp_broad_inplace$bNTIiRCa[,3:7]))

######

ASV_merged<-prune_samples(ASV_merged@sam_data$Sample[!ASV_merged@sam_data$Sample =="WRDW1_1"],ASV_merged)
ASV_merged_disease<-ASV_merged@tax_table@.Data

ASV_merge_high<-aggregate_taxa(ASV_merged,level="genus")
ASV_merge_high_com<-microbiome::transform(ASV_merge_high,transform = "compositional")
ASV_merge_high_com<-prune_taxa(paste0(c("Aquimarina", "Colwellia", "Gaetbulibacter", "Kordia", 
                                        "Pseudoalteromonas", "Vibrio", "Alteromonas", "Pseudomonas", 
                                        "Flavobacterium", "Cytophaga", "Nautella", "Phaeobacter", "Thalassospira", 
                                        "Pleurocapsa", "Planococcus", "Bacillus", "Halomonas", "Croceitalea", 
                                        "Aeromonas", "Moraxella", "Plectonema")),ASV_merge_high)

ASV_merge_high_com<-data.frame(MeanAbund=colMeans(ASV_merge_high_com@otu_table@.Data))
ASV_merge_high_com$Sample<-colnames(ASV_merge_high@otu_table@.Data)
ASV_merge_high_com$Spec_Sample<-metdat$Specific.Sample[match(ASV_merge_high_com$Sample,
                                                             metdat$Sample)]
ASV_merge_high_com$RaftTime<-RaftTimes_mins$RaftTime[match(ASV_merge_high_com$Spec_Sample,
                                                           RaftTimes_mins$Sample)]
ASV_merge_high_com$RaftTime[is.na(ASV_merge_high_com$RaftTime)]<-0
plot(ASV_merge_high_com$RaftTime,ASV_merge_high_com$MeanAbund)
ASV_merge_high_com$Percent<-100*(ASV_merge_high_com$MeanAbund/4000)
ASV_merge_high_com$Munida<-ifelse(ASV_merge_high_com$RaftTime==0,"InPlace","Munida")
ASV_merge_high_com$Perc<-100*(ASV_merge_high_com$MeanAbund/4000)
ASV_merge_high_com$Trans<-log10(ASV_merge_high_com$Perc+1)
wilcox.test(ASV_merge_high_com$Perc~ASV_merge_high_com$Munida)
disease_boxplot<-ggplot(ASV_merge_high_com) +
  aes(x = Munida, y = Perc,fill=Munida) +
  geom_violin(adjust=1L,scale="area") +stat_summary(fun.data=data_summary)+
  labs(
    x = " ",
    y = "Mean percentage of DAB"
  ) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 13L),legend.position="none") +
  geom_signif(comparisons = list(c("InPlace","Munida")), step_increase = 0.01,
              map_signif_level=TRUE,
              test="wilcox.test") +
  scale_y_continuous(trans = "sqrt") +
  scale_fill_manual(values = c(`InPlace` = "#5DBB63",`Munida` = "#3281a8"))
disease_boxplot
library(dysbiosisR)
ASV_merged@sam_data$Munida<-ifelse(sample_names(ASV_merged) %in% sample_names(ASV_InPlace),"Non-Raft","Raft")
ASV_merged@sam_data$disease<-ifelse(ASV_merged@sam_data$Munida=="Raft","Raft","Non-Raft")
merged_centroid_nonraft<-phyloseq::distance(ASV_merged,method="bray")
library(vegan)
ASV_merged_dist2cent<-betadiver(t(ASV_merged@otu_table@.Data),triangular=F,method=1)
ASV_merged_dist2cent<-(betadisper(ASV_merged_dist2cent,
                         ASV_merged@sam_data$Munida,type="centroid",bias.adjust=T))
bdiv_dist<-data.frame(Distance=ASV_merged_dist2cent$distances)
bdiv_dist$Raft<-ASV_merged@sam_data$Munida
write.csv(bdiv_dist,"../../../defensepres/bdiv_raft.csv")


dysbiosis_2 <- euclideanDistCentroids(ASV_merged,
                                      dist_mat = merged_centroid_nonraft,
                                      use_squared = TRUE,
                                      group_col = "disease",
                                      control_label = "Non-Raft",
                                      case_label = "Raft")

dysbiosis_2$disease <- factor(dysbiosis_2$disease, 
                              levels = c("Non-Raft", "Raft"))

roc_2 <- pROC::roc(as.factor(dysbiosis_2$disease),
                   dysbiosis_2$CentroidDist_score ,
                   #direction= ">",
                   plot=TRUE,
                   ci = TRUE,
                   auc.polygon=TRUE,
                   max.auc.polygon=TRUE,
                   print.auc=TRUE)


p3 <- plotDysbiosis(df=dysbiosis_2,
                    xvar="disease",
                    yvar="CentroidDist_score",
                    colors=c(CRC="brown3", N="steelblue"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score (Centroid)")


library(ggpubr)
xy.df<-as.data.frame(t(combn(c("Raft","NonRaft"), 2)))
xy.list <- split(xy.df, seq(nrow(xy.df)))
xy.list <- lapply(xy.list,unlist)
MunidaPlot<-ggplot(dysbiosis_2) +
  aes(x = Munida, y = CentroidDist_score,fill=Munida) +
  geom_violin(adjust = 1L, scale = "area") + 
  theme_classic2()+ theme(text = element_text(size=18,face="bold"),legend.position="none") + geom_signif(comparisons = xy.list, step_increase = 0.01,
                                                                                                         map_signif_level=TRUE,
                                                                                                         test="wilcox.test") +ylab("Dysbiosis Score") + 
  scale_fill_manual(
    values = c(`Non-Raft` = "#5DBB63",
               `Raft` = "#3281a8")
  ) + stat_summary(fun.data=data_summary)
MunidaPlot

dysbiosis_2$RaftTime<-rich_preds$RaftTime[match(dysbiosis_2$Specific.Sample,rich_preds$Spec_Sample)]
dysbiosis_2$SST<-rich_preds$AverageSST[match(dysbiosis_2$Specific.Sample,rich_preds$Spec_Sample)]
dysbiosis_2$SD<-rich_preds$MeanSD[match(dysbiosis_2$Specific.Sample,rich_preds$Spec_Sample)]

dysbiosis_2_nona<-dysbiosis_2[!is.na(dysbiosis_2$SD),]
dysbiosis_timeplot<-ggplot(dysbiosis_2_nona) +
  aes(x = SD, y = CentroidDist_score) +
  geom_point(
    shape = "circle filled",
    size = 3.45,
    aes(fill = RaftTime/24)
  ) +
  #scale_fill_gradient(low = "#FF0053", high = "#00B4FF") +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  geom_smooth(method="lm") +
  theme_minimal() + ylab("Dysbiosis Sore") + xlab("Standard Deviation in SST")
dysbiosis_timeplot
dis_icamp<-(cowplot::align_plots(disease_boxplot,icamp_props,dysbiosis_timeplot))
dis_icamp<-plot_grid(dis_icamp[[1]],dis_icamp[[2]],dis_icamp[[3]],ncol=3,rel_widths=c(2,2,3.5))##,
ggsave("dis_icamp.pdf",dis_icamp,width = 9,height=3)


################ ########## SVP Violin plots
#SVP Validation code
library(readr)
chunks<-list.files("../../SVP_Validation_nocollisions/RDS_Zipped/",pattern="*chunk_coords_sst.RDS",full.names = T)
sst_dat<-list.files("../../SVP_Validation_nocollisions/RDS_Zipped/",pattern="*sst_dats.RDS",full.names = T)
median(RaftTimes_mins$RaftTime)/24

sd_means_all<-list()
for(i in 1:21){
  print(chunks[i])
  svp_dat<-readr::read_rds(chunks[i])
  sst_dat_svp<-readr::read_rds(sst_dat[i])
  sd_means<-list()
  for(y in 1:length(sst_dat_svp)){
    svp_dat[[y]]$SST<-sst_dat_svp[[y]]
    dates_of_int<-as.Date(stringr::word(svp_dat[[y]]$CollectDate[1],1,sep=" "))
    date_vector <- seq(as.POSIXct(stringr::word((svp_dat[[y]]$CollectDate[1]-as.difftime(11, unit="days")),1,sep=" ")),svp_dat[[y]]$CollectDate[1],by = "days")
    svp_sub_int<-svp_dat[[y]][svp_dat[[y]]$Date %in% stringr::word(date_vector,1,sep=" "),]
    svp_sub_int$SST<-svp_sub_int$SST-273.15
    sd_means[[y]]<-data.frame(SVP_SD=sd(svp_sub_int$SST,na.rm = T),Mean_SST=mean(svp_sub_int$SST,na.rm=T))
  }
  sd_means<-do.call("rbind",sd_means)
#  svp_sub_int <- svp_sub_int %>%
#    arrange(Date) %>%     # Sort the data frame by the date column
#    group_by(Date) %>%    # Group by the date column
#    slice(c(1, n()))            # Retain the first and last instance in each group
#  sd_means[[i]]<-data.frame(SVP_SD=sd(svp_sub_int$SST,na.rm = T),Mean_SST=mean(svp_sub_int$SST,na.rm=T))
  sd_means_all[[i]]<-sd_means
}

inferred_svp_data_11day<-do.call("rbind",lapply(sd_means_all,function(x){colMeans(x,na.rm=T)}))#do.call("rbind",sd_means)
inferred_svp_data_11day<-as.data.frame(inferred_svp_data_11day)
inferred_svp_data_11day$Drifter<-as.character(tab4grant$drifter)
inferred_svp_data_11day$Type="Inferred"


drifters_actual<-read.csv("../../SVP_Validation_nocollisions/sub_drifters.csv")
drifters_actual<-split(drifters_actual,drifters_actual$drifter)
actual_sd_means<-list()
for(i in 1:21){
  dates_int_sub<-seq(as.Date(drifters_actual[[i]]$passage_date[1])-as.difftime(11,units="days"),
                     as.Date(drifters_actual[[i]]$passage_date[1]),by="days")
  svp_actual_sub_int<-drifters_actual[[i]][drifters_actual[[i]]$date %in%as.character(dates_int_sub),]
  # svp_actual_sub_int <- svp_actual_sub_int %>%
  #    arrange(date) %>%     # Sort the data frame by the date column
  #    group_by(date) %>%    # Group by the date column
  #    slice(c(1, n()))            # Retain the first and last instance in each group
  
  actual_sd_means[[i]]<-data.frame(SVP_SD=sd(svp_actual_sub_int$SST,na.rm = T),Mean_SST=mean(svp_actual_sub_int$SST,na.rm=T))
}
actual_svp_data11<-do.call("rbind",actual_sd_means)
actual_svp_data11$Drifter<-names(actual_sd_means)

tab4grant<-read.csv("../../SVP_Validation_nocollisions/tabforgrant.csv")

plot(inferred_svp_data_11day$SVP_SD~actual_svp_data11$SVP_SD)
plot(inferred_svp_data_11day$Mean_SST~actual_svp_data11$Mean_SST)

summary(lm(actual_svp_data11$SVP_SD~inferred_svp_data_11day$SVP_SD))
summary(lm(actual_svp_data11$Mean_SST~inferred_svp_data_11day$Mean_SST))

actual_svp_data11$Type="Actual"
inferred_svp_data_11day$Drifter<-NULL


svp_valid<-rbind(actual_svp_data11,inferred_svp_data_11day)

svp_valid_temp<-svp_valid 
#esquisse::esquisser(svp_valid)
t.test(svp_valid$SVP_SD~svp_valid$Type)
library(TOSTER)
res1<-t_TOST(formula=SVP_SD~Type,
             data=svp_valid,
             eqb=0.25,
             smd_ci="z",
             eqbound_type = "raw",
             hypothesis = "EQU",
             #rm_correction = T,
             paired=T)
describe(res1)

res3a = simple_htest(
 formula=SVP_SD~Type,
 data=svp_valid,
  paired = TRUE,
  mu = 0.3,
  alternative = "e"
)
describe_htest(res3a)
library(tidyverse)
library(ggplot2)
library(ggpubr)
svp_vilplot<-svp_valid %>%
  #filter(!is.na(Raft)) %>%
  ggplot() +
  aes(x = Type, y = SVP_SD,fill=Type) +
  theme_minimal() + geom_violin(adjust = 1L, scale = "area") +  
  theme_classic2()+ theme(text = element_text(size=16,face="bold"),legend.position="none") +
  geom_signif(comparisons = c("Actual","Inferred"), step_increase = 0.01,
              map_signif_level=TRUE,test = "t.test") +ylab("SST SD") +
   stat_summary(fun.data=data_summary) +
  scale_fill_manual(
    values = c(`Actual` = "#5DBB63",
               `Inferred` = "#3281a8")) 
svp_vilplot
ggsave(filename = "svp_violins.pdf",plot=svp_vilplot,width = 9,height = 4.5)

library(cowplot)
svp_violins<-(cowplot::align_plots(svp_vilplot,drift_SD_vilplot))
svp_violins<-plot_grid(svp_violins[[1]],svp_violins[[2]],
                       ncol=2,nrow=1,rel_widths = c(3,3),rel_heights = c(2,2))

svp_violins
ggsave(filename = "svp_violins.pdf",plot=svp_violins,width = 9,height = 4.5)

#esquisse::esquisser(svp_valid)
library(phyloseq)
inplace_core<-create_core_comm(ASV_InPlace,4000)
inplace_core[[1]]



