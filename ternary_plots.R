create_core_comm<-function(physeq,nReads){
  library(tidyverse)
  library(reshape2)
  library(vegan)
  library(ggsci)
  library(phyloseq)
  theme_set(theme_light())
  nReads=nReads
  otu<-(physeq@otu_table@.Data)
  otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
  otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance  
  occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
    rownames_to_column('otu')
  map<-sample_data(physeq)
  
  PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
    gather(Sample, abun, -otu) %>%
    dplyr::left_join(map, by = 'Sample') %>%
    group_by(otu, Location) %>%
    summarise(plot_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
              coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
              detect=ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if detected and 0 if not detected with specific genotype
    group_by(otu) %>%
    dplyr::summarise(sumF=sum(plot_freq),
                     sumG=sum(coreSite),
                     nS=length(Location)*2,
                     Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 
  
  otu_ranked <- occ_abun %>%
    dplyr::left_join(PresenceSum, by='otu') %>%
    transmute(otu=otu,
              rank=Index) %>%
    arrange(desc(rank))
  
  BCaddition <- NULL
  
  otu_start=otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start,])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCaddition <- rbind(BCaddition,df_s)
  
  for(i in 2:200){
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i 
    BCaddition <- dplyr::left_join(BCaddition, df_a, by=c('x_names'))  
  }
  
  x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  df_full <- data.frame(x_names,x)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- dplyr::left_join(BCaddition,df_full, by='x_names')
  
  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  
  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    dplyr::summarise(MeanBC=mean(BC)) %>%
    arrange(-desc(MeanBC)) %>%
    mutate(proportionBC=MeanBC/max(MeanBC))
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- dplyr::left_join(BC_ranked, increaseDF)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]
  
  fo_difference <- function(pos){
    left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
    right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
    return(left - right)
  }
  BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
  elbow <- which.max(BC_ranked$fo_diffs)
  
  lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])))
  
  graph<-ggplot(BC_ranked[1:50,], aes(x=factor(BC_ranked$rank[1:50], levels=BC_ranked$rank[1:50]))) +
    geom_point(aes(y=proportionBC)) +
    theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
    #geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
    geom_vline(xintercept=lastCall, lty=3, col='blue', cex=.5) +
    labs(x='ranked OTUs',y='Bray-Curtis similarity') +
    #annotate(geom="text", x=elbow+10, y=.15, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
    annotate(geom="text", x=lastCall-4, y=.08, label=paste("Last 2% increase (",lastCall,')',sep=''), color="blue")
  occ_abun$fill <- 'no'
  occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:lastCall]] <- 'core'
  
  spp=t(otu)
  taxon=as.vector(rownames(otu))
  
  #Models for the whole community
  obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
  sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
  sta.np.16S <- sta.np
  
  above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
  below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness
  
  sncmplot<-ggplot() +
    geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.2)+
    geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='blue', size=1.8) +
    geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=.25) +
    geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=.25)+
    geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=.25)+
    labs(x="log10(mean relative abundance)", y="Occupancy")
  
  
  otu_ranked<-otu_ranked
  last_call<-lastCall
  output<-list(graph,otu_ranked,lastCall,BC_ranked,sncmplot)
  output
}

sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
  
  ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}
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
#ggsave(file="../../coreplots_supp.pdf",core_plots,width=6,height=3)
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
