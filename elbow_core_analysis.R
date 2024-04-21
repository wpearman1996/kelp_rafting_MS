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
create_core_func_comm<-function(physeq,nReads){
  library(tidyverse)
  library(reshape2)
  library(vegan)
  library(ggsci)
  theme_set(theme_light())
  nReads=nReads
  otu<-physeq@otu_table@.Data
  otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
  otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance  
  occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
    rownames_to_column('otu')
  map<-sample_data(physeq)
  
  PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
    gather(Sample, abun, -otu) %>%
    left_join(map, by = 'Sample') %>%
    group_by(otu, Location) %>%
    summarise(plot_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
              coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
              detect=ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if detected and 0 if not detected with specific genotype
    group_by(otu) %>%
    summarise(sumF=sum(plot_freq),
              sumG=sum(coreSite),
              nS=length(Location)*2,
              Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 
  
  otu_ranked <- occ_abun %>%
    left_join(PresenceSum, by='otu') %>%
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
  
  for(i in 2:30){
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))  
  }
  
  x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  df_full <- data.frame(x_names,x)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- left_join(BCaddition,df_full, by='x_names')
  
  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  
  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    summarise(MeanBC=mean(BC)) %>%
    arrange(-desc(MeanBC)) %>%
    mutate(proportionBC=MeanBC/max(MeanBC))
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- left_join(BC_ranked, increaseDF)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]
  
  fo_difference <- function(pos){
    left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
    right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
    return(left - right)
  }
  BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
  elbow <- which.max(BC_ranked$fo_diffs)
  
  lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])))
  
  graph<-ggplot(BC_ranked[1:30,], aes(x=factor(BC_ranked$rank[1:30], levels=BC_ranked$rank[1:30]))) +
    geom_point(aes(y=proportionBC)) +
    theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
    geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
    geom_vline(xintercept=lastCall, lty=3, col='blue', cex=.5) +
    labs(x='ranked OTUs',y='Bray-Curtis similarity') +
    annotate(geom="text", x=elbow+10, y=.15, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
    annotate(geom="text", x=lastCall-4, y=.08, label=paste("Last 2% increase (",lastCall,')',sep=''), color="blue")
  

  
  
  otu_ranked<-otu_ranked
  last_call<-lastCall
  output<-list(graph,otu_ranked,lastCall,)
  output
}

core_raft<-create_core_comm(ASV_Munida_Kelp,4000)
library(phyloseq);library(ggplot2)
core_inplace<-create_core_comm(ASV_InPlace,4000)
core_macro[[4]]$MeanBC

ASV_function
core_func_poha<-create_core_func_comm(ASV_function_poha,3000)
core_func_ant<-create_core_comm(ASV_function_Ant,3000)

core_inplace[[5]]
core_raft[[1]]
core_inplace[[2]]$otu[1:10]
core_raft[[2]]$otu[1:4]
