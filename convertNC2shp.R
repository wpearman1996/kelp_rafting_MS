#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(tidync)
require(terra)
library(parallel)

pth <- args[1]

src <- tidync(pth)
filename <- tail(strsplit(pth, "/")[[1]], 1)
raftID <- substr(filename, 1, nchar(filename)-3)

status_dat <- hyper_tibble(activate(src, status)) #needed intially to subset 'good' data

#all_dat <- tibble::add_column(status_dat[status_dat[,1] >= 0,], hyper_tibble(activate(src, object_type))[status_dat[,1] >= 0,1], hyper_tibble(activate(src, lon))[status_dat[,1] >= 0,1], hyper_tibble(activate(src, lat))[status_dat[,1] >= 0,1]) #ditching superfluous data

all_dat <- tibble::add_column(status_dat[status_dat[,1] >= 0,], hyper_tibble(activate(src, object_type))[status_dat[,1] >= 0,1],
                              hyper_tibble(activate(src, lon))[status_dat[,1] >= 0,1], hyper_tibble(activate(src, lat))[status_dat[,1] >= 0,1],
                              hyper_tibble(activate(src, age_seconds))[status_dat[,1] >= 0,1]) #ditching superfluous data


rm(status_dat) #should free up some memory
 
all_dat_sub<-all_dat#[all_dat$trajectory %in% sample(all_dat$trajectory,1000),]

#polyline_shp <- lapply(unique(all_dat$trajectory), function(traj_num){
polyline_shp <- lapply(unique(all_dat_sub$trajectory), function(traj_num){

traj_dat <- all_dat_sub[all_dat_sub$trajectory == traj_num, ]
rsn <- max(tail(traj_dat[, "status"], 3)) #get status 
traj_dat <- traj_dat[1:tail(which(traj_dat[, "status"] == rsn), 1),] #drop final point if particle left roi
shp <- as.lines(vect(data.frame(traj_dat[, "lon"], traj_dat[, "lat"],traj_dat[,"age_seconds"]), crs = "+proj=longlat +datum=WGS84 +no_defs"))

shp$raft <- raftID 
shp$trajectory <- traj_num 
shp$object <- unique(traj_dat[, "object_type"])
shp$start <- head(traj_dat[, "time"], 1)
shp$end <- tail(traj_dat[, "time"], 1) #get final timestep
shp$age <- tail(traj_dat[, "age_seconds"], 1)
shp$reason <- rsn #get status 
return(shp) 
}) |> vect(crs = "+proj=longlat +datum=WGS84 +no_defs") #make a list of trajectories and combine them - if we do all at once it tries to link all the lines

writeVector(polyline_shp, paste0(raftID, ".shp"),overwrite=T) #output
#
