#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(tidync)
require(terra)
#library(parallel)
library(lubridate)
library(proj4)
library(doParallel)
shp_withage <- args[1]


ostia_sst<-rast("/nesi/nobackup/uoo03316/Rafting_Inference/ostia_sst_interped.tiff")

#for (i in 1:nlyr(ostia_sst)) {
#  layer <- ostia_sst[[i]]
#  layer_mean <- terra::focal(layer, fun=mean, na.rm=TRUE, w=matrix(1, nrow=3, ncol=3))#;plot(layer_mean)
#  layer_mean <- terra::focal(layer_mean, fun=mean, na.rm=TRUE, w=matrix(1, nrow=3, ncol=3))
#  layer_filled<-layer
#  terra::values(layer_filled) <- ifelse(is.na(terra::values(layer_filled)), terra::values(layer_mean), terra::values(layer_filled))
#  ostia_sst[[i]] <- layer_filled
#}
#terra::writeRaster(ostia_sst,"ostia_sst_interped.tiff")



extract_sst<-function(df){
  sst_raft<-sapply(1:nrow(df), function(i) {
    z<-terra::extract(x = terra::subset(ostia_sst, df$Date[i]), y = df[i,1:2],method="simple")
    z[,2]
  })
  sst_static<-sapply(1:nrow(df), function(i) {
    z<-terra::extract(x = terra::subset(ostia_sst, df$Date[i]), y = df[i,6:7],method="simple")
    z[,2]
  })
  data.frame(Raft=sst_raft,Static=sst_static)
  }
calculat_sst_fast <- function(shp) {
  munf7 <- vect(shp)
  munf7_coords <- lapply(1:(length(munf7$trajectory)), function(i) {
    coords <- crds(munf7[munf7$trajectory == munf7$trajectory[i]], df = TRUE)
    coords$CollectDate <- abs(as.numeric(munf7$age[i]))
    coords$DepartDate <- abs(as.numeric(munf7$end[i]))
    coords$CollectDate <- as.POSIXct(coords$DepartDate + coords$CollectDate, origin = "1970-01-01")
    coords$DepartDate <- as.POSIXct(coords$DepartDate, origin = "1970-01-01")
    coords$Date <- seq(coords$CollectDate[1], coords$DepartDate[1], length.out = nrow(coords))
    coords$Date<-stringr::word(coords$Date,1,sep=" ")
   # coords$Date_MD <- format(coords$Date, "%m-%d")
    #coords_transformed <- proj4::project(coords[,1:2],
    #                              "+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", degrees=T)
    #coords$x<-abs(as.vector(coords_transformed$x))
    #coords$y<-abs(as.vector(coords_transformed$y))
    coords$StaticX <- tail(coords$x, 1)
    coords$StaticY <- tail(coords$y, 1)
    coords
  })
  munf7_coords_sst <- lapply(munf7_coords, extract_sst)
  munf7_coords_sst
}


registerDoParallel(strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")))
print(strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")))
shp_withage<-vect(shp_withage)
chunks <- split(shp_withage, shp_withage$trajectory)
chunk_coords <- lapply(chunks, function(chunk) {
  coords <- crds(chunk, df = TRUE)
  coords$CollectDate <- abs(as.numeric(chunk$age[1]))
  coords$DepartDate <- abs(as.numeric(chunk$end[1]))
  coords$CollectDate <- as.POSIXct(coords$DepartDate + coords$CollectDate, origin = "1970-01-01")
  coords$DepartDate <- as.POSIXct(coords$DepartDate, origin = "1970-01-01")
  coords$Date <- seq(coords$CollectDate[1], coords$DepartDate[1], length.out = nrow(coords))
  coords$Date<-stringr::word(coords$Date,1,sep=" ")
#  coords$Date_MD <- format(coords$Date, "%m-%d")
  #coords_transformed <- proj4::project(coords[,1:2],
  #                              "+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", degrees=T)
 # coords$x<-abs(as.vector(coords_transformed$x))
  #coords$y<-abs(as.vector(coords_transformed$y))
  coords$StaticX <- tail(coords$x, 1)
  coords$StaticY <- tail(coords$y, 1)
  coords
})

print("Coords Generated, now calculating SST")
outputdir<-stringr::word(args[1],1,sep="_")
outputdir<-paste0(outputdir,"_ostia_interp")
dir.create(path=outputdir)


#sst_dats <- foreach(i = 1:length(chunk_coords), .combine = "c") %dopar% {
#    extract_sst(chunk_coords[[i]])
#}
#save(sst_dats,file=args[2])
log_file <- file(paste0(outputdir,"_log.txt"), open = "w")
already_finish_chunks<- list.files(path=outputdir,pattern="*txt")
already_done<-as.numeric(stringr::word(stringr::word(already_finish_chunks,3,sep="_"),1,sep="[.]"))

# Run the foreach loop and write output to files
sst_files <- foreach(i = 1:length(chunk_coords)) %dopar% {
  if(i %in% already_done) {
    cat(paste("Chunk", i, "is already done.\n"))
    next  # Skip iteration if i is in already_done
  }
  # Create a unique filename for this chunk
  file_name <- paste0(outputdir,"/sst_chunk_", i, ".txt")
  
  # Call extract_sst() and write output to file
  sst_data <- tryCatch(
    extract_sst(chunk_coords[[i]]),
    error = function(e) e
  )
  if (inherits(sst_data, "error")) {
    cat(paste("Error processing chunk", i, ":", sst_data$message, "\n"))
  } else {
    write.table(sst_data, file_name, row.names = FALSE)
    cat(paste("Processed chunk", i, "and saved to", file_name, "\n"))
    file_name  # Return filename for later use
  }
}

# Close the file connection
sink()
close(log_file)
files2gzip <- dir(outputdir, full.names = TRUE)

zip(paste0(outputdir,".zip"),files=files2gzip)
unlink(outputdir,recursive=T)
