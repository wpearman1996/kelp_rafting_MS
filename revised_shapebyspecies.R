pacman::p_load(terra, ggplot2, tidyterra, ape, ggtree, patchwork,raster)
library(patchwork)
library(terra)
library(ggplot2)
library(tidyterra)
library(ape)
library(ggtree)
setwd("~/Dropbox/PhD Project/Kelp_Rafting_MS/Figure_1_creationGAD/Fig1/")

NZGD2000 <- "+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
NZ_projext <- ext(1000000, xmax = 1800000, ymin = 4600000, ymax = 5600000)
NZ_shp <- crop(project(vect("data/NZL_adm0.shp"), NZGD2000), NZ_projext)


kelp_sites <- project(vect(read.csv("data/kelp_coords.csv", row.names = 1), geom = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs"), NZGD2000)  #as.numeric(coords[coords$site_name == site_name, c("lon", "lat")] ##PROJECT COORDS

#placeholder - SST range
sigmaSST <- crop(project(rast("../ostia_sd.tiff.asc"), NZGD2000), NZ_projext)


##raster to poly contour function (terra -> raster -> terra; at least for now)
# Modified (years ago!) from: https://stackoverflow.com/questions/14379828/how-does-one-turn-contour-lines-into-filled-contour
raster2contourPolys <- function(r, levels = NULL) {
  require(maptools)
  require(rgeos)
  require(raster)
  
  ## set-up levels
  r <- raster(r)
  levels <- sort(levels)
  plevels <- c(min(values(r), na.rm=TRUE), levels, max(values(r), na.rm=TRUE)) # pad with raster range
  llevels <- paste(plevels[-length(plevels)], plevels[-1], sep=" - ")  
  llevels[1] <- paste("<", min(levels))
  llevels[length(llevels)] <- paste(">", max(levels))
  
  ## convert raster object to matrix so it can be fed into contourLines
  xmin <- extent(r)@xmin
  xmax <- extent(r)@xmax
  ymin <- extent(r)@ymin
  ymax <- extent(r)@ymax
  rx <- seq(xmin, xmax, length.out=ncol(r))
  ry <- seq(ymin, ymax, length.out=nrow(r))
  rz <- t(as.matrix(r))
  rz <- rz[,ncol(rz):1] # reshape
  
  ## get contour lines and convert to SpatialLinesDataFrame
  cat("Converting to contour lines...\n")
  cl <- contourLines(rx,ry,rz,levels=levels) 
  cl <- ContourLines2SLDF(cl)
  
  ## extract coordinates to generate overall boundary polygon
  xy <- coordinates(r)[which(!is.na(values(r))),]
  i <- chull(xy)
  b <- xy[c(i,i[1]),]
  b <- SpatialPolygons(list(Polygons(list(Polygon(b, hole = FALSE)), "1")))
  
  ## add buffer around lines and cut boundary polygon
  cat("Converting contour lines to polygons...\n")
  bcl <- gBuffer(cl, width = 0.0001) # add small buffer so it cuts bounding poly
  cp <- gDifference(b, bcl)
  
  ## restructure and make polygon number the ID
  polys <- list() 
  for(j in seq_along(cp@polygons[[1]]@Polygons)) {
    polys[[j]] <- Polygons(list(cp@polygons[[1]]@Polygons[[j]]),j)
  }
  cp <- SpatialPolygons(polys)
  cp <- SpatialPolygonsDataFrame(cp, data.frame(id=seq_along(cp)))
  
  ## cut the raster by levels
  rc <- cut(r, breaks=plevels)
  
  ## loop through each polygon, create internal buffer, select points and define overlap with raster
  cat("Adding attributes to polygons...\n")
  l <- character(length(cp))
  for(j in seq_along(cp)) {
    p <- cp[cp$id==j,] 
    bp <- gBuffer(p, width = -max(res(r))) # use a negative buffer to obtain internal points
    if(!is.null(bp)) {
      xy <- SpatialPoints(coordinates(bp@polygons[[1]]@Polygons[[1]]))[1]
      l[j] <- llevels[extract(rc,xy)]
    } 
    else { 
      xy <- coordinates(gCentroid(p)) # buffer will not be calculated for smaller polygons, so grab centroid
      l[j] <- llevels[extract(rc,xy)]
    } 
  }
  
  ## assign level to each polygon
  cp$level <- factor(l, levels=llevels)
  cp$min <- plevels[-length(plevels)][cp$level]
  cp$max <- plevels[-1][cp$level]  
  cp <- cp[!is.na(cp$level),] # discard small polygons that did not capture a raster point
  df <- unique(cp@data[,c("level","min","max")]) # to be used after holes are defined
  df <- df[order(df$min),]
  row.names(df) <- df$level
  llevels <- df$level
  
  ## define depressions in higher levels (ie holes)
  cat("Defining holes...\n")
  spolys <- list()
  p <- cp[cp$level==llevels[1],] # add deepest layer
  p <- gUnaryUnion(p)
  spolys[[1]] <- Polygons(p@polygons[[1]]@Polygons, ID=llevels[1])
  for(i in seq(length(llevels)-1)) {
    p1 <- cp[cp$level==llevels[i+1],] # upper layer
    p2 <- cp[cp$level==llevels[i],] # lower layer
    x <- numeric(length(p2)) # grab one point from each of the deeper polygons
    y <- numeric(length(p2))
    id <- numeric(length(p2))
    for(j in seq_along(p2)) {
      xy <- coordinates(p2@polygons[[j]]@Polygons[[1]])[1,]
      x[j] <- xy[1]; y[j] <- xy[2]
      id[j] <- as.numeric(p2@polygons[[j]]@ID)
    }
    xy <- SpatialPointsDataFrame(cbind(x,y), data.frame(id=id))
    holes <- over(xy, p1)$id
    holes <- xy$id[which(!is.na(holes))]
    if(length(holes)>0) {
      p2 <- p2[p2$id %in% holes,] # keep the polygons over the shallower polygon
      p1 <- gUnaryUnion(p1) # simplify each group of polygons
      p2 <- gUnaryUnion(p2)
      p <- gDifference(p1, p2) # cut holes in p1      
    } else { p <- gUnaryUnion(p1) }
    spolys[[i+1]] <- Polygons(p@polygons[[1]]@Polygons, ID=llevels[i+1]) # add level 
  }
  cp <- SpatialPolygons(spolys, pO=seq_along(llevels), proj4string=CRS(proj4string(r))) # compile into final object
  cp <- SpatialPolygonsDataFrame(cp, df)
  cat("Done!")
  terra::vect(cp)
  
}


tol_sunset_pal <- c("#364B9A", "#5385BC", "#83B8D7", "#B7DDEB", "#EAECCC", "#FDD081", "#F99858", "#E34D34", "#A50026")
grey_pal  <- gray.colors(5)

#buffer to fill NAs around coastline - for plotting only not to be used for analysis 
sigmaSST_buffered <- focal(sigmaSST, w=3, fun="median", na.policy="only",na.rm=T) 
sigmaSST_buffered[is.na(sigmaSST_buffered)] <- -999

SST_levels <- c(-999.01,0,1,2,3,4)
SST_contour <- raster2contourPolys(sigmaSST_buffered, levels = SST_levels)
SST_contour$cols <- c("white", grey_pal[1:nrow(SST_contour)-1])


#trajectory data




#target zone
#crop coast using each polygo - then buffer and simplify?
target_shp <- vect(lapply(list.files("data/shapefiles_2km_buffer/", full.names = TRUE, pattern = "shp"), function(x){
  shp_crop <- crop(as.lines(NZ_shp), project(vect(x), NZGD2000))
  return(buffer(simplifyGeom(buffer(shp_crop, 1000), tol = 50), 5000))}))

target_shp$target_name <- list.files("data/shapefiles_2km_buffer/", full.names = FALSE, pattern = "shp")

grat <- project(vect(sf::st_graticule(lon=seq(160, 180, 2), lat = seq(-50, -30, 2), ndiscr = 5000)), NZGD2000) #graticle
# phylo
popmap<-read.delim("../popmap_allsamples(1).tsv",head=F)
raft_data<-read.csv("raft_data.csv")
raft_data$Pop<-popmap$V2[match(raft_data$sample_name,popmap$V1)]
raft_data$Munida<-ifelse(raft_data$Pop=="Munida","Raft","Non-Raft")

tree_kelp <- read.tree("../../all_samples.contree")
library(tidyverse)
tips2keep<-unique(c(tree_kelp$tip.label[grepl("Tau",tree_kelp$tip.label)],raft_data$sample_name))
tree_kelp<-drop.tip(tree_kelp,tree_kelp$tip.label[!tree_kelp$tip.label %in% tips2keep])
tree_kelp$edge.length<-sqrt(tree_kelp$edge.length)
phylo_kelp <- ggtree(tree_kelp, layout="equal_angle") + geom_tiplab() 

tree_kelp$tip.label[grepl("Tau",tree_kelp$tip.label)]
raft_data_poha<-raft_data[raft_data$Species=="Poha",]
a<-MRCA(phylo_kelp, raft_data_poha$sample_name[raft_data_poha$Pop %in% c("All Day Bay","Karitane","Taiaroa_Heads")])
b<-MRCA(phylo_kelp, raft_data_poha$sample_name[raft_data_poha$Pop %in% c("Doubtful_Sound","Snares","StewartIsland")])

raft_data_north<-raft_data[raft_data$Species=="Ant North",]
c<-MRCA(phylo_kelp, raft_data_north$sample_name)

d<-MRCA(phylo_kelp, raft_data_poha$sample_name[raft_data_poha$Pop %in% c("Doubtful_Sound","Snares","StewartIsland")])
e<-MRCA(phylo_kelp, c(raft_data_poha$sample_name))

raft_data_south<-raft_data[raft_data$Species=="Ant South",]
f<-MRCA(phylo_kelp, raft_data_south$sample_name[raft_data_south$Pop %in% c("Snares")])
g<-MRCA(phylo_kelp, raft_data_south$sample_name[raft_data_south$Pop %in% c("All Day Bay","Moeraki")])
h<-MRCA(phylo_kelp, raft_data_south$sample_name[raft_data_south$Pop %in% c("Tautuku","Cannibal Bay")])
i<-MRCA(phylo_kelp, raft_data_south$sample_name[raft_data_south$Pop %in% c("Orepuki", "Charleston", "Westport")])
j<-MRCA(phylo_kelp, raft_data_south$sample_name[raft_data_south$Pop %in% c("Orepuki", "Charleston", "Westport","Smails Beach")])

k<-MRCA(phylo_kelp, raft_data_south$sample_name[raft_data_south$Pop %in% c("All Day Bay","Moeraki")])

p <- ggtree(tree_kelp) + geom_tiplab()
viewClade(p, e)
tree_kelp
tree_kelp<-unroot(tree_kelp)
#Smails Beach
raft_data$Munida
tip_labs<-data.frame(Sample=tree_kelp$tip.label)
tip_labs$Munida<-ifelse(tip_labs$Sample %in% raft_data$sample_name[raft_data$Munida=="Raft"],"Raft","Non-Raft")
tree_kelp$edge.length<-sqrt(tree_kelp$edge.length)

pastel_palette <- c("#FFB3D9", "#7b40bc", "#CCE6C8", "#F4CE97", "#B3FFDE", 
                    "#C3AED6", "#FFB399", "#AAD2FF", "#FF8C94")

phylo_kelp <- ggtree(tree_kelp, layout="equal_angle") + 
  geom_tippoint(fill=ifelse(tip_labs$Munida=="Raft","#3281A8","#5DBB63"), shape= ifelse(tip_labs$Munida=="Raft",21,23), size=3) +
  geom_hilight(node = a, type = "encircle", expand = 0.05, alpha = 0.7, to.bottom = TRUE,
               fill = "#da8c00", col = NA)  + 
  geom_hilight(node = b, type = "encircle", expand = 0.05, alpha = 0.7, to.bottom = TRUE, 
               fill = "#7b40bc", col = NA)  + 
  geom_hilight(node = c, type = "encircle", expand = 0.05, alpha = 0.7, to.bottom = TRUE,
               fill = "#CCE6C8", col = NA)  +  
  geom_hilight(node = d, type = "encircle", expand = 0.05, alpha = 0.7, to.bottom = TRUE,
               fill = "#7b40bc", col = NA)  + 
  geom_hilight(node = e, type = "encircle",
               expand = 0.05, alpha = 0.7, to.bottom = TRUE, fill = "#B3FFDE", col = NA)  + 
  geom_hilight(node = f, type = "encircle", expand = 0.05, alpha = 0.7, to.bottom = TRUE,
               fill = "#F4CE97", col = NA) +
  geom_hilight(node = g, type = "encircle", expand = 0.05, alpha = 0.7, to.bottom = TRUE,
               fill = "#C3AED6", col = NA) +
  geom_hilight(node = h, type = "encircle", expand = 0.05, alpha = 0.7, to.bottom = TRUE,
               fill ="#AAD2FF", col = NA) +
  geom_hilight(node = i, type = "encircle", expand = 0.05, alpha = 0.7, to.bottom = TRUE,
               fill = "#FF8C94", col = NA) +
  geom_hilight(node = j, type = "encircle", expand = 0.05, alpha = 0.7, to.bottom = TRUE,
               fill = "#B3FFDE", col = NA) #+ geom_tiplab()

phylo_kelp
viewClade(p, g)
#function to convert alpha colours to 'flat' solid colours
rgba2rgb <- function(x){
  
  background_RGB <- col2rgb("white")
  color_RGB <- col2rgb(as.character(x))
  alpha <- 0.7
  # get new color  
  new_col=matrix(c(
    (1 - alpha) * background_RGB[1] + alpha * color_RGB[1],
    (1 - alpha) * background_RGB[2] + alpha * color_RGB[2],
    (1 - alpha) * background_RGB[3] + alpha * color_RGB[3]),
    nrow=3,ncol=1,dimnames=list(c("red","green","blue"))
  )
  return(rgb(t(new_col/255)))
}

library(phyloseq)
ASV_InPlace<-readr::read_rds("./ASV_InPlace.RDS")
ASV_Munida_Kelp<-readr::read_rds("ASV_Munida_Kelp.RDS")
all_samp_dist_phy<-merge_phyloseq(ASV_Munida_Kelp,ASV_InPlace)
kelp_sampdat<-rbind(data.frame(Sample=ASV_Munida_Kelp@sam_data$Specific.Sample,Raft="Raft",Location="Raft",
                               Species=ASV_Munida_Kelp@sam_data$GeneticSpecies),
                    data.frame(Sample=ASV_InPlace@sam_data$Specific.Sample,Raft="Non-Raft",Location=as.vector(ASV_InPlace@sam_data$Location), Species = ASV_InPlace@sam_data$GeneticSpecies))
kelp_sampdat<-unique(kelp_sampdat)
library(phyloseq)
all_samp_dist_phy<-phyloseq::merge_samples(all_samp_dist_phy,
                                           group=all_samp_dist_phy@sam_data$Specific.Sample,
                                           fun=mean)

all_samp_dist<-phyloseq::distance(all_samp_dist_phy,method="bray")
tree_micro<-(ape::nj(all_samp_dist))
tree_micro
#pastel_palette <- c("#FFB3D9", "#7b40bc", "#CCE6C8", "#F4CE97", "#B3FFDE", "#C3AED6", "#FFB399", "#AAD2FF", "#FF8C94")
kelp_sampdat$RGB_Hex <- pastel_palette[match(kelp_sampdat$Location, unique(kelp_sampdat$Location))]
kelp_sampdat$Shape<-ifelse(kelp_sampdat$Raft=="Raft",23,21)
#kelp_sampdat<-read.csv("kelp_sampdat.csv")
kelp_sampdat_cols<-data.frame(Sample=tree_micro$tip.label)
kelp_sampdat_cols$Col<-kelp_sampdat$RGB_Hex[match(kelp_sampdat_cols$Sample,kelp_sampdat$Sample)]
kelp_sampdat_cols$Col[kelp_sampdat_cols$Col=="#00008b"]<-"#3281A8"
kelp_sampdat_cols$Raft<-ifelse(kelp_sampdat_cols$Col=="#3281A8","Raft","Non-Raft")
phylo_micro <- ggtree(tree_micro, layout="equal_angle") + 
  #geom_tippoint(fill=ifelse(kelp_sampdat$Raft=="Raft","#3281A8","#5DBB63"), 
  #              shape= ifelse(kelp_sampdat$Raft=="Raft",21,23), size=3) +  geom_tiplab()# +
  geom_tippoint(col= kelp_sampdat_cols$Col,fill=kelp_sampdat_cols$Col,
                shape= ifelse(kelp_sampdat_cols$Raft=="Raft",21,23), size=3) #+ 
#  geom_tiplab(kelp_sampdat_cols$Sample)
phylo_micro
#params for legend
boxsize <- 10000
kelp_sites$fill<-ifelse(kelp_sites$type=="raft","#00003b","white")
leg_loc <- data.frame(x = rep(NZ_projext[1] + 50000, length(SST_levels)-1), y = seq(NZ_projext[4] - 300000, NZ_projext[4] - 80000, length.out = length(SST_levels)-1), fill = rev(grey_pal), text = c(paste(SST_levels[3:length(SST_levels)-1]-0.01, SST_levels[3:length(SST_levels)]-0.01, sep = "-"), paste0(SST_levels[length(SST_levels)]-0.01, "+")))
target_shp$col<-c("#da8c00", "#C3AED6", "#7b40bc", "#AAD2FF", "#7b40bc", "#FF8C94", "#da8c00", 
                  "#da8c00", "#CCE6C8", "#7b40bc", "#F4CE97", "#7b40bc", "#da8c00", "#B3FFDE")
traj_shps<-list.files("./data/")
traj_shp <- vect("data/rafts/MunP1/MunP2_reason1_att3_withage.shp")#[sample(1:length(traj_shp), 10)] #random subset of 100 traj
traj_shp<-traj_shp[sample(1:length(traj_shp), 10)]

traj_shp_muf7 <- vect("data/rafts/MunF7/MunF7_reason1_att3_withage.shp")
traj_shp_muf7<-traj_shp_muf7[sample(1:length(traj_shp_muf7), 10)]

traj_shp_mu6 <- vect("data/rafts/MunP1/Mu6_reason1_att3_withage.shp")
traj_shp_mu6<-traj_shp_mu6[sample(1:length(traj_shp_mu6), 10)]



traj_shp<-rbind(traj_shp_muf7,traj_shp,traj_shp_mu6)


map_plot <- ggplot(data = SST_contour) +
  geom_spatvector(aes(fill = cols), fill = as.character(SST_contour$cols), alpha = 0.5, col = NA) +
  geom_spatvector(data = grat, col = "grey60") +
  geom_spatvector(data = traj_shp, alpha = 0.8, aes(col = as.factor(reason)), col = sapply(traj_shp$reason, function(x){ifelse(x == 1, "grey10", "grey40")}), linetype = sapply(traj_shp$reason, function(x){ifelse(x == 1, "solid", "12")})) +
  geom_spatvector(data = NZ_shp, fill = "grey97", col = "grey50") +
  geom_spatvector(data = target_shp, aes(col = target_shp$col, fill = target_shp$col), 
                  alpha = 0.9, linewidth = 0.6, linetype = "12", fill = target_shp$col, col = target_shp$col) +
  geom_spatvector(data = kelp_sites, size = 3, aes(pch = type,fill=fill), 
                  stroke = 2,inherit.aes=F) +
  scale_shape_manual(values = c(21, 23)) +
  scale_fill_manual(values = c("#3281A8","#5DBB63")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), legend.position = "none") +
  coord_sf(crs = NZGD2000, expand = FALSE, xlim = c(NZ_projext[1] + 10000, NZ_projext[2] -10000), ylim = c(NZ_projext[3] + 10000, NZ_projext[4] -10000)) +
  annotate(geom = "rect", xmin = leg_loc$x[1]-boxsize*2, xmax = leg_loc$x[1]+boxsize*18, ymin = leg_loc$y[1]-boxsize*2, ymax = leg_loc$y[nrow(leg_loc)]+boxsize*5, fill = "white", col = "black", alpha = 0.7) + 
  annotate(geom = "rect", alpha=0.5,xmin = leg_loc$x-boxsize, xmax = leg_loc$x+boxsize, ymin = leg_loc$y-boxsize, ymax = leg_loc$y+boxsize, fill = leg_loc$fill, col = "black") + 
  annotate(geom = "text", x = leg_loc$x+boxsize*3, y = leg_loc$y-500, label = leg_loc$text, hjust = 0) + 
  annotate(geom = "text", x = leg_loc$x[1]-boxsize, y = leg_loc$y[nrow(leg_loc)]+boxsize*3, label = "σSST", hjust = 0)
map_plot


fig1_panels<-map_plot + phylo_kelp + phylo_micro + plot_layout(ncol = 2, design = c("
AB
AC")) + plot_annotation(tag_levels = list("A"))
fig1_panels
#ggsave(filename = "xfig1_trajs.pdf",fig1_panels)
