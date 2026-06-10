########################################## 
###### LiDAR Pre-processing ###### 
###### This part extract .las for each delineated crown 
###### and denoise the .las

library(data.table)
library(sf)
library(canopyLazR)
library(fields)
library(spam)
library(dotCall64)
library(grid)
library(rlas)
library(plyr)
library(ggplot2)
library(lidR)
library(leafR)
library(viewshed3d)
library(fs)
library(tidyverse)
library(stringr) 

setwd("F:/NEON/AOP/UKFS")

############ 1. Split MultiPolygon shapefile into separate Polygons ############ 
# Read shp
# This is the crown delineation shapefile file provided by NEON detailed in metadata  
singlepoly_UKFS = st_read("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/Crown_polygons/Foliar_traits_DP1.10026.001/UKFS/filesToStack10026/NEON.D06.UKFS.DP1.10026.001.2023-06.expanded.20231031T220656Z.PROVISIONAL/UKFS-2023-polygons/UKFS-2023-polygons.shp")
# Set coordinate system using EPGS code
singlepoly_UKFS <- spTransform(singlepoly_UKFS,"+init=epsg:32617")
#get the names of the attribute table
names(singlepoly_UKFS)

#select the column of the attribute table that will determine the split of the shp
#which is "crownPolygonID"
unique <- unique(singlepoly_UKFS$crownPolyg)
#create new polygons based on the determined column
for (i in 1:length(unique)) {
  tmp <- singlepoly_UKFS[singlepoly_UKFS$crownPolyg == unique[i], ]
  # readOGR(tmp, dsn="G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/Crown_polygons/Foliar_traits_DP1.10026.001/UKFS/Split_poly", unique[i], driver="ESRI Shapefile",
  #          overwrite_layer=TRUE)
  st_write(tmp, dsn="G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/Crown_polygons/Foliar_traits_DP1.10026.001/UKFS/Split_poly", unique[i], driver="ESRI Shapefile",
           overwrite_layer=TRUE)
}

########### 2. Extract lidar point cloud for each crown polygon ############ 
las_UKFS <- readLAS("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/LiDAR/UKFS/las_clipped/las_UKFS.las")

treepoly <- list.files("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/Crown_polygons/Foliar_traits_DP1.10026.001/UKFS/Split_poly",
                       pattern = "shp$",full.names = T,recursive=TRUE)
# Loop through all shapefiles to extract point cloud one by one
for (i in 1:length(treepoly)) {
  # Extract polygon filename for naming the individual .las file to be extracted 
  filename <- basename(treepoly[i])
  filename <- substring(filename, 1, nchar(filename)-4)
  treepolygon  <- st_read(treepoly[i]) # Read in polygon shapefile one at a time
  split_las <- clip_roi(las = las_UKFS, treepolygon) # Extract .las within the polygon
  # Filter out polygon where there are no .las
  if (!is.null(split_las) == T) {
    writeLAS(split_las,file = paste0("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/LiDAR/UKFS/2022/las_clipped/Split_las/",filename,".las"))
  }
}

########### 3. Lidar pre-processing for crown traits ########### 
########### 3.1 Generate CHM las ########### 
las_treepoly <- list.files("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/LiDAR/UKFS/2022/las_clipped/Split_las",
                           pattern = "las$",full.names = T,recursive=TRUE)
for (i in 1:length(las_treepoly)) {
  ## (1) Denoise 
  # Classify noise signals using sor algorithm
  las_treepolygon <- readLAS(las_treepoly[i])
  las <- classify_noise(las_treepolygon,  sor(15,7))
  # Remove outliers 
  las_denoise <- filter_poi(las, Classification != LASNOISE)
  ## (2) Classify ground returns  
  parameters_pmf <- util_makeZhangParam(b = 1,dh0 = 0.1, dhmax = 3,s=, max_ws = 12, exp = F)
  las_class <- classify_ground(las_denoise, pmf(parameters_pmf$ws, parameters_pmf$th)) 
  ## (3) Normalize by DTM to generate CHM
  dtm <- grid_terrain(las_class, algorithm = knnidw(k = 6L, p = 2))
  chm_las <-normalize_height(las_class, dtm, method = "bilinear")
  chm <- grid_canopy(chm_las,res = 1,p2r()) 
  chm_las <- filter_poi(chm_las, Z >= 0, Z <= 50)
  # Grab the name and write to the output
  filename <- substring(basename(las_treepoly[i]), 1, nchar(basename(las_treepoly[i]))-4)
  writeLAS(chm_las, file = paste0("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/LiDAR/UKFS/2022/Product/CHM_",filename,".las"))
  # Use lidar points only up to 99% height to remove uncertainty. See Appendix table A1 for details.
  chm99_las <- filter_poi(chm_las, Z<= max(chm_las$Z) * 0.99) 
  chm99 <- grid_canopy(chm99_las,res = 1,p2r()) # Generate a CHM raster
  # Save the 99% cleaned-up CHM 
  writeLAS(chm99_las, file = paste0("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/LiDAR/UKFS/2022/Product/CHM99/CHM99_", filename,".las"))
}


