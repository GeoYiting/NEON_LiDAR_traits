########################################## 
###### 4. All other crown traits ###### 
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


mytheme <- list(
  theme_classic()+ 
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=13))
)

########### PAI
las_chm_treepoly <- list.files("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/LiDAR/UKFS/2022/Product/CHM99",
                               pattern = "las$",full.names = T,recursive=TRUE)

PAI_AOP <- c()
polyID <- c()

for (i in 1:length(las_chm_treepoly)) {
  a <- readLAS(las_chm_treepoly[i]) 
  a <- filter_poi(a, Z >= 5, Z <= 50) # Adjust base height accordingly
  Zs <- a@data$Z
  Zs <- Zs[!is.na(Zs)]
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=5) 
  PAI <- sum(LADen$lad, na.rm=TRUE) 
  PAI_AOP <- append(PAI_AOP,PAI) # Append PAI to the empty list
  
  # Save the results
  filename <- substring(basename(las_chm_treepoly[i]), 7, nchar(basename(las_chm_treepoly[i]))-4)
  polyID <- append(polyID, filename) # Append filename to the empty list
  # Combine global ID and PAI for all crowns
  PAI_all <- cbind(polyID, PAI_AOP) 
  write.csv(PAI_all,file = "F:/NEON/AOP/UKFS/2022/Outputs/PAI_all.csv")
}

########### LAI, top rugosity, crown rugosity, MAXCH, MOCH, vertical variance, rumple
# Read in the 99% CHM .las as a list
las_chm_treepoly <- list.files("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/LiDAR/UKFS/2022/Product/CHM99",
                               pattern = "las$",full.names = T,recursive=TRUE)
polyID <- c()
LAI <- c()
mean_top_rug_list <- c()
mean_rug_list <- c()
MAXCH_list <- c()
MOCH_list <- c()
vert_var_list <- c()
rumple_list <- c()

for (i in 1:length(las_chm_treepoly)) {
  ########### LAI
  # Read in the saved .las
  chm99_las <- readLAS(las_chm_treepoly[i])
  chm99_las <- filter_poi(chm99_las, Z >= 5, Z <= 50) # Adjust base height accordingly
  # Creates a data frame of the 3D voxels information (xyz) 
  # with PAD values from the .las file
  VOXELS_PAD99 <- lad.voxels(las_chm_treepoly[i], grain.size = 1,k=0.5)
  # Calculate PADP
  PADP99 <- lad.profile(VOXELS_PAD99, relative = F)
  # Retrieve only overstory PADP 
  # Adjust base height accordingly
  PADP99 <- subset(PADP99,PADP99$height >4.5)
  # Total LAI 
  LAI99 <- lai(PADP99)
  # Real trees
  LAI99_overstory <- lai(PADP99,min = 5)# Adjust base height accordingly
  LAI99_overstory_df <- as.data.frame(LAI99_overstory)
  LAI <- append(LAI,LAI99_overstory_df) # Append LAI to the empty list
  LAI_all <- as.data.frame(LAI_all)
  
  filename <- substring(basename(las_chm_treepoly[i]), 7, nchar(basename(las_chm_treepoly[i]))-4)
  polyID <- append(polyID, filename) # Append filename to the empty list
  
  # Combine global ID and PAI for all crowns
  LAI_all <- cbind(polyID, LAI) 
  
  ########### Top rugosity
  # Create CHM raster from the .las
  chm <- grid_canopy(readLAS(las_chm_treepoly[i]),res = 1,p2r())
  # Calculate top rugosity from the raster
  top_rugosity <- toc.rugosity(chm.raster = chm,xy.res = 1,z.res = 1) 
  mean_top_rug <- cellStats(top_rugosity,stat = "mean") 
  mean_top_rug_list <- append(mean_top_rug_list,mean_top_rug)
  
  ########### Total rugosity
  # Reads in each .las file from the.las list and
  # converts each one into a voxelized array based on the x,y,z resolutions you specify
  # We specify 1 x 1 x 1 m throughout our project, which is more than 20 points per square meters
  laz_array <- laz.to.array(las_chm_treepoly[i],voxel.resolution = 1, z.resolution = 1, use.classified.returns = F)
  leveld_lidar_array <- canopy.height.levelr(laz_array)
  PAD_array <- machorn.lad(leveld_lidar_array, voxel.height =  1, beer.lambert.constant = 0.5)
  # Calculate crown rugosity from the PAD array
  rugosity <- rugosity.within.canopy(lad.array = PAD_array,laz.array = laz_array,ht.cut = 3,epsg.code = 32615)
  mean_rug <- cellStats(rugosity$rugosity.raster,stat = "mean")
  mean_rug_list <- append(mean_rug_list,mean_rug)
  
  
  ########### Vertical variance
  # Retrieve PAD profile from voxelization 
  PAD99 <- PAD_array$rLAD 
  # Vertical variance of PAD
  vertical_sd_PAD99 <- apply(PAD99, c(1,2), var,na.rm=TRUE)  
  # Vertical variance - mean of the vertical PAD variation
  vert_var <- mean(vertical_sd_PAD99,na.rm=TRUE)  
  # Horizontal variation of the vertical variation, which is crown rugusity
  # rugosity99 <- var(vertical_sd_PAD99,na.rm=TRUE) 
  vert_var_list <- append(vert_var_list,vert_var)
  
  ########### Structure complexity
  # MEAN OUTER CANOPY HEIGHT (MOCH) 
  MOCH <- mean(chm@data@values, na.rm = TRUE) 
  MOCH_list <- append(MOCH_list,MOCH)
  #MAXCH
  MAXCH <- max(chm@data@values, na.rm=TRUE) 
  MAXCH_list <- append(MAXCH_list,MAXCH)
  #RUMPLE
  #calculate rumple, a ratio of outer canopy surface area to 
  #ground surface area (1600 m^2)
  rumple <- rumple_index(chm) 
  rumple_list <- append(rumple_list,rumple)
  
  # OUTPUT all CALCULATED METRICS INTO A TABLE
  CT <- 
    cbind(LAI_all, mean_top_rug_list, mean_rug_list, MAXCH_list, MOCH_list, vert_var_list, rumple_list) 
  
  colnames(CT) <- 
    c("crownPolygonID", "LAI", "Top rugosity","Crown rugosity", "MAXCH","MOCH", "Vertical variance",
      "Rumple") 
  
}

# Merge with PAI we generated before
# Join by the common column (global ID)
CT_df <- as.data.frame(CT)
PAI_all_df <- as.data.frame(PAI_all) 
# Rename the polygon ID column to have a same column name with the CT table for joining 
colnames(PAI_all_df) <- c("crownPolygonID", "PAI")
# Specify all.x = F if you have removed some crowns 
# since they are not removed from the PAI file
CT_final <- merge(PAI_all_df,CT_df,by = "crownPolygonID", all.x= F) 


########### Append species info and APAD ########### 
CT_df <-  read.csv("F:/NEON/AOP/UKFS/2022/Outputs/Crown_traits.csv")
# Pull up the attributes file, which is from NEON foliar sampling dataset
# Look for file name with "mappingandtagging"
attr <- read.csv("F:/NEON/AOP/UKFS/NEON.D06.UKFS.DP1.10026.001.vst_mappingandtagging.expanded.20231031T220656Z.CSV")
# find out the attributes we want to append to CT file, as well as the common column (global ID) the join is based on
names(attr) 
# create a new dataframe
species <- attr[,c("individualID","domainID","siteID","taxonID","scientificName","plotID","subplotID","nestedSubplotID")]
# Join the two tables using the column "crownPolygonID" from the crown traits .csv 
# and "individualID" from the attributes .csv 
# Extract the 5 or 6 digits from the last part of the "individualID" column 
# and from the "crownPolygonID" column
# That is the unique ID we use for table join
species$crownID <-substr(species$individualID, start = 19,stop =max(nchar(species$individualID)))
CT_df$crownID <- substr(CT_df$crownPolygonID,start = 6, stop =max(nchar(species$individualID)))
CT_df$crownID <- str_sub(CT_df$crownID, 1, -6) # remove the year part
Crown_traits <- join(CT_df,species, by= "crownID")
Crown_traits <- as.data.frame(Crown_traits)

# Add in PAD10/50/75
PAD_accumulative_10_50_75 <- read.csv("F:/NEON/AOP/UKFS/2022/Outputs/PAD_accumulative_10_50_75.csv")
PAD_accumulative_10_50_75 = as.data.frame(PAD_accumulative_10_50_75)
Crown_traits <- join(Crown_traits,PAD_accumulative_10_50_75,by = "crownPolygonID")
# Drop Duplicated Columns:
Crown_traits <- Crown_traits[!duplicated(as.list(Crown_traits))]
# Remove NA rows
Crown_traits=Crown_traits[!is.na(Crown_traits$APAD50),] 
