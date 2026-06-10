########################################## 
###### 3. Calculate LAI, Aint10/50/75, APAD 10/50/75 ###### 
###### Use base height defined by PAR to cut off understoreys by
###### manually setting the threshold based on the step 1 graphs
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
library(gridExtra)    
library(dplyr)
library(MASS)
library(tidyverse)
library(forcats)
library(psych)
library(multcomp)
library(report)
library(moments)
library(ggpubr)
library(nortest)
library(lattice)
library(car)
library(FSA)
library(see)

setwd("F:/NEON/AOP/UKFS")

mytheme <- list(
  theme_classic()+ 
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=13))
)

las_chm_treepoly <- list.files("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/LiDAR/UKFS/2022/Product/CHM99",
                               pattern = "las$",full.names = T,recursive=TRUE)

par_all_diurnal = read.csv("F:/NEON/AOP/UKFS/Outputs/par_all_diurnal.csv")

# Create an empty LAI list and associated crown polygon global ID to store all PAI values later 
LAI <- c()
polyID <- c()
PAD10_list <- c()
PAD50_list <- c()
PAD75_list <- c()
species_list <- c()
filename_list <- c()
plotID_list <- c()
subplotID_list <-c()
AInt10_list <- c()
AInt50_list <- c()
AInt75_list <- c()

for (i in 1:length(las_chm_treepoly)) {
  ########### Voxelizaion and LAI
  # Read in the saved .las
  chm99_las <- readLAS(las_chm_treepoly[i])
  chm99_las <- filter_poi(chm99_las, Z >= 0, Z <= 50)
  # Creates a data frame of the 3D voxels information (xyz) 
  # with PAD values from the .las file
  VOXELS_PAD99 <- lad.voxels(las_chm_treepoly[i], grain.size = 1,k=0.5)
  # Calculate PADP 
  PADP99 <- lad.profile(VOXELS_PAD99, relative = F)
  lai_raster <- lai.raster(VOXELS_PAD99,min = 3)
  # Set 3m as the initial threshold for understory
  PADP99 <- subset(PADP99,PADP99$height >2.5)
  
  # Calculate LAI from PAD
  PAI99 <- lai(PADP99) 
  PAI99_overstory <- lai(PADP99,min = 3)
  PAI99_overstory_df <- as.data.frame(PAI99_overstory)
  # Save the results
  filename <- substring(basename(las_chm_treepoly[i]), 7, nchar(basename(las_chm_treepoly[i]))-4)
  polyID <- append(polyID, filename) # Append filename to the empty list
  LAI <- append(LAI,PAI99_overstory_df) # Append LAI to the empty list
  LAI_all <- cbind(polyID, LAI) # Combine global ID and PAI for all crowns
  LAI_all <- as.data.frame(LAI_all)
  fwrite(LAI_all,file = "F:/NEON/AOP/UKFS/2022/Outputs/LAI_all.csv")
  
  ########### Prepare: PAD and PADP   ########### 
  # Calculate PADP and store as a dataframe
  Zs <- chm99_las@data$Z
  Zs <- Zs[!is.na(Zs)]
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
  colnames(LADen) <- c("height_bin","lad")
  
  # Create a list representing all height bins
  countHeight_list <- setNames(vector(mode = "list", length(LADen$height_bin)), 1:length(LADen$height_bin))
  for (i in 1:length(countHeight_list)) {
    h <- countHeight_list[i]
    h <- c()
  }
  
  # count the number of height bins of the same height from the voxelization
  for (i in 1:length(countHeight_list)) {
    countHeight_list[i] <- sum(Zs >i - 1 & Zs<i+1)
  }
  countHeight_list <- t(as.data.frame(countHeight_list))
  countHeight_list <- as.data.frame(countHeight_list)
  colnames(countHeight_list) <- "height_bin_count"
  countHeight_list$height_bin <- c(1:nrow(countHeight_list))
  countHeight_list$height_bin <- countHeight_list$height_bin +2.5
  
  # count the number of lidar returns of the same height from the voxelization
  den_height <-cbind(chm99_las@data$Z, chm99_las@data$NumberOfReturns)
  den_height[den_height == 0 | den_height <0] <- NA
  den_height <- na.omit(den_height)
  den_height <- as.data.frame(den_height)
  colnames(den_height) <- c("height","point_count")
  den_height$height_bin <- trunc(den_height$height)+0.5
  
  # Calculate sd of returns count per voxel by group (each height bin)
  den_height_stats <- den_height %>% 
    group_by(height_bin) %>%
    mutate(sd = sd(point_count)) %>% 
    as.data.frame()
  
  # Combine all teh above info together
  stats <- cbind(den_height_stats$height_bin,den_height_stats$sd)
  stats <- stats[!duplicated(stats), ]
  stats <- as.data.frame(stats)
  colnames(stats) <- c("height_bin","point_count_sd")
  stats_sum <- join(stats, countHeight_list,by="height_bin")
  stats_sum <- na.omit(stats_sum)
  
  # Calculate standard error based on the indo in the dataframe
  stats_sum$se <- stats_sum$point_count_sd / (2.5*sqrt(stats_sum$height_bin_count))
  
  # Combine the stats dataframe with PADP dataframe based on height 
  PAD_sum <- join(stats_sum,LADen,by ="height_bin")
  
  # Remove na and Inf se rows
  PAD_sum <- PAD_sum[!is.infinite(PAD_sum$se),]
  PAD_sum <- na.omit(PAD_sum)
  
  # Sort the dataframe by height
  PAD_sum <- PAD_sum[order(PAD_sum$height_bin,decreasing=T),]
  
  # Calculate APAD from crown top and graph it by top % of the height
  PAD_sum[,"accummulative_PAD"] <- cumsum(PAD_sum$lad)
  PAD_sum[,"accummulative_PAD"] <- (PAD_sum$accummulative_PAD - min(PAD_sum$accummulative_PAD))/ (max(PAD_sum$accummulative_PAD) - min(PAD_sum$accummulative_PAD))
  PAD_sum[,"height_percentile_from_top"]<- (PAD_sum$height_bin - max(PAD_sum$height_bin))/(max(PAD_sum$height_bin) - min(PAD_sum$height_bin))
  PAD_sum$height_percentile_from_top = PAD_sum$height_percentile_from_top * -1
  
  
  # Retrieve species info so we can put on the graphs
  attr <- read.csv("F:/NEON/AOP/UKFS/NEON.D06.UKFS.DP1.10026.001.vst_mappingandtagging.expanded.20231031T220656Z.csv")
  info <- attr[,c("individualID","taxonID","scientificName","plotID","subplotID")]
  info$crownID <- substr(info$individualID, start = 19,stop =max(nchar(info$individualID)))
  crownID <- str_sub(filename,start = 6, -6)
  species <- info$taxonID[which(info$crownID == crownID) ]  
  plotID <- info$plotID[which(info$crownID == crownID) ]  
  subplotID <- info$subplotID[which(info$crownID == crownID) ]  
  
 
  ######## Prepare: Intensity profile   ########
  Intensity = chm99_las$Intensity
  Z = chm99_las$Z
  NumberOfReturns = chm99_las$NumberOfReturns
  LidIntensity = cbind(Z,Intensity,NumberOfReturns)
  LidIntensity = as.data.frame(LidIntensity)
  LidIntensity = LidIntensity[order(-LidIntensity$Z), ] 
  LidIntensity =na.omit(LidIntensity)
  LidIntensity = LidIntensity[LidIntensity$Z >= 0,]
  LidIntensity = LidIntensity[LidIntensity$Z <= 50,]
  LidIntensity = LidIntensity[LidIntensity$Intensity <= 256,]
  LidIntensity[,"height_percentile_from_top"]<- -1 * (  (LidIntensity$Z - max(LidIntensity$Z))/(max(LidIntensity$Z) - min(LidIntensity$Z)) )
  LidIntensity$height_percentile_from_top = LidIntensity$height_percentile_from_top * 100
  LidIntensity[,"accumulative_intensity"] = cumsum(LidIntensity$Intensity)
  LidIntensity[,"accumulative_intensity"] <- (LidIntensity$accumulative_intensity - min(LidIntensity$accumulative_intensity))/ (max(LidIntensity$accumulative_intensity) - min(LidIntensity$accumulative_intensity))
  LidIntensity$accumulative_intensity = LidIntensity$accumulative_intensity * 100
  
  ### Extract Accumulative Intensity at top 10% 50% 75% height 
  # Find the value closest to 50 since there's no 50 in the dataframe
  closest<-function(closest_val,my_value){
    closest_val[which(abs(closest_val-my_value)==min(abs(closest_val-my_value)))] }
  
  h_closest50_intensity <- closest(LidIntensity$height_percentile_from_top,50)
  Intensity50 <- LidIntensity$accumulative_intensity[which(LidIntensity$height_percentile_from_top == h_closest50_intensity[1])]
  
  h_closest10_intensity <- closest(LidIntensity$height_percentile_from_top,10)
  Intensity10 <- LidIntensity$accumulative_intensity[which(LidIntensity$height_percentile_from_top == h_closest10_intensity[1])]
  
  h_closest75_intensity <- closest(LidIntensity$height_percentile_from_top,75)
  Intensity75 <- LidIntensity$accumulative_intensity[which(LidIntensity$height_percentile_from_top == h_closest75_intensity[1])]
  
  ######## Prepare: Overlay of PAD, intensity and PAR   ########
  ### OVERLAY PROFILE LINES OF INTENSITY, PAD AND PAR
  # We must make the variables to have the same names. So we'll just call everything X and Y
  # we also need to scale them in order to graph in a same plot
  LidIntensity_copy = LidIntensity
  PAD_sum_copy = PAD_sum
  par_all_diurnal_copy = par_all_diurnal
  
  #PAD_sum_copy[,"norm_PAD"] <- (PAD_sum_copy$lad - min(PAD_sum_copy$lad))/ (max(PAD_sum_copy$lad) - min(PAD_sum_copy$lad))
  PAD_sum_copy$y = PAD_sum_copy$lad
  PAD_sum_copy$x = PAD_sum_copy$height_bin
  
  #norm_Intensity =(LidIntensity_copy$Intensity - min(LidIntensity_copy$Intensity))/ (max(LidIntensity_copy$Intensity) - min(LidIntensity_copy$Intensity))
  LidIntensity_copy$y = LidIntensity_copy$Intensity / 195
  LidIntensity_copy$x = LidIntensity$Z
  
  #norm_PAR<- (par_all_diurnal_copy$PARMean- min(par_all_diurnal_copy$PARMean))/ (max(par_all_diurnal_copy$PARMean) - min(par_all_diurnal_copy$PARMean))
  coeff = max(par_all_diurnal_copy$PARMean) # Value used to transform the data
  par_all_diurnal_copy$y = par_all_diurnal_copy$PARMean / coeff
  par_all_diurnal_copy$x = par_all_diurnal_copy$height
  
  ######## Prepare: PADP and Intensity on either side of a common Y axis ######## 
  # lidar intensity of the same height from the voxelization
  intensity_height <-cbind(chm99_las@data$Z, chm99_las@data$Intensity)
  intensity_height[intensity_height <0] <- NA
  intensity_height <- na.omit(intensity_height)
  intensity_height <- as.data.frame(intensity_height)
  colnames(intensity_height) <- c("height","intensity")
  intensity_height$height_bin <- trunc(intensity_height$height)+0.5
  
  # Calculate sd of intensity by group (each height bin)
  intensity_height_stats <- intensity_height %>% 
    group_by(height_bin) %>%
    mutate(sd = sd(intensity)) %>% 
    as.data.frame()
  
  # Combine all the above info together
  # intensity_stats <- cbind(intensity_height_stats$height_bin,intensity_height_stats$sd)
  intensity_stats <- intensity_height_stats
  intensity_stats <- intensity_stats[!duplicated(intensity_stats), ]
  intensity_stats <- as.data.frame(intensity_stats)
  colnames(intensity_stats) <- c("height","intensity","height_bin","intensity_sd")
  intensity_stats_sum <- join(intensity_stats, countHeight_list,by="height_bin")
  intensity_stats_sum <- na.omit(intensity_stats_sum)
  intensity_stats_sum <- subset(intensity_stats_sum,intensity_stats_sum$intensity <= 195)
  
  # Calculate standard error based on the info in the dataframe
  intensity_stats_sum$se <- intensity_stats_sum$intensity_sd / sqrt(intensity_stats_sum$height_bin_count)
  
  # Remove na and Inf se rows and unreasonable values
  intensity_stats_sum <- intensity_stats_sum[!is.infinite(intensity_stats_sum$se),]
  intensity_stats_sum <- na.omit(intensity_stats_sum)
  
  # Sort the dataframe by height
  intensity_stats_sum <- intensity_stats_sum[order(intensity_stats_sum$height_bin,decreasing=T),]
  # Remove duplicated height_bin rows
  intensity_stats_sum <- intensity_stats_sum[!duplicated(intensity_stats_sum$height_bin), ]
  
  # Prepare graph
  intensity_stats_sum$y = ((intensity_stats_sum$intensity) / 195) * -1
  intensity_stats_sum$x = intensity_stats_sum$height_bin
  intensity_stats_sum$se1 = ((intensity_stats_sum$se) / 195) * -1
  
  ########  Use base height defined by PAR to cut off understoreys   ########  
  #### Do this by manually check the PAR profile !!!
  #### The base height for UKFS is 5m 
  PAD_sum_copy1 = subset(PAD_sum_copy, PAD_sum_copy$height_bin>4.5)
  intensity_stats_sum1 = subset(intensity_stats_sum, intensity_stats_sum$height_bin>5)
  
  par_all_diurnal_copy1 = par_all_diurnal_copy
  
  # Jump out of the loop if the crown maximum height is less than the threshold
  if (length(PAD_sum_copy1$height_bin) < 4 ) {next}
  
  ########### Prepare: PAD and PADP   ########### 
  # Calculate PADP and store as a dataframe
  Zs <- chm99_las@data$Z
  Zs <- Zs[!is.na(Zs)]
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
  colnames(LADen) <- c("height_bin","lad")
  
  # Create a list representing all height bins
  countHeight_list <- setNames(vector(mode = "list", length(LADen$height_bin)), 1:length(LADen$height_bin))
  for (i in 1:length(countHeight_list)) {
    h <- countHeight_list[i]
    h <- c()
  }
  
  # count the number of height bins of the same height from the voxelization
  for (i in 1:length(countHeight_list)) {
    countHeight_list[i] <- sum(Zs >i - 1 & Zs<i+1)
  }
  countHeight_list <- t(as.data.frame(countHeight_list))
  countHeight_list <- as.data.frame(countHeight_list)
  colnames(countHeight_list) <- "height_bin_count"
  countHeight_list$height_bin <- c(1:nrow(countHeight_list))
  countHeight_list$height_bin <- countHeight_list$height_bin +2.5
  
  # count the number of lidar returns of the same height from the voxelization
  den_height <-cbind(chm99_las@data$Z, chm99_las@data$NumberOfReturns)
  den_height[den_height == 0 | den_height <0] <- NA
  den_height <- na.omit(den_height)
  den_height <- as.data.frame(den_height)
  colnames(den_height) <- c("height","point_count")
  den_height$height_bin <- trunc(den_height$height)+0.5
  
  # Calculate sd of returns count per voxel by group (each height bin)
  den_height_stats <- den_height %>% 
    group_by(height_bin) %>%
    mutate(sd = sd(point_count)) %>% 
    as.data.frame()
  
  # Combine all teh above info together
  stats <- cbind(den_height_stats$height_bin,den_height_stats$sd)
  stats <- stats[!duplicated(stats), ]
  stats <- as.data.frame(stats)
  colnames(stats) <- c("height_bin","point_count_sd")
  stats_sum <- join(stats, countHeight_list,by="height_bin")
  stats_sum <- na.omit(stats_sum)
  
  # Calculate standard error based on the indo in the dataframe
  stats_sum$se <- stats_sum$point_count_sd / (2.5*sqrt(stats_sum$height_bin_count))
  
  # Combine the stats dataframe with PADP dataframe based on height 
  PAD_sum <- merge(stats_sum,LADen,by ="height_bin" )
  
  # Remove na and Inf se rows
  PAD_sum <- PAD_sum[!is.infinite(PAD_sum$se),]
  PAD_sum <- na.omit(PAD_sum)
  
  # Sort the dataframe by height
  PAD_sum <- PAD_sum[order(PAD_sum$height_bin,decreasing=T),]
  PAD_sum = subset(PAD_sum,PAD_sum$height_bin > 4.5)
  
  
  # Calculate accumulative PAD from crown top and graph it by top % of the height
  PAD_sum[,"accummulative_PAD"] <- cumsum(PAD_sum$lad)
  PAD_sum[,"accummulative_PAD"] <- (PAD_sum$accummulative_PAD - min(PAD_sum$accummulative_PAD))/ (max(PAD_sum$accummulative_PAD) - min(PAD_sum$accummulative_PAD))
  PAD_sum[,"height_percentile_from_top"]<- (PAD_sum$height_bin - max(PAD_sum$height_bin))/(max(PAD_sum$height_bin) - min(PAD_sum$height_bin))
  PAD_sum$height_percentile_from_top = PAD_sum$height_percentile_from_top * -1
  
  PADP_accumulative <- ggplot(PAD_sum, aes(x = height_percentile_from_top*100, y = accummulative_PAD*100))+
    mytheme +
    geom_smooth(method = "loess",se= T) +
    scale_x_continuous(breaks=seq(0, 100, 10),expand = c(0, 0)) +
    scale_y_continuous(breaks=seq(0, 100, 10),expand = c(0, 0)) +
    ylab("Accumulative PAD (%)") +
    xlab("Top Height Percentile of the Canopy (%)") +
    ggtitle(paste0("PADP 2022: ", species, " at ",plotID, " (plot) - ", subplotID, " (subplot)")) +
    theme(plot.title = element_text(size = 9),
          axis.text=element_text(size=8), 
          axis.title=element_text(size=8)) +
    annotate(geom = 'text', label = paste0(species," ", plotID," - ",subplotID ), 
             x = 95, y= 8.5, size = 4) +
    coord_flip()
  
  ### Extract PAD at 50% height from the lowess fit curve
  bb <- ggplot_build(PADP_accumulative)
  ## extract the right component, just the x/y coordinates
  out <- bb$data[[1]]
  out <- out[,c("x","y")]
  # Find the value closest to 50 since there's no 50 in the dataframe
  closest<-function(closest_val,my_value){
    closest_val[which(abs(closest_val-my_value)==min(abs(closest_val-my_value)))] }
  
  h_closest <- closest(out$x,50)
  # PAD50 <- out$y[which(out$x == h_closest[1])]
  PAD50 <- median(out$y)
  
  h_closest10 <- closest(out$x,10)
  PAD10 <- out$y[which(out$x == h_closest10[1])]
  
  h_closest75 <- closest(out$x,75)
  PAD75 <- out$y[which(out$x == h_closest75[1])]
  
  # Save APAD10/50/75 and AInt10/50/75
  species_list <- append(species_list,species)
  filename_list <- append(filename_list,filename)
  plotID_list <- append(plotID_list,plotID)
  subplotID_list <- append(subplotID_list,subplotID)
  PAD10_list <- append(PAD10_list,PAD10)
  PAD50_list <- append(PAD50_list,PAD50)
  PAD75_list <- append(PAD75_list,PAD75)
  AInt10_list <- append(AInt10_list,Intensity10[1])
  AInt50_list <- append(AInt50_list,Intensity50[1])
  AInt75_list <- append(AInt75_list,Intensity75[1])
  
  PAD_accumulative_10_50_75 <- cbind(filename_list,species_list,PAD10_list,PAD50_list,PAD75_list,
                                     AInt10_list,AInt50_list,AInt75_list,
                                     plotID_list,subplotID_list)
  
  colnames(PAD_accumulative_10_50_75) <-
    c("crownPolygonID","taxonID","APAD10","APAD50","APAD75",
      "AInt10","AInt50","AInt75", "plotID","subplotID")
  
  write.csv(PAD_accumulative_10_50_75,file = "F:/NEON/AOP/UKFS/2022/Outputs/PAD_accumulative_10_50_75.csv")
  
}




