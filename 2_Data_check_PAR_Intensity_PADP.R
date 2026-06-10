########################################## 
###### Data cleaning and quality check ###### 
###### plot vertical profiles of PAR, LiDAR intensity, plant area density (PAD) in one graph
###### compare distribution and exclude crowns where 
###### the three vertical profiles markedly deviated from each other.
###### Determine crown base height based on PAR profile
library(lidR)
library(leafR)
library(rlas)
library(canopyLazR)
library(data.table)
library(sf)
library(fields)
library(spam)
library(dotCall64)
library(grid)
library(plyr)
library(ggplot2)
library(viewshed3d)
library(fs)
library(tidyverse)
library(stringr) 
library(gridExtra)    
library(dplyr)
library(MASS)
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

setwd("G:/My Drive")


mytheme <- list(
  theme_classic()+ 
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=13))
)

las_chm_treepoly <- list.files("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/LiDAR/UKFS/2022/Product/CHM99",
                               pattern = "las$",full.names = T,recursive=TRUE)

# This data are provided by NEON detailed in the metadata
par_all_diurnal = read.csv("F:/NEON/AOP/UKFS/Outputs/par_all_diurnal.csv")

# Create an empty LAI list and associated crown polygon global ID to store all PAI values later 
LAI <- c()
polyID <- c()
PAD25_list <- c()
PAD50_list <- c()
PAD75_list <- c()
species_list <- c()
filename_list <- c()
plotID_list <- c()
subplotID_list <-c()

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
  #lai_raster <- lai.raster(VOXELS_PAD99,min = 3)
  # Set 3m as the initial understory threshold
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
  #fwrite(LAI_all,file = "F:/NEON/AOP/UKFS/Outputs/LAI_all.csv")
  
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
  
  # Calculate standard error based on the info in the dataframe
  stats_sum$se <- stats_sum$point_count_sd / (2.5*sqrt(stats_sum$height_bin_count))
  
  # Combine the stats dataframe with PADP dataframe based on height 
  PAD_sum <- join(stats_sum,LADen,by ="height_bin")
  
  # Remove na and Inf se rows
  PAD_sum <- PAD_sum[!is.infinite(PAD_sum$se),]
  PAD_sum <- na.omit(PAD_sum)
  
  # Sort the dataframe by height
  PAD_sum <- PAD_sum[order(PAD_sum$height_bin,decreasing=F),]
  
  # Calculate accumulative PAD and graph it by %heiht - PADP
  PAD_sum[,"accummulative_PAD"] <- cumsum(PAD_sum$lad)
  PAD_sum[,"accummulative_PAD"] <- (PAD_sum$accummulative_PAD - min(PAD_sum$accummulative_PAD))/ (max(PAD_sum$accummulative_PAD) - min(PAD_sum$accummulative_PAD))
  PAD_sum[,"percent_height"]<- (PAD_sum$height_bin - min(PAD_sum$height_bin))/(max(PAD_sum$height_bin) - min(PAD_sum$height_bin))
  
  #### KOLMOGOROV-SMIRNOV TEST for all PADP distribution 
  sample1 <- PAD_sum$lad
  sample2 <- rnorm(length(PAD_sum$lad))
  group <- c(rep(paste0(filename), length(sample1)), rep("Normal distribution", length(sample2)))
  dat <- data.frame(KSD = c(sample1,sample2), group = group)
  ks <- ks.test(PAD_sum$lad,"pnorm")
  
  skewness <- skewness(PAD_sum$lad)
  kurtosis <- kurtosis(PAD_sum$lad)
  # skewness2 <- skewness * (sqrt(length(PAD_sum$lad) * (length(PAD_sum$lad) - 1)) / (length(PAD_sum$lad) - 2))
  # skewness3 <- skewness * ((1 - 1 / length(PAD_sum$lad))^1.5)
  
  # Retrieve species info for graphing
  # This is provided as metadata within the NEON foliar sampling dataset, where crown polygons are included.
  attr <- read.csv("F:/NEON/AOP/UKFS/NEON.D06.UKFS.DP1.10026.001.vst_mappingandtagging.expanded.20231031T220656Z.CSV")
  info <- attr[,c("individualID","taxonID","scientificName","plotID","subplotID")]
  info$crownID <- substr(info$individualID, start = 19,stop =max(nchar(info$individualID)))
  crownID <- str_sub(filename,start = 6, -6)
  species <- info$taxonID[which(info$crownID == crownID) ]  
  plotID <- info$plotID[which(info$crownID == crownID) ]  
  subplotID <- info$subplotID[which(info$crownID == crownID) ]  
  
  ######## Prepare: Intensity  ########
  Intensity = chm99_las$Intensity
  Z = chm99_las$Z
  NumberOfReturns = chm99_las$NumberOfReturns
  LidIntensity = cbind(Z,Intensity,NumberOfReturns)
  LidIntensity = as.data.frame(LidIntensity)
  LidIntensity = LidIntensity[order(LidIntensity$Z), ] 
  LidIntensity =na.omit(LidIntensity)
  # Exclude potential outliers
  LidIntensity = LidIntensity[LidIntensity$Z >= 0,]
  LidIntensity = LidIntensity[LidIntensity$Z <= 50,]
  LidIntensity = LidIntensity[LidIntensity$Intensity <= 256,]
  LidIntensity[,"percent_height"]<- (LidIntensity$Z - min(LidIntensity$Z))/(max(LidIntensity$Z) - min(LidIntensity$Z))
  LidIntensity$percent_height = LidIntensity$percent_height * 100
  
  ######## Prepare: COMBINED LINE GRAPH INTENSITY, PAD AND PAR ########
  ### COMBINE PROFILES OF INTENSITY, PAD AND PAR
  # We must make the variables to have the same names. So we'll just call everything X and Y
  # we also need to scale them in order to graph in a same plot
  LidIntensity_copy = LidIntensity
  PAD_sum_copy = PAD_sum
  par_all_diurnal_copy = par_all_diurnal
  
  PAD_sum_copy$y = PAD_sum_copy$lad
  PAD_sum_copy$x = PAD_sum_copy$height_bin
  
  LidIntensity_copy$y = LidIntensity_copy$Intensity / 195 # normalize by max intensity 
  LidIntensity_copy$x = LidIntensity$Z
  
  coeff = max(par_all_diurnal_copy$PARMean) # Value used to transform the data
  par_all_diurnal_copy$y = par_all_diurnal_copy$PARMean / coeff
  par_all_diurnal_copy$x = par_all_diurnal_copy$height
  
  ######## Graph: smoothed lines of PAD, intensity and PAR in one plot ######## 
  All_profiles = ggplot(PAD_sum_copy, aes(x = x, y = y),fill="blue")+
    mytheme +
    geom_smooth(method = "loess",se= F,span=0.7,alpha=0.2 ) +
    scale_x_continuous(breaks=seq(0, 60,10),expand = c(0, 0)) +
    #scale_y_continuous(breaks=seq(0, 1.3,0.2),expand = c(0, 0), limits = c(0,1.3)) +
    geom_smooth(method = "loess",data = par_all_diurnal_copy,se= F,span=0.25,alpha=0.2,color="red")+
    geom_smooth(method = "loess",data = LidIntensity_copy,se= F,span=0.3,alpha=0.2,color="orange")+
    coord_flip()+
    scale_y_continuous(
      # Features of the first axis
      name = "PAD (m2/m3) and point cloud intensity (lm/m2)",
      breaks=seq(0, 1,0.2),expand = c(0, 0), limits = c(0,1),
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*coeff, name="Mean Diurnal PAR (?mol m-2 s-1)")
    ) + 
    xlab("Canopy height (m)") +
    ggtitle(paste0("2022 Profiles of /n PAD, Diurnal PAR & LiDAR Intensity of: /n", 
                   species, " at ",plotID, " (plot) - ", subplotID, " (subplot)")) +
    theme(plot.title = element_text(size=9),text = element_text(size=8)) +
    geom_segment(aes(x = 60,y = 0.05, xend = 60,yend = 0.12),size = 1,colour = "blue") +
    annotate(geom = 'text', label = "PAD",
             hjust = 0,x = 60, y= 0.18, size = 3,colour="blue" ) +
    geom_segment(aes(x = 58,y = 0.05, xend = 58,yend = 0.12),size = 1,colour = "red") +
    annotate(geom = 'text', label = "PAR",
             hjust = 0,x = 58, y= 0.18, size = 3,colour="red" ) +
    geom_segment(aes(x = 56,y = 0.05, xend = 56,yend = 0.12),size = 1,colour = "orange") +
    annotate(geom = 'text', label = "Intensity",
             hjust = 0, x = 56, y= 0.18, size = 3,colour="orange" ) 
  
  # Save to local for data quality check
  ggsave(plot = All_profiles,filename = paste0("F:/NEON/AOP/UKFS/Outputs/Data_check/all_profiles/all_profiles_",filename,"_",species,"_",
                                               plotID,"_", subplotID, ".jpeg"))
   
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
  intensity_stats <- intensity_height_stats
  intensity_stats <- intensity_stats[!duplicated(intensity_stats), ]
  intensity_stats <- as.data.frame(intensity_stats)
  colnames(intensity_stats) <- c("height","intensity","height_bin","intensity_sd")
  intensity_stats_sum <- join(intensity_stats, countHeight_list,by="height_bin")
  # Remove potential na and outliers
  intensity_stats_sum <- na.omit(intensity_stats_sum)
  intensity_stats_sum <- subset(intensity_stats_sum,intensity_stats_sum$intensity <= 6272)
  
  # Calculate standard error based on the info in the dataframe
  intensity_stats_sum$se <- intensity_stats_sum$intensity_sd / sqrt(intensity_stats_sum$height_bin_count)
  
  # Remove na and Inf se rows and unreasonable values
  intensity_stats_sum <- intensity_stats_sum[!is.infinite(intensity_stats_sum$se),]
  intensity_stats_sum <- na.omit(intensity_stats_sum)
  
  # Sort the dataframe by height
  intensity_stats_sum <- intensity_stats_sum[order(intensity_stats_sum$height_bin,decreasing=F),]
  # Remove duplicated height_bin rows
  intensity_stats_sum <- intensity_stats_sum[!duplicated(intensity_stats_sum$height_bin), ]
  
  # Normalize data and prepare for graphing
  intensity_stats_sum$y = ((intensity_stats_sum$intensity) / 195) * -1
  intensity_stats_sum$x = intensity_stats_sum$height_bin
  intensity_stats_sum$se1 = ((intensity_stats_sum$se) / 195) * -1
  
  ######## Graph: PADP and Intensity on either side  of a common Y axis  ######## 
  ######## Place PAR profile on the right   ######## 
  PADP_Intensity <- ggplot(PAD_sum_copy, aes(x = x, y = y))+
    mytheme +
    geom_bar(stat = "identity",position=position_dodge(),alpha = 0.8) +
    geom_smooth(color= "darkgreen", se= F,alpha=0.2,span = 0.42) +
    # geom_ribbon(aes(ymin = 0,ymax = predict(loess(lad ~ height_bin,span = 0.42 ))),
    #             alpha = 0.1,fill = 'green') +
    scale_x_continuous(breaks=seq(0, max(par_all_diurnal_copy$x),5),expand = c(0, 0), limits = c(0,max(par_all_diurnal_copy$x))) +
    #scale_x_continuous(breaks=seq(5, max(PAD_sum_copy$height_bin)+5,5),expand = c(0, 0), limits = c(5,max(PAD_sum_copy$height_bin)+5)) +
    #scale_y_continuous(breaks=seq(-1, max(PAD_sum_copy$lad + PAD_sum_copy$se),0.1),expand = c(0, 0), limits = c(-1,max(PAD_sum_copy$lad + PAD_sum_copy$se)+0.1)) +
    theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
    geom_errorbar(data = intensity_stats_sum,aes(ymin=y, ymax=y+se1), alpha = 0.7, width=.2,position=position_dodge(.9)) +
    ylab("LiDAR intensity (lm/m2)   PAD (m2/m3)") +
    xlab("Canopy Height (m)") +
    #ggtitle(paste0("PADP 2022: ", species, " at ",plotID, " (plot) - ", subplotID, " (subplot)")) +
    theme(plot.title = element_text(size = 9),
          axis.text=element_text(size=8), 
          axis.title=element_text(size=8),
          axis.title.x = element_text(hjust = 0.5, size = 7)) +
    geom_bar(data = intensity_stats_sum,stat = "identity",position=position_dodge(),alpha = 0.8) +
    geom_smooth(data = intensity_stats_sum,color= "orange", se= F,alpha=0.2,span = 0.42) +
    #theme(axis.text.x = element_text(angle = 95, vjust = 1, hjust=1))+
    coord_flip()+
    geom_segment(aes(x = 0,y = 0, xend = max(par_all_diurnal_copy$x),yend = 0),size = 0.5,color="black") +
    # geom_ribbon(data = intensity_stats_sum,aes(ymin = 0,ymax = predict(loess(y ~ x,span = 0.42 ))),
    #             alpha = 0.1,fill = 'yellow')+
    # scale_y_continuous(
    #   # Features of the first axis
    #   name = "PAD (m2/m3) and point cloud intensity (lm/m2)",
    #   breaks=seq(0, 1,0.2),expand = c(0, 0), limits = c(-1,1),
    #   # Add a second axis and specify its features
    #   sec.axis = sec_axis((~.-1)*coeff, name="LiDAR intensity (lm/m2)",breaks=seq(-300, 0,50))
    # )
    scale_y_continuous(breaks=seq(-1.2, 1,0.2),expand = c(0, 0), limits = c(-1.2,1),
                       label = c("300", "250","200", "150", "100","50","0.0","0.2","0.4","0.6","0.8"," "))+
    theme(plot.margin=unit(c(0,0,0,0.2), "cm")) # c(top, right, bottom, left)
  
  PAR = ggplot(par_all_diurnal, aes(x = height, y = PARMean))+
    mytheme +
    geom_smooth(method = "loess",se= F,span=0.3,alpha=0.3,color="red") +
    scale_x_continuous(breaks=seq(0, max(par_all_diurnal$height), 5),expand=c(0,0),limits = c(0,NA)) +
    scale_y_continuous(breaks=seq(0, max(par_all_diurnal$PARMean), 100),expand=c(0,0),limits = c(0,NA)) +
    ylab("Mean PAR (?mol m-2 s-1)") +
    xlab(" ") +
    theme(plot.title = element_text(size = 9),
          axis.text=element_text(size=8), 
          axis.title=element_text(size=8),axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
    coord_flip()+
    theme( plot.margin=unit(c(0,0.2,0,0), "cm"))
  
  grid.arrange(PADP_Intensity,PAR,ncol=2, top=textGrob(paste0("PADP 2022: ", species, " at ",plotID, " (plot) - ", subplotID, " (subplot)")))
  ggsave(plot = grid.arrange(PADP_Intensity,PAR,ncol=2, top=textGrob(paste0("PADP 2022: ", species, " at ",plotID, " (plot) - ", subplotID, " (subplot)"))),filename = paste0("F:/NEON/AOP/UKFS/Outputs/Data_check/all_profiles/PADP_Intensity_PAR_",filename,"_",species,"_",
                                                                                                                                                                               plotID,"_", subplotID, ".jpeg"))
  
}



