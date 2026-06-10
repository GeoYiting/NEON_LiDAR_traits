########################################## 
###### Compute imaging spectrometer indices ###### 
library(raster)
library(rhdf5)
library(terra)
library(stars)
library(glcm)  
############ PART 1: FIRST, choose the bands we want and ############ 
############ convert them to .tif from HDF5 ############ 
wd = "D:/Huge dataset/UKFS/2022/DP3.30006.002/neon-aop-provisional-products/2022/FullSite/D06/2022_UKFS_6/L3/Spectrometer/Reflectance"
setwd(wd)
hyper_im <- list.files("D:/Huge dataset/UKFS/2022/DP3.30006.002/neon-aop-provisional-products/2022/FullSite/D06/2022_UKFS_6/L3/Spectrometer/Reflectance",pattern=".h5$",full.names = F)


for (j in 1:length(hyper_im))
{
  myEPSG <- h5read(hyper_im[j],"/UKFS/Reflectance/Metadata/Coordinate_System/EPSG Code" )
  myCRS <- crs(paste0("+init=epsg:",myEPSG))
  reflInfo <- h5readAttributes(hyper_im[j],"/UKFS/Reflectance/Reflectance_Data" )
  # Grab the UTM coordinates of the spatial extent
  xMin <- reflInfo$Spatial_Extent_meters[1]
  xMax <- reflInfo$Spatial_Extent_meters[2]
  yMin <- reflInfo$Spatial_Extent_meters[3]
  yMax <- reflInfo$Spatial_Extent_meters[4]
  # define the extent (left, right, top, bottom)
  rasExt <- extent(xMin,xMax,yMin,yMax)
  myNoDataValue <- as.integer(reflInfo$Data_Ignore_Value)
  
  band2Raster <- function(file, band, noDataValue, extent, CRS){
    out <- h5read(file,"/UKFS/Reflectance/Reflectance_Data",index=list(band,NULL,NULL))
    # Convert from array to matrix
    out <- (out[1,,])
    # transpose data to fix flipped row and column order 
    # depending upon how your data are formatted 
    out <- t(out)
    # assign data ignore values to NA
    out[out == myNoDataValue] <- NA
    
    # turn the out object into a raster
    outr <- raster(out,crs=CRS)
    
    # assign the extents to the raster
    extent(outr) <- extent
    
    return(outr)
  }
  
  # create a list of the bands we want in our stack
  # Below numbers stand for RGB,Red edge, NIR and SWIR
  rgb <- list(58,34,19,68,90,250) 
  
  # lapply tells R to apply the function to each element in the list
  rgb_rast <- lapply(rgb,FUN=band2Raster, file = hyper_im[j],
                     noDataValue=myNoDataValue, 
                     extent=rasExt,
                     CRS=myCRS)

  # check out the properties or rgb_nir_rast

  # create a raster stack from our list of rasters
  rgbStack <- stack(rgb_rast)
  # Create a list of band names
  bandNames <- paste("Band_",unlist(rgb),sep="")
  # set the rasterStack's names equal to the list of bandNames created above
  names(rgbStack) <- bandNames
  
  filename <- substring(basename(hyper_im[j]), 1, nchar(basename(hyper_im[j]))-3)
  rgbStack <- st_as_stars(rgbStack)
  write_stars(rgbStack, paste0("D:/Huge dataset/UKFS/2022/DP3.30006.002/Spectrometer/RGB_IR_RE_",filename,".tif"))

}


############ PART 3: Compute VIs  #############
treepoly <- list.files("G:/Shared drives/NSF_MSB-NES_ Tree_Crowns/AOP/Crown_polygons/Foliar_traits_DP1.10026.001/UKFS/Split_poly",
                       pattern = "shp$",full.names = T,recursive=TRUE)

# Read in the mosaiced .tif
im <- stack("D:/Huge dataset/UKFS/2022/DP3.30006.002/mosaic_masked.tif")
# Extract each band
RED = im$mosaic_masked_1/10000
GREEN =im$mosaic_masked_2/10000
BLUE = im$mosaic_masked_3/10000
REDEDGE = im$mosaic_masked_4/10000
NIR = im$mosaic_masked_5/10000
SWIR = im$mosaic_masked_6/10000
# Remove erroneous values
values(RED)[values(RED) < 0] = NA
values(RED)[values(RED) > 1] = NA
values(GREEN)[values(GREEN) < 0] = NA
values(GREEN)[values(GREEN) > 1] = NA
values(BLUE)[values(BLUE) < 0] = NA
values(BLUE)[values(BLUE) > 1] = NA
values(NIR)[values(NIR) < 0] = NA
values(NIR)[values(NIR) > 1] = NA
values(REDEDGE)[values(REDEDGE) < 0] = NA
values(REDEDGE)[values(REDEDGE) > 1] = NA
values(SWIR)[values(SWIR) < 0] = NA
values(SWIR)[values(SWIR) > 1] = NA

###### Calculate VIs
ndvi = (NIR - RED) / (NIR + RED)
values(ndvi)[values(ndvi) > 1] = NA
values(ndvi)[values(ndvi) < -1] = NA

nirv = ndvi * NIR 
values(nirv)[values(nirv) > 1] = NA
values(nirv)[values(nirv) < -1] = NA

evi =  2.5 * ((NIR - RED) / (NIR + 6*RED - 7.5*BLUE + 1)) 
values(evi)[values(evi) > 1] = NA
values(evi)[values(evi) < -1] = NA

# Red edge chlorophyll Index
CIred.edge = NIR/REDEDGE - 1
values(CIred.edge)[values(CIred.edge) < 0] = NA
values(CIred.edge)[values(CIred.edge) > 15] = NA

# green chlorophyll Index
CIgreen = NIR/GREEN -1
values(CIgreen)[values(CIgreen) < 0] = NA
values(CIgreen)[values(CIgreen) > 15] = NA

# Normalized Difference Red-Edge Index 
NDRE = (NIR - REDEDGE) / (NIR +REDEDGE)
values(NDRE)[values(NDRE) > 1] = NA
values(NDRE)[values(NDRE) < -1] = NA

###### Calculate WIs
ndwi = (GREEN - NIR) / (GREEN + NIR)
values(ndwi)[values(ndwi) > 1] = NA
values(ndwi)[values(ndwi) < -1] = NA

ndmi = (GREEN - SWIR) / (GREEN + SWIR)
values(ndmi)[values(ndmi) > 1] = NA
values(ndmi)[values(ndmi) < -1] = NA


##### Extract VI values for each tree crown
library(sf)

polyID <- c()
NDVI_mean <- c()
NIRv_mean <- c()
EVI_mean<- c()
rededgeCI_mean <- c()
greenCI_mean <- c()
NDRE_mean <- c()
NDWI_mean <- c()
NDMI_mean <- c()

for (i in 1:length(treepoly)) {
  treecrown <- st_read(treepoly[i])
  filename <- substring(basename(treepoly[i]), 1, nchar(basename(treepoly[i]))-4)
  polyID <- append(polyID, filename)
  
  # ndvi = rast(ndvi)
  # nirv = rast(nirv)
  # evi = rast(evi)
  ndvi.vals <- raster::extract(ndvi, treecrown,df = TRUE)
  nirv.vals <- raster::extract(nirv, treecrown,df = TRUE)
  evi.vals <- raster::extract(evi, treecrown,df = TRUE)
  CIred.edge.vals <- raster::extract(CIred.edge, treecrown,df = TRUE)
  CIgreen.vals <- raster::extract(CIgreen, treecrown,df = TRUE)
  NDRE.vals <- raster::extract(NDRE, treecrown,df = TRUE)
  ndwi.vals <- raster::extract(ndwi, treecrown,df = TRUE)
  ndmi.vals <- raster::extract(ndmi, treecrown,df = TRUE)
    
  ndvi.mean = mean(ndvi.vals$layer,na.rm =T )
  nirv.mean = mean(nirv.vals$layer,na.rm =T)
  evi.mean = mean(evi.vals$layer,na.rm =T)
  CIre.mean = mean(CIred.edge.vals$layer,na.rm =T)
  CIg.mean = mean(CIgreen.vals$layer,na.rm =T)
  NDRE.mean = mean(NDRE.vals$layer,na.rm =T)
  ndwi.mean = mean(ndwi.vals$layer,na.rm =T)
  ndmi.mean = mean(ndmi.vals$layer,na.rm =T)
  
  NDVI_mean <- append(NDVI_mean,ndvi.mean)
  NIRv_mean <- append(NIRv_mean,nirv.mean)
  EVI_mean<- append(EVI_mean,evi.mean)
  rededgeCI_mean <- append(rededgeCI_mean,CIre.mean)
  greenCI_mean <- append(greenCI_mean,CIg.mean)
  NDRE_mean <- append(NDRE_mean,NDRE.mean)
  NDWI_mean <- append(NDWI_mean,ndwi.mean)
  NDMI_mean <- append(NDMI_mean,ndmi.mean)
  
  hyper_VI_all <- cbind(polyID, NDVI_mean,
                        NIRv_mean,
                        EVI_mean,
                        rededgeCI_mean,
                        greenCI_mean,
                        NDRE_mean,
                        NDWI_mean,
                        NDMI_mean)
  
  hyper_VI_all= as.data.frame(hyper_VI_all)
  
  colnames(hyper_VI_all) = c("crownPolygonID", "NDVI_mean",
                             "NIRv_mean",
                             "EVI_mean",
                             "rededgeCI_mean",
                             "greenCI_mean",
                             "NDRE_mean",
                             "NDWI_mean",
                             "NDMI_mean")
  
}    







