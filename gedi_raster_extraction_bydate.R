###############################################################################
###                                                                         ###
### GEDI Processing Script 2:                                               ###
### Polygon Value Extraction Script                                         ###
### Re-Configured by Neal Swayze 5/3/2022                                   ###
### Modified by Steven Filippelli 8/31/2022                                 ###
###                                                                         ###
###############################################################################

###_________________________###
### INSTALL / LOAD PACKAGES ###
###_________________________###

library(pacman)
p_load(parallel, foreach, data.table, tidyverse, exactextractr, sf, terra, rgeos, rgdal, mapview, doParallel, sfarrow, lubridate)

###______________________________
###
### CONFIGURE Input/Output ###
###_______________________________###

### ESTABLISH YOUR ROOT DIRECTORY
rootDir <- ("D:/ECOFOR")
setwd(rootDir)

point_path <- r"(J:\projects\ECOFOR\gedi\gedi_data\04_gedi_filtered_data_shp\GEDI_2AB_2019to2023_leafon_sampy500m.parquet)"
pred_basedir <- file.path(rootDir, "lt")


outbasename <- tools::file_path_sans_ext(basename(point_path))
spectral_data_dir = file.path(rootDir, "gedi", "extracted")
dir.create(spectral_data_dir)


###___________________________________________________________________________________________________###
### READ IN THE GEDI SAMPLES FOR EACH YEAR, BUFFER BY FOOTPRINT SIZE (12.5m) ###
###___________________________________________________________________________________________________###

### GEDI Data
df <- st_read_parquet(point_path)
st_crs(df) <- 32636 #4326       # Setting CRS since it was lost somehow
# df <- sample_n(df, 1000) # Test subset

df$date <- as.POSIXlt(df$delta_time, format = "%y-%m-%d")
df$doy <- yday(df$date)
df$year <- year(df$date)

# Rain year is defined as the year beginning with the start of the dry season (121-273) and the following wet season (274-120)
# e.g., rain year 2018 is May 1, 2018 - April 30 2019
df$rainyear <- df$year
df[df$doy<121, "rainyear"] <- as.data.frame(df)[df$doy<121,"rainyear"] - 1

# df <- st_transform(df, crs=CRS("+init=epsg:32636")) # Reproject to match spectral data if loaded points are in WGS84
df_buff <-  st_buffer(df[,c("shot_number", "year", "rainyear")], 12.5,
                      endCapStyle="ROUND", joinStyle = "ROUND", nQuadSegs = 20)

years <- sort(unique(df$rainyear))
seasons = c("wet", "dry")
combos = expand.grid(years, seasons)

### Auto detect the cores and register parallel session
cores = detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

df_ext <- foreach(i=1:nrow(combos), .combine =rbind, .inorder = FALSE, .errorhandling="remove") %dopar% {
  ### Get the necessary packages
  library(terra)
  library(rgdal)
  library(sf)
  library(rgeos)
  library(exactextractr)
  
  terraOptions(memfrac = 0.15, todisk = FALSE)
  
  # Set up inputs
  year <- combos[i,1]
  season <- combos[i,2]
  rast_path <-  file.path(pred_basedir, season, paste0("lt_",season,"_",as.character(year),".vrt"))
  raster_layer <- rast(rast_path)
  ydf = df_buff[df_buff$rainyear==year,]
  ydf$season = season
  
  ### Extract average pixel values for the intersection of GEDI polygons in the full raster layer
  keep_cols = c("shot_number", "year", "season")
  extracted_values <- exactextractr::exact_extract(raster_layer, ydf, fun='mean', append_cols = keep_cols) #, max_cells_in_memory = 2.4e+09) # actually faster to read per feature than load whole raster at once
  rm(raster_layer, ydf)
  gc()
  return(extracted_values)
}

### Stop the processing cluster
stopCluster(cl)

### Check the extracted spectral dataset
str(df_ext)

# split into wet and dry season cols
df_wet <- df_ext[df_ext$season=="wet",]
colnames(df_wet) = gsub("mean.", "wet_", colnames(df_wet))

df_dry <- df_ext[df_ext$season=="dry",]
colnames(df_dry) = gsub("mean.", "dry_", colnames(df_dry))

### Merge the two seasons for export
dfout <- left_join(select(df_dry,-c("season", "year")),
                   select(df_wet,-c("season", "year")), by="shot_number")

# dfout_geo <- left_join(df[c('shot_number', 'geometry')],
#                     dfout, by='shot_number')

# Write out data as csv, gpkg, and parquet
setwd(spectral_data_dir)
fwrite(as.data.frame(dfout), file = paste0(outbasename,"_lt.csv"), sep = ",")
# st_write_parquet(dfout_geo, paste0(outbasename,"_lt.parquet"))
# st_write(dfout_geo, paste0(outbasename,"_lt.gpkg"), append=FALSE)
