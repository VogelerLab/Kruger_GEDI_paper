###############################################################################
###                                                                         ###
### GEDI Processing Script 2:                                               ###
### Polygon Value Extraction Script                                         ###
### Re-Configured by Neal Swayze 5/3/2022                                   ###
### POC: nealswayze1@gmail.com    
###
### Modified by Steven Filippelli 8/24/2023
###                                                                         ###
###############################################################################

###_________________________###
### INSTALL / LOAD PACKAGES ###
###_________________________###

library(pacman)
p_load(data.table, tidyverse, exactextractr, sf, terra, rgeos, rgdal, mapview, sfarrow, lubridate)

terraOptions(memfrac = 0.15, todisk = FALSE)

###_______________________________###
### CONFIGURE Input/Output ###
###_______________________________###
point_path <- r"(C:\scratch\ECOFOR\gedi\GEDI_2AB_2019to2023_leafon_sampy500m.parquet)"
# rast_path <- r"(J:\projects\ECOFOR\topo\topo_all.vrt)"
# rast_path <- r"(J:\projects\ECOFOR\soils\soil_all.vrt)"
rast_path <- r"(J:\projects\ECOFOR\climate\worldclim_bio_all.vrt)"

### Setup output
rootDir <- ("H:/ECOFOR")
setwd(rootDir)
outbasename <- tools::file_path_sans_ext(basename(point_path))
spectral_data_dir = file.path(rootDir, "gedi", "extracted")
dir.create(spectral_data_dir)


###___________________________________________________________________________________________________###
### READ IN THE GEDI SAMPLES FOR EACH YEAR, BUFFER BY FOOTPRINT SIZE (12.5m) ###
###___________________________________________________________________________________________________###

### GEDI Data
df <- st_read_parquet(point_path)
# df <- sample_n(df, 1000) # Test subset

# st_crs(df) <- 4326       # Setting CRS since it was lost somehow
# df <- st_transform(df, crs=CRS("EPSG:32636")) # Reproject to match spectral data

df_buff <-  st_buffer(df[,c("shot_number")], 12.5,
                      endCapStyle="ROUND", joinStyle = "ROUND", nQuadSegs = 20)

rast_layer <- rast(rast_path)

### Extract average pixel values for the intersection of GEDI polygons in the full raster layer
keep_cols = c("shot_number")
df_ext <- exactextractr::exact_extract(rast_layer, df_buff, fun='mean', append_cols = keep_cols) #, max_cells_in_memory = 2.4e+09)

rm(rast_layer, df_buff)
gc()

### Check the extracted spectral dataset
colnames(df_ext) = sub("mean.", "", colnames(df_ext))
str(df_ext)


# Write out data as csv, gpkg, and parquet
setwd(spectral_data_dir)
fwrite(as.data.frame(df_ext), file = paste0(outbasename,"_climate.csv"), sep = ",")
# st_write_parquet(df_ext, paste0(outbasename,"_topo.parquet"))
# st_write(dfout, paste0(outbasename,"_topo.gpkg"), append=FALSE)
