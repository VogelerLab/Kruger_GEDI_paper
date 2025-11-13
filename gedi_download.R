###############################################################################
###                                                                         ###
### Turnkey GEDI Downloading                                                ###
### rGEDI_engine_V1.1                                                       ###
### Re-Configured by Neal Swayze 5/3/20212                                  ###
### POC: nealswayze1@gmail.com                                              ###
###                                                                         ###
### Modified 8/29/22 by Steven Filippelli to not use rGEDI and to handle    ###
### granules with no shots intersecting the study area after base filtering ###
###############################################################################

###_________________________###
### INSTALL / LOAD PACKAGES ###
###_________________________###

### Run pacman to load all required packages
library(pacman)
pacman::p_load(rgdal, mapview, sp, sf, raster, data.table, dplyr, foreach, doParallel, 
               rlas, ggplot2, httr, getPass, sys, curl, jsonlite, hdf5r)#, rGEDI)


###______________________###
### INSTRUCTIONS FOR USE ###
###______________________###

### Step 1: Download and install aria2c to your machine from https://github.com/aria2/aria2/releases/tag/release-1.35.0 -> aria2-1.35.0-win-64bit-build1.zip

### Step 2: Change line 42 to reflect the location of the aria2c.exe executable ie ("C:\\Users\\neal\\aria2-1.35.0-win-64bit-build1\\aria2c.exe") 

### Step 3: Create a folder on your desired storage drive as a base directory (example : "GEDI_processing"), this is your rootDir. Update line 59 to reflect this

### Step 4: Run LINES 49-51 of this script to generate a study_area_input subfolder in your root directory

### Step 5: Create a shapefile of your desired study area in WGS 1984 (ESPG: 4326) and move the shapefile to your root directory -> 01_study_area_input folder

### Step 6: Update LINE 54 to reflect the name of the shapefile in the 01_study_area_input folder

### Step 7: You will need to generate a .netrc file if you dont have one, the script will prompt you. Create an account at: urs.earthdata.nasa.gov

###_______________________________###
### STUDY AREA AND DATE SELECTION ###
###_______________________________###

### POINT THIS LINE AT YOUR ARIA 2 BATCH DOWNLOAD EXECUTABLE IN THE ARIA2 FOLDER
aria_exe <- "C:\\aria2-1.36.0-win-64bit-build1\\aria2c.exe"

### ESTABLISH YOUR ROOT DIRECTORY (use / instead of \\ otherwise the script fails at readLevel2A_custom)
rootDir <- ("D:/ecofor/gedi")     
setwd(rootDir)

###_______________________________###
### AUTOGENERATE BASE DIRECTORIES ###
###_______________________________###

dir.create(file.path(rootDir, "gedi_data"))
base_dir = file.path(rootDir, "gedi_data")

### GENERATE A DIRECTORY FOR YOUR STUDY AREA SHAPEFILE
dir.create(file.path(base_dir, "/01_study_area_input"))
study_area_input <- file.path(base_dir, "/01_study_area_input/")
setwd(study_area_input)

### INPUT THE NAME OF YOUR DESIRED STUDY AREA SHAPEFILE
study_area <- readOGR("greaterkruger_buf1000simp_wgs84.gpkg")
study_area <- spTransform(study_area, CRS("+init=epsg:4326"))
viewExtent(study_area, map.type = "Esri.WorldImagery")

### SET DESIRED DATE TIMEFRAME FOR GEDI FOOTPRINT EXTRACTION
daterange=c("2019-01-01","2023-12-31")

### DEFINE THE BEAMS THAT YOU WANT TO KEEP FROM THE DATA (FULL POWER BEAMS INCLUDED)
target_beams <- c("BEAM0101", "BEAM0110", "BEAM1000", "BEAM1011", "BEAM0000", "BEAM0001", "BEAM0010", "BEAM0011")

### Define the number of cores to use when processing GEDI granules (Use a fast SSD or you will be sad)
num_cores <- 18

### Starting Time 
start_time <- Sys.time()

###__________________________________________________###
### Creating a set of processing folders for outputs ###
###    MODIFICATION OF LINES BELOW NOT REQURIRED     ###
###__________________________________________________###

### CREATE DATA DOWNLOAD DIRECTORIES
dir.create(file.path(base_dir, "/02_download_files"))
dir.create(file.path(base_dir, "/02_download_files/01_GEDI_downloads_2A"))
gedi_downloads_2A <- file.path(base_dir, "/02_download_files/01_GEDI_downloads_2A")
dir.create(file.path(base_dir, "/02_download_files/01_GEDI_downloads_2B"))
gedi_downloads_2B <- file.path(base_dir, "/02_download_files/01_GEDI_downloads_2B")

### CREATE PROCESSING DIRECTORIES
gedi_processing_stats_output <- file.path(base_dir)
dir.create(file.path(base_dir, "/03_gedi_filtered_data_csv"))
gedi_processed_output <- file.path(base_dir, "/03_gedi_filtered_data_csv")
dir.create(file.path(base_dir, "/04_gedi_filtered_data_shp"))
gedi_processed_output_shp <- file.path(base_dir, "/04_gedi_filtered_data_shp")

###_________________________###
### PROCESSING STARTS BELOW ###
###_________________________###

### EXTRACTING COORDINATES FROM STUDY AREA SHAPEFILE
study_coords <- extent(study_area)
study_coords <- coordinates(study_coords)
study_coords <- as.data.frame(study_coords)
colnames(study_coords) <- c("long", "lat")
study_coords_spdf <- SpatialPointsDataFrame(cbind(study_coords$long,study_coords$lat), data=study_coords)
proj4string(study_coords_spdf) = CRS("+init=epsg:4326")
mapview(study_coords_spdf, map.type = "Esri.WorldImagery", xcol = "long", ycol = "lat")

### AUTOMATICALLY CONFIGURING BOUNDING BOX FOR GEDIFINDER FUNCTION 
lower_left_long <- study_coords$long[1]
lower_left_lat <- study_coords$lat[1]
upper_right_long <- study_coords$long[3]
upper_right_lat <- study_coords$lat[3] 
bbox <- paste0(lower_left_long, ",", lower_left_lat, ",", upper_right_long, ",", upper_right_lat)

###_____________________________###
### DEFINE GEDI FINDER FUNCTION ### 
###_____________________________###

gedifinder_custom <- function(product, bbox) {
  
  # Define the base CMR granule search url, including LPDAAC provider name and max page size (2000 is the max allowed)
  cmr <- "https://cmr.earthdata.nasa.gov/search/granules.json?pretty=true&provider=LPDAAC_ECS&page_size=2000&concept_id="
  
  # Set up dictionary where key is GEDI shortname + version and value is CMR Concept ID
  concept_ids <- list('GEDI01_B.002'='C1908344278-LPDAAC_ECS', 
                      'GEDI02_A.002'='C1908348134-LPDAAC_ECS', 
                      'GEDI02_B.002'='C1908350066-LPDAAC_ECS')
  
  # CMR uses pagination for queries with more features returned than the page size
  page <- 1
  bbox <- sub(' ', '', bbox)  # Remove any white spaces
  granules <- list()          # Set up a list to store and append granule links to
  
  # Send GET request to CMR granule search endpoint w/ product concept ID, bbox & page number
  cmr_response <- GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page))
  cmr_response
  
  # Verify the request submission was successful
  if (cmr_response$status_code==200){
    
    # Send GET request to CMR granule search endpoint w/ product concept ID, bbox & page number, format return as a list
    cmr_response <- content(GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page)))$feed$entry
    
    # If 2000 features are returned, move to the next page and submit another request, and append to the response
    while(length(cmr_response) %% 2000 == 0){
      page <- page + 1
      cmr_response <- c(cmr_response, content(GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page)))$feed$entry)
    }
    
    # CMR returns more info than just the Data Pool links, below use for loop to go through each feature, grab DP link, and add to list
    for (i in 1:length(cmr_response)) {
      granules[[i]] <- cmr_response[[i]]$links[[1]]$href
    }
    
    # Return the list of links
    return(granules)
  } else {
    
    # If the request did not complete successfully, print out the response from CMR
    print(content(GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page)))$errors)
  }
}

###___________________________________________________________________________###
### DEFINE FUNCTION TO FILTER THE DETECTED GEDI GRANULES TO DESIRED TIMEFRAME ###
###___________________________________________________________________________###

filter_date <- function(granules) {
  final_granules <- list() # Create empty list for final granules
  start_date <- daterange[1] #Get the start date from the daterange object
  end_date <- daterange[2] #Get the end date from the daterange object
  ### Loop to check if the dates for each granule are within the input date range
  for(i in seq(granules)){
    list_id <- granules[[i]] # Get the first granule in the list
    list_id_info <- unlist(strsplit(list_id, split = "/")) #Split it up 
    list_id_info <- list_id_info[8] #Get the granule date
    list_id_info <- gsub("\\.", "-", list_id_info) #change the period to a -
    result <- (list_id_info <= end_date) && (list_id_info >= start_date)
    #Conditional statement to append the granule file to the final granule list if it is within the daterange
    if (isTRUE(result)) {final_granules <- append(final_granules, list_id)}
    if (isFALSE(result)) {print("Date out of bounds")}
  }
  
  final_granules
  return(final_granules)
}

###________________________________________________________________________###
### USE THE GEDIFINDER FUNCTION TO SEARCH FOR GRANULES OVER THE STUDY AREA ###
###________________________________________________________________________###

### QUERY FOR 2A products
granules <- gedifinder_custom(product="GEDI02_A.002", bbox)
granules <- filter_date(granules)
print(sprintf("%s Version 2 2A granules found.", length(granules)))
outName = "desired_granules_2A.txt"
setwd(base_dir)
write.table(granules, outName, row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\n')

### QUERY FOR 2B products
granules <- gedifinder_custom(product="GEDI02_B.002", bbox)
granules <- filter_date(granules)
print(sprintf("%s Version 2 2B granules found.", length(granules)))
outName = "desired_granules_2B.txt"
setwd(base_dir)
write.table(granules, outName, row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\n')

###______________________________###
### CONFIGURE PATH TO NETRC FILE ###
###______________________________###

# Retrieve user directory (for netrc file)
usr <- file.path(Sys.getenv("USERPROFILE")) 

# If no user profile exists, use home directory
if (usr == "") {usr = Sys.getenv("HOME")}    

# Path to netrc file
netrc <- file.path(usr,'.netrc', fsep = .Platform$file.sep) 

###_____________________________________________###
### AUTOMATICALLY QUERY IF YOU HAVE .NETRC FILE ###
###_____________________________________________###

### If you do not have a  .netrc file with your Earthdata Login credentials stored in your home dir,
### below you will be prompted for your NASA Earthdata Login Username and Password and a netrc file
### will be created to store your credentials (home dir). Create an account at: urs.earthdata.nasa.gov

if (file.exists(netrc) == FALSE || grepl("urs.earthdata.nasa.gov", readLines(netrc)) == FALSE) {
  netrc_conn <- file(netrc)
  
  # User will be prompted for NASA Earthdata Login Username and Password below
  writeLines(c("machine urs.earthdata.nasa.gov",
               sprintf("login %s", getPass(msg = "Enter NASA Earthdata Login Username \n (or create an account at urs.earthdata.nasa.gov):")),
               sprintf("password %s", getPass(msg = "Enter NASA Earthdata Login Password:"))), netrc_conn)
  close(netrc_conn)
}
netrc <- read.table(netrc)

### Grab Username
username <- as.data.frame(netrc$V2)
username <- username[2,]
username <- noquote(as.character(username))
username

### Grab Password
password <- as.data.frame(netrc$V2)
password <- password[3,]
password <- noquote(as.character(password))
password

###____________________________________________________________________________###
### Batch Downloading Method - THIS WILL DOWNLOAD DESIRED GRANULES IN PARALLEL ###
###____________________________________________________________________________###

### Generate and execute command for batch downloading of 2A GEDI granules
aria_command_2A = paste0(aria_exe, " --max-concurrent-downloads=10 --optimize-concurrent-downloads --auto-file-renaming=false --file-allocation=none --http-auth-challenge=true --http-user=", username, " --http-passwd=", password, " --dir ", file.path(base_dir, "02_download_files", "01_GEDI_downloads_2A"), " -i ", file.path(base_dir, "desired_granules_2A.txt")," --force-save")
print(aria_command_2A)    # Print and paste to an independent command window instead
# shell(aria_command_2A)

### Generate and execute command for batch downloading of 2B GEDI granules
aria_command_2B = paste0(aria_exe, " --max-concurrent-downloads=10 --optimize-concurrent-downloads --auto-file-renaming=false --file-allocation=none --http-auth-challenge=true --http-user=", username, " --http-passwd=", password," --dir ", file.path(base_dir, "02_download_files", "01_GEDI_downloads_2B"), " -i ", file.path(base_dir, "desired_granules_2B.txt")," --force-save") #--max-concurrent-downloads=5
print(aria_command_2B)   # Print and paste to an independent command window instead
# shell(aria_command_2B)


###_______________________________________________________________________________________###
### DEFINE CUSTOM FUNCTIONS FOR READING IN GEDI GRANULES, AND COMPULE DATA FROM .H5 FILES ###
###_______________________________________________________________________________________###


### Define read level 2A Function from .h5 files
readLevel2A_custom = function (level2Apath) {
  
  #setClass("gedi.level1b", representation(h5="H5File",level1b.spdf='SpatialPointsDataFrame'))
  #' @importFrom hdf5r H5File
  setRefClass("H5File")
  requireNamespace("data.table")
  
  #' Class for GEDI level2A
  #'
  #' @slot h5 Object of class H5File from `hdf5r` package containing the
  #'GEDI level2A products: ground elevation, canopy top height, and relative heights (RH).
  #'
  #' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
  #' \url{https://lpdaac.usgs.gov/products/gedi02_av002/}
  #'
  #' @import methods
  #' @export
  gedi.level2a <- setClass(
    Class="gedi.level2a",
    slots = list(h5 = "H5File")
  )
  
  h5closeall = function(con, ...) {
    try(con@h5$close_all(), silent=TRUE)
  }
  
  #'Close hdf5 connections from gedi* objects
  #'
  #' @description 
  #' Closing files will avoid locking HDF5 GEDI files.
  #' 
  #'@param con An object of class `gedi.level*`
  #'@param ... Inherited from base
  #'
  #' @export
  #' @rdname close
  #' @method close gedi.level1b
  setGeneric("close", function(con, ...)
    standardGeneric("close"))
  
  #' Handles the [`rGEDI::gedi.level2a-class`].
  #'@rdname close
  setMethod("close", signature = c("gedi.level2a"), h5closeall)
  
  level2a_h5 <- hdf5r::H5File$new(level2Apath, mode = "r")
  level2a <- new("gedi.level2a", h5 = level2a_h5)
  return(level2a)}


### Define read level 2B Function from .h5 files
readLevel2B_custom = function (level2Bpath) {
  
  #setClass("gedi.level1b", representation(h5="H5File",level1b.spdf='SpatialPointsDataFrame'))
  #' @importFrom hdf5r H5File
  setRefClass("H5File")
  requireNamespace("data.table")
  
  #' Class for GEDI level2B
  #'
  #' @slot h5 Object of class [`H5File`][hdf5r::H5File-class] from `hdf5r` package containing the
  #'GEDI level2B products: canopy cover, Plant Area Index (PAI), Plant Area Volume Density (PAVD),
  #'and Foliage Height Diversity (FHD).
  #'
  #' @seealso [`H5File`][hdf5r::H5File-class] in the `hdf5r` package and
  #' \url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
  #'
  #' @import methods
  #' @export
  gedi.level2b <- setClass(
    Class="gedi.level2b",
    slots = list(h5 = "H5File")
  )
  
  h5closeall = function(con, ...) {
    try(con@h5$close_all(), silent=TRUE)
  }
  
  #'Close hdf5 connections from gedi* objects
  #'
  #' @description 
  #' Closing files will avoid locking HDF5 GEDI files.
  #' 
  #'@param con An object of class `gedi.level*`
  #'@param ... Inherited from base
  #'
  #' @export
  #' @rdname close
  #' @method close gedi.level1b
  setGeneric("close", function(con, ...)
    standardGeneric("close"))

  #' Handles the [`rGEDI::gedi.level2b-class`].
  #'@rdname close
  setMethod("close", signature = c("gedi.level2b"), h5closeall)
  
  
  level2b_h5 <- hdf5r::H5File$new(level2Bpath, mode = "r")
  level2b <- new("gedi.level2b", h5 = level2b_h5)
  return(level2b)}

### Define custom processing function to get desired 2A data
getLevel2A_custom = function (level2a) {
  level2a <- level2a@h5
  groups_id <- grep("BEAM\\d{4}$", gsub("/", "", 
                                        hdf5r::list.groups(level2a, recursive = F)), value = T)
  rh.dt <- data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), 
                              style = 3)
  i.s = 0
  for (i in groups_id) {
    i.s <- i.s + 1
    utils::setTxtProgressBar(pb, i.s)
    level2a_i <- level2a[[i]]
    if (any(hdf5r::list.datasets(level2a_i) == "shot_number")) {
      if (length(level2a_i[["rh"]]$dims) == 2) {
        rh = t(level2a_i[["rh"]][, ])
      }
      else {
        rh = t(level2a_i[["rh"]][])
      }
      rhs <- data.table::data.table(beam <- rep(i, length(level2a_i[["shot_number"]][])), 
                                    shot_number = level2a_i[["shot_number"]][], 
                                    degrade_flag = level2a_i[["degrade_flag"]][], 
                                    quality_flag = level2a_i[["quality_flag"]][], 
                                    deltatime = level2a_i[["delta_time"]][], 
                                    sensitivity = level2a_i[["sensitivity"]][], 
                                    solar_elevation = level2a_i[["solar_elevation"]][], 
                                    lat_lowestmode = level2a_i[["lat_lowestmode"]][], 
                                    lon_lowestmode = level2a_i[["lon_lowestmode"]][], 
                                    elev_highestreturn = level2a_i[["elev_highestreturn"]][], 
                                    elev_lowestmode = level2a_i[["elev_lowestmode"]][], 
                                    rh)
      rh.dt <- rbind(rh.dt, rhs)
    }
  }
  colnames(rh.dt) <- c("beam", "shot_number", "degrade_flag", "quality_flag", "delta_time", "sensitivity", "solar_elevation", 
                       
                       "lat_lowestmode", "lon_lowestmode", "elev_highestreturn", "elev_lowestmode", 
                       
                       paste0("rh", seq(0, 100)))
  close(pb)
  return(rh.dt)
}

### Define custom processing function to get desired 2B data
getLevel2B_custom = function(level2b){
  
  ### Function to get the full set of 2B data, not including Z profiles
  #####################################################################
  level2b<-level2b@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","", hdf5r::list.groups(level2b, recursive = F)), value = T)
  m.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0
  var.map = data.table::data.table(t(data.frame(list(
    # COL_NAMES              # H5_ADDRESS
    c("shot_number",         "shot_number"),
    c("algorithmrun_flag",   "algorithmrun_flag"),
    c("degrade_flag",        "geolocation/degrade_flag"),
    c("l2b_quality_flag",    "l2b_quality_flag"),
    c("stale_return_flag",   "stale_return_flag"),
    c("surface_flag",        "surface_flag"),
    c("solar_elevation",     "geolocation/solar_elevation"),
    c("delta_time",          "geolocation/delta_time"),
    c("sensitivity",         "sensitivity"),
    c("lat_lowestmode",      "geolocation/lat_lowestmode"),
    c("lon_lowestmode",      "geolocation/lon_lowestmode"),
    c("elev_highestreturn",  "geolocation/elev_highestreturn"),
    c("elev_lowestmode",     "geolocation/elev_lowestmode"),
    c("local_beam_elevation","geolocation/local_beam_elevation"),
    c("fhd_normal",          "fhd_normal"),
    c("pgap_theta",          "pgap_theta"),
    c("rh100",               "rh100"),
    c("pai",                 "pai"),
    c("rhov",                "rhov"),
    c("rhog",                "rhog"),
    c("omega",               "omega"),
    c("cover",               "cover")
  ))))
  
  colnames(var.map) = c("COL_NAMES", "H5_ADDRESS")
  
  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2b_i<-level2b[[i]]
    m <-data.table::data.table(beam=rep(i,length(level2b_i[["shot_number"]][])))
    
    for (row_index in 1:nrow(var.map)) {
      colname = var.map$COL_NAMES[row_index]
      h5.address = var.map$H5_ADDRESS[row_index]
      m[[colname]] <- level2b_i[[h5.address]][]
    }
    
    m.dt<-rbind(m.dt,m)
    colnames(m.dt)<-c("beam", var.map$COL_NAMES)
    main_data_table <- m.dt
    
  }
  
  close(pb)
  main_data_table
  
  ### New function to get the cover profiles (modified the PAI to get this function running)
  ##########################################################################################
  groups_id<-grep("BEAM\\d{4}$", gsub("/","", hdf5r::list.groups(level2b, recursive = F)), value = T)
  groups_id
  m.dt<-data.table::data.table()
  m.dt
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0
  
  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2b_i<-level2b[[i]]
    m<-data.table::data.table(
      beam<-rep(i,length(level2b_i[["shot_number"]][])),
      shot_number=level2b_i[["shot_number"]][],
      height_lastbin=level2b_i[["geolocation/height_lastbin"]][],
      height_bin0=level2b_i[["geolocation/height_bin0"]][],
      cover_z=t(level2b_i[["cover_z"]][,1:level2b_i[["cover_z"]]$dims[2]]))
    m.dt<-rbind(m.dt,m)
    level2b_i[["cover_z"]][,1:level2b_i[["cover_z"]]$dims[2]]
    
  }
  
  colnames(m.dt)<-c("beam","shot_number","height_lastbin","height_bin0",paste0("cover_z_",seq(0,30*5,5)[-31],"_",seq(5,30*5,5),"m"))
  cover_data_table <- m.dt
  close(pb)
  
  ### Function to get the PAI profiles
  ####################################
  groups_id<-grep("BEAM\\d{4}$",gsub("/","", hdf5r::list.groups(level2b, recursive = F)), value = T)
  m.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0
  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2b_i<-level2b[[i]]
    m<-data.table::data.table(
      beam<-rep(i,length(level2b_i[["shot_number"]][])),
      shot_number=level2b_i[["shot_number"]][],
      height_lastbin=level2b_i[["geolocation/height_lastbin"]][],
      height_bin0=level2b_i[["geolocation/height_bin0"]][],
      pai_z=t(level2b_i[["pai_z"]][,1:level2b_i[["pai_z"]]$dims[2]]))
    m.dt<-rbind(m.dt,m)
  }
  colnames(m.dt)<-c("beam","shot_number","height_lastbin",
                    "height_bin0",paste0("pai_z",seq(0,30*5,5)[-31],"_",seq(5,30*5,5),"m"))
  
  pai_data_table <- m.dt
  close(pb)
  
  ### Function to get the PAVD profiles
  #####################################
  groups_id<-grep("BEAM\\d{4}$",gsub("/","", hdf5r::list.groups(level2b, recursive = F)), value = T)
  m.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0
  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2b_i<-level2b[[i]]
    m<-data.table::data.table(
      beam<-rep(i,length(level2b_i[["shot_number"]][])),
      shot_number=level2b_i[["shot_number"]][],
      pavd_z=t(level2b_i[["pavd_z"]][,1:level2b_i[["pavd_z"]]$dims[2]]))
    m.dt<-rbind(m.dt,m)
  }
  colnames(m.dt)<-c("beam","shot_number",paste0("pavd_z",seq(0,30*5,5)[-31],"_",seq(5,30*5,5),"m"))
  pavd_data_table <- m.dt
  close(pb)
  
  ### Select the desired data from each data table
  cover_data_table <- dplyr::select(cover_data_table, subset = -c(beam, height_lastbin, height_bin0))
  pai_data_table <- dplyr::select(pai_data_table, subset = -c(beam, height_lastbin, height_bin0))
  pavd_data_table <- dplyr::select(pavd_data_table, subset = -c(beam))
  
  ### Master merge the data into a single dataframe by shot number
  master_dt <- merge(main_data_table, cover_data_table, by = "shot_number")
  master_dt <- merge(master_dt, pai_data_table, by = "shot_number")
  master_dt <- merge(master_dt, pavd_data_table, by = "shot_number")
  master_dt
  
  ### Clean up afteryourself and save memory
  rm(m.dt, main_data_table, cover_data_table, pai_data_table, pavd_data_table)
  gc()
  
  return(master_dt)
}

###______________________________________________________________###
### PROCESSING AND FILTERING LEVEL 2A DATA - OUTPUT TO DATAFRAME ###
###______________________________________________________________###

### Create temporary directory for storing 2A data  
dir.create(file.path(base_dir, "/03_gedi_filtered_data_csv/01_gedi_2A_temp"))
gedi_processed_output_2A <- file.path(base_dir, "/03_gedi_filtered_data_csv/01_gedi_2A_temp")

### Generating a list of granule files to process from the GEDI downloads folder
setwd(gedi_downloads_2A)
list <- list.files(setwd(gedi_downloads_2A), pattern = "h5$")
list

### Initiate the processing cluster
cl <- makeCluster(num_cores)
registerDoParallel(cl)

### Process the granules in parallel, return .csv for each granule to be merged later
foreach(i=1:length(list)) %dopar% {
  
  # library(rGEDI)
  library(hdf5r)  
  library(dplyr)
  library(data.table)
  library(sp)
  library(raster)
  
  setwd(gedi_downloads_2A) #set the working directory to the folder of downloaded 2A data
  granule_name <- list[[i]]
  out_name <- tools::file_path_sans_ext(granule_name)
  out_name <- paste0(out_name, ".csv")
  file_1 = readLevel2A_custom(granule_name)
  file_1 = getLevel2A_custom(file_1) #Convert the .h5 file to a dataframe in R memory for Plant area index profiles
  file_1 = subset(file_1, select = c(shot_number, beam, degrade_flag, quality_flag, sensitivity, solar_elevation, lat_lowestmode, lon_lowestmode, rh0, rh5, rh10, rh15, rh20, rh25, rh30, rh35, rh40, rh45, rh50, rh55, rh60, rh65, rh70, rh75, rh80, rh85, rh90, rh95, rh96, rh97, rh98, rh99, rh100) ) # subset data down to what we need
  file_1 = na.omit(object = file_1) #Omit missing values from the dataframe, some gedi data is missing coordinate data
  file_1$shot_number = as.character(file_1$shot_number) #Convert shot number from integer 64 to character  
  file_1 = filter(file_1, ((degrade_flag == 0) & (quality_flag == 1)) ) #Base filtering, other filtering later # (solar_elevation < 0) &  & (beam %in% target_beams) & (sensitivity >= 0.8) 
  num_rows <- nrow(file_1)
  num_rows
  
  if(num_rows > 0){
    file_1 = SpatialPointsDataFrame(cbind(file_1$lon_lowestmode,file_1$lat_lowestmode), data=file_1) #Converting dataframe to spatial points dataframe
    proj4string(file_1) = CRS("+init=epsg:4326")
    instudyarea <- !is.na(over(file_1, study_area)$id)
    if (sum(instudyarea) > 0){
      file_1 = file_1[instudyarea,]  #crop(file_1, study_area) #Neal originally cropped to the study area whereas I'm using the intersection. Crop throws an error if no points intersect the study area.
      file_1 = as.data.frame(file_1) #Converting cropped spatial points dataframe back to dataframe for further processing
      setwd(gedi_processed_output_2A) #Set the working directory to the temporary processing folder
      fwrite(file_1, file = out_name, sep = ",")
      rm(file_1) #remove file_1 for next iteration
      gc() #clear R memory for deleted files 
      
    } else {
      setwd(gedi_processed_output_2A)
      empty_df <- data.frame()
      write.csv(empty_df, file = out_name)
      rm(file_1)
      gc
    }
  }
  
  if(num_rows < 1){
    setwd(gedi_processed_output_2A)
    empty_df <- data.frame()
    write.csv(empty_df, file = out_name)
    rm(file_1)
    gc
  }
  
  out_name
}

stopCluster(cl)

### Setting up an empty storage dataframe .csv file to append level 2A data to
storage_dataframe_2A <- data.frame()                          

### Generating a list of processed 2A .csv granule files
setwd(gedi_processed_output_2A)
list <- list.files(setwd(gedi_processed_output_2A), pattern = "csv$")
list

### For loop to bind list of 2A granule csv files together
for(i in seq(list)) {
  csv_id <- list[[i]]
  granule_data <- fread(csv_id)
  num_row <- nrow(granule_data)
  if (num_row > 1){
    storage_dataframe_2A <- rbind(storage_dataframe_2A, granule_data)
    rm(granule_data)
    gc()
  }
  
}

### Write out the 2A merged GEDI data to a .csv file, remove temp directory
storage_dataframe_2A
setwd(gedi_processed_output)
fwrite(storage_dataframe_2A, "GEDI_2A_raw.csv", row.names=FALSE)
unlink(gedi_processed_output_2A, recursive = TRUE) #Remove the files in the temporary directory

###______________________________________________________________###
### PROCESSING AND FILTERING LEVEL 2B DATA - OUTPUT TO DATAFRAME ###
###______________________________________________________________###

### Create temporary directory for storing 2B data 
dir.create(file.path(base_dir, "/03_gedi_filtered_data_csv/01_gedi_2B_temp"))
gedi_processed_output_2B <- file.path(base_dir, "/03_gedi_filtered_data_csv/01_gedi_2B_temp")

### Generating a list of granule files to process from the GEDI donwloads folder
setwd(gedi_downloads_2B)
list <- list.files(setwd(gedi_downloads_2B), pattern = "h5$")
list

### Initiate the processing cluster
cl <- makeCluster(num_cores)
registerDoParallel(cl)

### Process the granules in parallel, return .csv for each granule to be merged later
foreach(i=1:length(list)) %dopar% {
  
  # library(rGEDI)
  library(dplyr)
  library(data.table)
  library(sp)
  library(raster)
  
  setwd(gedi_downloads_2B) #set the working directory to the folder of downloaded 2B data
  granule_name <- list[[i]]
  out_name <- tools::file_path_sans_ext(granule_name)
  out_name <- paste0(out_name, ".csv")
  input_file = readLevel2B_custom(list[[i]]) # Read in the first file from the list
  # TODO: getLevel2B_custom reads in a lot of columns. Why do this if much (or all) of the rows could be filtered out anyway? Try just reading coords and filter columns first, then data for the remaining rows.
  file_1 <- getLevel2B_custom(input_file) #Run custom function to get main 2B, cover Z, pai z, and pavd Z data
  file_1 = na.omit(object = file_1) #Omit missing values from the dataframe, some gedi data is missing coordinate data
  file_1$delta_time = as.Date(  as.POSIXct(file_1$delta_time,  tz = "", origin = '2018-01-01')) #converting delta time (s) to date (year-month-day)
  file_1 = filter(file_1, ( (l2b_quality_flag == 1) & (degrade_flag == 0) )) # Base filtering, solar_elevation < 0) &  & (beam %in% target_beams) & (sensitivity >= 0.8) 
  num_rows <- nrow(file_1)
  
  if(num_rows > 0){
    file_1 = SpatialPointsDataFrame(cbind(file_1$lon_lowestmode,file_1$lat_lowestmode), data=file_1) #Converting dataframe to spatial points dataframe
    proj4string(file_1) = CRS("+init=epsg:4326")
    instudyarea <- !is.na(over(file_1, study_area)$id)
    if (sum(instudyarea) > 0){
      file_1 = file_1[instudyarea,]  #crop(file_1, study_area) #Neal originally cropped to the study area whereas I'm using the intersection. Crop throws an error if no points intersect the study area.
      file_1 = as.data.frame(file_1) #Converting cropped spatial points dataframe back to dataframe for further processing
      setwd(gedi_processed_output_2B) #Set the working directory to the temporary processing folder
      fwrite(file_1, file = out_name, sep = ",")
      rm(file_1) #remove file_1 for next iteration
      gc() #clear R memory for deleted files 
      
    } else {
      setwd(gedi_processed_output_2B)
      empty_df <- data.frame()
      write.csv(empty_df, file = out_name)
      rm(file_1)
      gc
    }   
  }
  
  if(num_rows < 1){
    setwd(gedi_processed_output_2B)
    empty_df <- data.frame()
    write.csv(empty_df, file = out_name)
    rm(file_1)
    gc
  }
  
  granule_name
  
}

stopCluster(cl)

### Setting up an empty storage dataframe for appending all data to
storage_dataframe_2B <- data.frame()

### Generating a list of processed 2B .csv granule files
setwd(gedi_processed_output_2B)
list <- list.files(setwd(gedi_processed_output_2B), pattern = "csv$")
list

### For loop to bind list of 2B granule csv files together
for(i in seq(list)) {
  csv_id <- list[[i]]
  granule_data <- fread(csv_id)
  num_row <- nrow(granule_data)
  if (num_row > 1){
    storage_dataframe_2B <- rbind(storage_dataframe_2B, granule_data)
    rm(granule_data)
    gc()
  }
  
}

### Write out the 2B merged GEDI data to a .csv file, remove temp directory
storage_dataframe_2B
setwd(gedi_processed_output)
fwrite(storage_dataframe_2B, "GEDI_2B_raw.csv", row.names=FALSE)
unlink(gedi_processed_output_2B, recursive = TRUE) #Remove the files in the temporary directory

###____________________________________________________________________________###
### Reading in the two temporary un-filtered storage data frames into R Memory ###
###____________________________________________________________________________###

### Read in the 2A dataframe, subset to remove unneeded columns, check headers and rows and file size in memory
setwd(gedi_processed_output)
storage_dataframe_2A <- fread("GEDI_2A_raw.csv")
storage_dataframe_2A
storage_dataframe_2A <- subset(storage_dataframe_2A, select = -c(lat_lowestmode, lon_lowestmode, beam, degrade_flag, quality_flag, sensitivity, solar_elevation, coords.x1, coords.x2) )
storage_dataframe_2A$shot_number = as.character(storage_dataframe_2A$shot_number)
gc()
print(object.size(x=lapply(ls(), get)), units="Mb")

### Read in the 2B dataframe, subset to remove unneeded columns, reformat column type, check headers and rows and file size in memory 
setwd(gedi_processed_output)
storage_dataframe_2B <- fread("GEDI_2B_raw.csv")
storage_dataframe_2B <- subset(storage_dataframe_2B, select = -c(coords.x1, coords.x2, rh100) )
storage_dataframe_2B$delta_time = as.Date(storage_dataframe_2B$delta_time)
storage_dataframe_2B$shot_number = as.character(storage_dataframe_2B$shot_number)
gc()
print(object.size(x=lapply(ls(), get)), units="Mb")

###___________________________________________________________________###
### Merging the subset 2A and 2B data frames into a single data frame ###
###___________________________________________________________________###

### Checking the 2A and 2B dataframes
names(storage_dataframe_2A)
names(storage_dataframe_2B)
nrow_2A<-nrow(storage_dataframe_2A)
nrow_2B<-nrow(storage_dataframe_2B)

### Merging the two dataframes together conditionally, matching based on shot_number column and dropping any mismatches
gedi_data_2A_2B <- merge(storage_dataframe_2B, storage_dataframe_2A, by="shot_number")
nrow(gedi_data_2A_2B)
gedi_data_2A_2B <- na.omit(gedi_data_2A_2B)
nrow(gedi_data_2A_2B)
rm(storage_dataframe_2A, storage_dataframe_2B)
gc()

names(gedi_data_2A_2B)

### Calculate the number of mismatched shots that were dropped, export informational result to text file
diff_2A_2B <- nrow_2A - nrow_2B
diff_2A_2B <- abs(diff_2A_2B)
setwd(gedi_processing_stats_output)
cat(diff_2A_2B, "mismatched shots dropped when merging spatially subset 2A and 2B Files", file="01_Number of Mismatched GEDI Shots.txt")

### Calculate the final number of shots after conditional filtering, export informational result to text file
filtered_gedi_data_rows <- nrow(gedi_data_2A_2B)
setwd(gedi_processing_stats_output)
cat(filtered_gedi_data_rows, "GEDI shot observations after conditional filtering", file="02_Number of Final Filtered GEDI Shots.txt")

###________________________________________________________________________________________###
### Converting merged 2A_2B dataframe to spatial points, exporting csv and ESRI shapefiles ###
###________________________________________________________________________________________###

### Generating a custom output name for the shapefile object and csv object
date_start <- noquote(as.character(daterange[1]))
date_end <- noquote((as.character(daterange[2])))
out_basename <- paste0("GEDI_2A_2B_merged_filtered_", date_start, "_", date_end)
out_gpkg_name <- paste0(out_basename, ".gpkg")
out_parquet_name <- paste0(out_basename, ".parquet")
out_csv_name <- paste0(out_basename, ".csv")

### Converting to spatial points dataframe
gedi_data_2A_2B_spdf <- SpatialPointsDataFrame(cbind(gedi_data_2A_2B$lon_lowestmode,gedi_data_2A_2B$lat_lowestmode), data=gedi_data_2A_2B)
proj4string(gedi_data_2A_2B_spdf) = CRS("+init=epsg:4326")

### writing out CSV file for storage
gedi_data_2A_2B = as.data.frame(gedi_data_2A_2B_spdf)
setwd(gedi_processed_output)
fwrite(gedi_data_2A_2B, file = out_csv_name, sep = ",")
rm(gedi_data_2A_2B)
gc()

## Writing out the simple features file as a geopackage
out_sf <- st_as_sf(gedi_data_2A_2B_spdf)
setwd(gedi_processed_output_shp)
st_write(out_sf, out_gpkg_name, append=FALSE)

## Write out as an Apache Paraquet file which can be efficiently read into geopandas
library(sfarrow)
st_write_parquet(out_sf, out_parquet_name)

