###############################
### Data wrangling for QGIS ###
###############################

rm(list = ls())

require(dplyr)
require(raster)
require(terra)



# Loading data #
data <- as.data.frame(read.csv("Myctobase/event_edited2.csv", header = TRUE, stringsAsFactors = F))
str(data)
length(unique(data$eventID))

group <- as.data.frame(read.csv("Myctobase/groupOccurrence.csv", header = TRUE, stringsAsFactors = F))
dim(size)

size <- as.data.frame(read.csv("Myctobase/individualOccurrence.csv", header = TRUE, stringsAsFactors = F))
dim(size)



##############
## Net type ##
##############

csv <- data[,c( "lon", "lat", 
                "netType", "cruiseCode")]
head(csv)

size_data <- unique(size$eventID)
size_data

size_data <- unique(data[data$eventID %in% size_data,]$cruiseCode)
size_data

csv$size <- NA
csv[csv$cruiseCode %in% size_data,]$size <- "y"
csv[is.na(csv$size),]$size <- "n"
head(csv)

write.csv(csv, "Net type.csv", row.names = FALSE)



###################################
### Southern Ocean surface area ###
###################################

# Load netcdf data
raster_ocean <- terra::rast("GEBCO_05_Dec_2023_c97a092c1373/gebco_2023_n-40.0_s-80.0_w-180.0_e180.0.nc")
raster_ocean
res(raster_ocean)


# Define crs
crs(raster_ocean) <- "EPSG:4326"
crs(raster_ocean)
terra::plot(raster_ocean)
res(raster_ocean)


# Lower resolution to 0.5x0.5
ocean_rough <- terra::aggregate(raster_ocean, fact=120, fun = mean)
res(ocean_rough)
terra::plot(ocean_rough)
ocean_rough


# Crop Southern Ocean  
e <- as(extent(-180, 180, -65, -48), 'SpatialPolygons')
crs(e) <- "EPSG:4326"
ocean_cropped <- crop(ocean_rough, e)
ocean_cropped
plot(ocean_cropped)


# Mask land region and ocean shallower than 1000m
ocean_only <- terra::ifel(ocean_cropped >= 0, NA, ocean_cropped)
plot(ocean_only)
ocean_only <- terra::ifel(ocean_cropped > -1000, NA, 1)
ocean_only
plot(ocean_only)


# Reproject to equal area mapping 
transformed <- terra::project(ocean_only, "epsg:6932")
transformed
crs(transformed)
plot(transformed)
linearUnits(transformed)


# Exporting as polygon 
polygon <- as.polygons(transformed)
plot(polygon)
writeVector(polygon, "Ocean_extent", filetype = 'ESRI Shapefile', overwrite = T)



############
# CPUE map #
############

table <- table(og_group$scientificName) 
table

order <- order(table, decreasing = TRUE)
table[order]

species_choice <- c(1,3,5,7,8,10,11,12)
table[order[species_choice]]

species <- list()

sp_list <- vector()

for(i in species_choice){
  name <- names(table[order[i]])
  sp_list <- c(sp_list, name)
  print(name)
  
  dataframe <- data.frame(og_group[(og_group$scientificName == name),])
  dataframe <- dataframe[rowSums(is.na(dataframe)) != ncol(dataframe),]
  
  species[[i]] <- dataframe
}

for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  print(length(species[[i]]$eventID))
  new <- merge(data, species[[i]], by = "eventID", all.x = TRUE, all.y = FALSE)
  
  print(length(new$eventID)) 
  
  new$scientificName <- name
  
  species[[i]] <- new
  
}  

for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  sub_data <- species[[i]]
  
  print(length(sub_data$eventID))
  print(sum(is.na(sub_data$individualCount)))
  
  sub_data[is.na(sub_data$individualCount),]$individualCount <- 0
  sub_data$CPUE <- as.integer(sub_data$individualCount/sub_data$volume *10000)
  sub_data$logCPUE <- log(sub_data$CPUE + 1)
  
  species[[i]] <- sub_data
  
}


for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  sub_data <- species[[i]]
  
  sp_csv <- sub_data[,c( "lon", "lat", 
                         "individualCount", "CPUE")]
  
  sp_csv <- sp_csv %>%
    group_by(lat, lon) %>%
    summarise(mean_CPUE = mean(CPUE),
              totalCount = sum(individualCount))
  
  sp_csv$presence <- NA
  sp_csv[sp_csv$totalCount > 0, ]$presence <- "y"
  sp_csv[sp_csv$totalCount == 0, ]$presence <- "n"
  
  write.csv(sp_csv, paste0(name, "_CPUE.csv"), row.names = F)
  
}
  