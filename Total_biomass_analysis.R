##############################
### Total biomass analysis ###
##############################

rm(list = ls())

require(dplyr)
require(mgcv)
require(raster)
require(terra)



####################
## Data wrangling ##
####################

# Loading data #
data <- as.data.frame(read.csv("Myctobase/event_edited2.csv", header = TRUE, stringsAsFactors = F))
str(data)
length(unique(data$eventID))

og_group <- as.data.frame(read.csv("Myctobase/groupOccurrence.csv", header = TRUE, stringsAsFactors = F))
group <- og_group[c("eventID", "family", "scientificName", "organismQuantity", "individualCount")]
str(group)
length(unique(group$eventID))

og_ind <- as.data.frame(read.csv("Myctobase/individualOccurrence.csv", header = TRUE, stringsAsFactors = F))
str(og_ind)
length(unique(og_ind$eventID))



# Making total abundance data of study species #
species <- group[which(group$scientificName %in% c("Electrona antarctica", 
                                                   "Electrona carlsbergi", 
                                                   "Gymnoscopelus braueri", 
                                                   "Gymnoscopelus fraseri", 
                                                   "Gymnoscopelus nicholsi", 
                                                   "Krefftichthys anderssoni", 
                                                   "Protomyctophum bolini",
                                                   "Protomyctophum tenisoni")),]

species_sum <- species %>%
  group_by(eventID) %>%
  summarise(sum = sum(individualCount))
species_sum

new <- merge(data, species_sum, by = "eventID", all.x = TRUE, all.y = FALSE)

new["sum"][is.na(new["sum"])] <- 0

new$CPUE <- (new$sum/new$volume)*1000
new$logCPUE <- log(new$CPUE + 1)



#####################
## All species GAM ## 
#####################

knots <- list(diel_num = c(0.5,1.5,2.5,3.5,4.5))
knots

# Mean raw abundance 
raw <- new %>% group_by(lat, diel_num, depth) %>% 
  summarise(mean = mean(CPUE), sum = sum(sum), count = n())

# Crop by latitude and net 
new <- new[new$lat >= -65,]
min(new$lat)
max(new$lat)

new_IKMT <- new[new$netType != "IKMT",]
dim(new_IKMT)


# Modelling
All_model1 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                        s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                      select = TRUE, data = new_IKMT,
                      family = "tw", knots = knots, method = "ML")

summary(All_model1)
plot(All_model1, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)
gam.check(All_model1)


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(All_model1, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_SE <- sum(day_predict_data$SE)
day_abundance_SE


# Nighttime abundance
night_predict_data <- data.frame(diel_num = 3.5, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(All_model1, night_predict_data, 
                             type = "response", se.fit = T)
night_predict_data <- cbind(night_predict, night_predict_data)
night_predict_data
night_predict_data$night_predict <- exp(night_predict_data$fit) -1
night_predict_data$SE <- exp(night_predict_data$se.fit) -1

night_predict_data$upr <- night_predict_data$night_predict + (2 * night_predict_data$SE)
night_predict_data$lwr <- night_predict_data$night_predict - (2 * night_predict_data$SE)
night_predict_data

night_abundance <- sum(night_predict_data$night_predict)
night_abundance
night_abundance_SE <- sum(night_predict_data$SE)
night_abundance_SE


# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
inflated_day_abundance <- sum(day_predict_data$inflated_day)
inflated_day_abundance


# Intercept
day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)
night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)

which.mins <- function(x, mins=3) {
  head(order(x), mins)
}

intercept <- day_prop$x[which.mins(abs(day_prop$y - night_prop$y))][1]
intercept


# Proportion of DVM population 
DVM_proportion <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
DVM_proportion



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


# Calculate area 
ocean_area <- cellSize(transformed, unit = "m", mask = T)
ocean_area
total_area <- sum(values(ocean_area), na.rm = T)
total_area/10^12


# Total myctophid abundance 
total_abundance <- (night_abundance/10^3)*total_area
total_abundance/10^12

total_abundance_SE <- (night_abundance_SE/10^3)*total_area
total_abundance_SE/10^12



##########################################
# Proportion of species in the community #
##########################################

# Making individual species abundance data #
table <- table(group$scientificName) 
table

order <- order(table, decreasing = TRUE)
table[order]

species_choice <- c(1,3,5,7,8,10,11,12)
sp_order <- names(table[order[species_choice]])

species <- list()

for(i in species_choice){
  name <- names(table[order[i]])
  print(name)
  
  dataframe <- data.frame(group[(group$scientificName == name),])
  dataframe <- dataframe[rowSums(is.na(dataframe)) != ncol(dataframe),]
  
  species[[i]] <- dataframe
}

# Combining data 
for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  print(length(species[[i]]$eventID))
  combined <- merge(data[data$lat >= -65,], species[[i]], by = "eventID", all.x = TRUE, all.y = FALSE)
  
  print(length(combined$eventID)) 
  
  combined$scientificName <- name
  
  species[[i]] <- combined
  
}  

# Filling 0 counts 
for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  data <- species[[i]]
  
  print(length(data$eventID))
  print(sum(is.na(data$individualCount)))
  
  data[["individualCount"]][is.na(data[["individualCount"]])] <- 0
  data$CPUE <- (data$individualCount/data$volume) *1000
  data$logCPUE <- log(data$CPUE + 1)
  
  species[[i]] <- data
  
}



# Calculating mean abundance per latitude #

all_catch <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), 
                      c("lat", "mean", "species"))
all_catch

name_order <- vector()

for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  name_order <- c(name_order, name)
  
  data <- species[[i]]
  
  data <- data[, c("CPUE", "lat")]
  print(plot(CPUE ~ lat, data))
  
  data <- data %>%
    group_by(lat) %>%
    summarise(mean = mean(CPUE))
  data$scientificName <- name
  
  print(head(data))
  all_catch <- rbind(all_catch, data)
  
}


# Average abundance across latitudes
proportion <- all_catch %>%
  group_by(scientificName) %>% summarise(mean = mean(mean))
proportion


# Proportion within the community 
total <- sum(proportion$mean)
total

proportion$freq <- proportion$mean/total
proportion
sum(proportion$freq)


# Abundance of individual species 
proportion$abundance <- (proportion$freq*night_abundance)/1000
proportion$abundance_SE <- (proportion$freq*night_abundance_SE)/1000
proportion



##################################
### Average individual biomass ###
##################################

# Extracting size data for study species 
ind_species <- og_ind[og_ind$scientificName == "Electrona antarctica" |
                  og_ind$scientificName == "Electrona carlsbergi"|
                  og_ind$scientificName == "Gymnoscopelus braueri"|
                  og_ind$scientificName == "Gymnoscopelus fraseri"|
                  og_ind$scientificName == "Gymnoscopelus nicholsi"|
                  og_ind$scientificName == "Krefftichthys anderssoni"|
                  og_ind$scientificName == "Protomyctophum bolini"|
                  og_ind$scientificName == "Protomyctophum tenisoni",]
head(ind_species)
species_list <- unique(ind_species$scientificName)
species_list


# Filling data from literature 
ind_species$standard_length[ind_species$standard_length == "SC 1"] <- 10
ind_species$standard_length[ind_species$standard_length == "SC  1"] <- 10
ind_species$standard_length[ind_species$standard_length == "SC 2"] <- 30
ind_species$standard_length[ind_species$standard_length == "SC 3"] <- 50

ind_species$standard_length <- as.numeric(ind_species$standard_length)


# Mean length and ww from empirical data 
mean_size <- ind_species %>%
  group_by(scientificName) %>%
  summarise(mean_legnth = round(mean(standard_length, na.rm = T), digits = 1),
            mean_weight = round(mean(weight/individualCount, na.rm = T), digits = 1),
            n.length = sum(!is.na(standard_length)),
            n.weight = sum(!is.na(weight)))
mean_size


# Modelling average fish weight #
model_output <- data.frame(scientificName = NA, 
                           lat.coef = NA, depth.coef = NA, pca.coef = NA,
                           lat = NA, depth = NA, pca = NA, average_fish = NA, 
                           a = NA, b = NA,
                           average_weight = NA)
model_output

species_list <- unique(ind_species$scientificName)
species_list

for(i in species_list){
  sp_length <- ind_species[ind_species$scientificName == i,]
  sp_length$ww <- sp_length$weight/sp_length$individualCount
  
  print(i)  
  print(colSums(!is.na(sp_length[,c("standard_length", "ww")])))

  sp_length$logww <- log(sp_length$ww)
  sp_length$logSL <- log(sp_length$standard_length)
  
  # Average fish length 
  sp_all <- merge(sp_length, new, by = "eventID", all.x = T, all.y = F)
  hist(sp_all$standard_length, main = i)
  length_model <- lm(standard_length ~ lat + depth + pca, sp_all)
  print((unique(sp_all$lat)))
  print(median(unique(sp_all$lat), na.rm = T))
  print(sort(unique(sp_all$depth)))
  print(median(unique(sp_all$depth), na.rm = T))
  length_coefficient <- data.frame(lat.coef = length_model$coefficients[[2]],
                                   depth.coef = length_model$coefficients[[3]],
                                   pca.coef = length_model$coefficients[[4]])
  average_point <- data.frame(lat = median(unique(sp_all$lat), na.rm = T), 
                              depth = round(median(unique(sp_all$depth), na.rm = T), digits = 0), 
                              pca = 0)
  average_fish <- predict(length_model, average_point)
  print(average_fish)
  
  # Weight ~ SL allometric relationship 
  allometric_model <- lm(logww ~ logSL, sp_length)
  weight_coefficient <- data.frame(a = round(allometric_model$coefficients[[1]], digits = 1),
                                   b = round(allometric_model$coefficients[[2]], digits = 1))
  
  average_weight <- predict(allometric_model, data.frame(logSL = log(average_fish)))
  average_weight <- exp(average_weight)
  print(average_weight)
  
  average_weight <- cbind(scientificName = i, length_coefficient, average_point, average_fish, 
                          weight_coefficient, average_weight)
  print(average_weight)
  
  model_output <- rbind(model_output, average_weight)
}

model_output <- na.omit(model_output)
model_output
model_output$average_fish <- round(model_output$average_fish, digits = 1)
model_output$average_weight <- round(model_output$average_weight, digits = 1)


# Export size summary table 
size_table <- merge(mean_size, model_output, by = "scientificName")
size_table
write.csv(size_table, "Size_summary.csv", row.names = F)



#####################
### Total biomass ###
#####################

predicted_size <- model_output[,c("scientificName", "average_fish", "average_weight")]
predicted_size
names(predicted_size) <- c("scientificName", "SL", "weight")

total_weight <- merge(predicted_size, proportion, by = "scientificName", all.x = T, all.y = T)
total_weight

# Biomass per metre squared per species 
total_weight$abundance_weight <- total_weight$abundance*as.numeric(total_weight$weight)
total_weight$abundance_weight_SE <- total_weight$abundance_SE*as.numeric(total_weight$weight)
total_weight

# Total abundance per species 
total_weight$total_abundance <- (total_weight$abundance)*total_area
total_weight$total_abundance_SE <- (total_weight$abundance_SE)*total_area
total_weight

# Total biomass per species 
total_weight$total_biomass <- total_weight$total_abundance*as.numeric(total_weight$weight)
total_weight$total_biomass_SE <- total_weight$total_abundance_SE*as.numeric(total_weight$weight)
total_weight

# 95% CI for total biomass
(sum(total_weight$total_biomass) + sum(2*total_weight$total_biomass_SE))/(10^12) 
(sum(total_weight$total_biomass) - sum(2*total_weight$total_biomass_SE))/(10^12)

#DVM biomass 
sum(total_weight$total_biomass)*DVM_proportion/(10^12)


# All species sum #
total_weight <- total_weight %>%
  mutate(across(2:last_col(), ~ as.numeric(.))) 

total_weight[nrow(total_weight) + 1,] <- c("All species",
                                          NA, NA,
                                          sum(total_weight$mean), NA,
                                          sum(total_weight$abundance),
                                          sum(total_weight$abundance_SE),
                                          sum(total_weight$abundance_weight),
                                          sum(total_weight$abundance_weight_SE),
                                          sum(total_weight$total_abundance),
                                          sum(total_weight$total_abundance_SE),
                                          sum(total_weight$total_biomass),
                                          sum(total_weight$total_biomass_SE))
total_weight


# Export table 
total_weight_final <- total_weight %>%
  mutate(across(2:last_col(), as.numeric)) %>%
  mutate(across(10:last_col(), ~ ./10^12)) %>%
  mutate(across(2:last_col(), round, 3)) %>%
  mutate(across(c(2,3,5,10:last_col()), round, 2))
total_weight_final

total_weight_final$abundance <- paste0(total_weight_final$abundance, " (",
                                       total_weight_final$abundance_SE, ")")
total_weight_final$abundance_weight <- paste0(total_weight_final$abundance_weight, " (",
                                              total_weight_final$abundance_weight_SE, ")")
total_weight_final$total_abundance <- paste0(total_weight_final$total_abundance, " (",
                                              total_weight_final$total_abundance_SE, ")")
total_weight_final$total_biomass <- paste0(total_weight_final$total_biomass, " (",
                                              total_weight_final$total_biomass_SE, ")")
total_weight_final <- total_weight_final %>%
  arrange(factor(scientificName, levels = sp_order))
total_weight_final

write.csv(total_weight_final, "Total_biomass_estimate.csv", row.names = F)


