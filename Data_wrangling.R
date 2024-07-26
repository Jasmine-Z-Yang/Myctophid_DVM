######################
### Data wrangling ###
######################

rm(list = ls())

require(dplyr)
require(factoextra)
require(ggplot2)
require(Hmisc)
require(lubridate)
require(patchwork)



# Loading data #
og_event <- as.data.frame(read.csv("Myctobase/event_edited.csv", header = TRUE, stringsAsFactors = F))
str(og_event)
length(unique(og_event$eventID))

event <- og_event[c("eventID", "cruiseCode", "sector", "institutionCode",
                    "meanTowSpeed", 
                    "start_latitude", "end_latitude", "start_longitude", "end_longitude", 
                    "diel", "solarPosition", "haulDuration", "haulType",
                    "minimumDepthInMeters", "maximumDepthInMeters",
                    "validation", "netType", "samplingProtocol", "netSurfaceArea", "codendMesh",
                    "haulDuration", "start_eventTime", "end_eventTime",
                    "volumeFiltered", "volumeFiltered2")] 
str(event)



# Removing invalid data #
data <- event[event$validation == "1",]



# Excluding target trawls #
data$haulType <- replace(data$haulType, is.na(data$haulType), "Unknown")

data <- data[data$haulType != "Target", ]
data <- data[data$haulType != "Failed", ]



# Latitude simplification #

# Checking latitude range doesn't exceed 1
data$lat_diff <- abs(data$start_latitude - data$end_latitude) 
sum((data$lat_diff > 1) & (!is.na(data$lat_diff)))

# Rounding latitude to whole number
data$lat <- trunc(rowMeans(data[,c("start_latitude", "end_latitude")], na.rm = TRUE))
sum(is.na(data$lat))
data <- data[!is.na(data$lat), ]

max(data$lat)
min(data$lat)



# Longitude simplification #

# Checking longitude range doesn't exceed 1
data$lon_diff <- abs(data$start_longitude - data$end_longitude) 
head(data$lon_diff)
sum((data$lon_diff > 1) & (!is.na(data$lon_diff)))

# Rounding longitude to whole number
data$lon <- trunc(rowMeans(data[,c("start_longitude", "end_longitude")], na.rm = TRUE))
sum(is.na(data$lon))

max(data$lon)
min(data$lon)



# Depth selection #

# Remove data without depth info 
dim(data)
data <- data[!is.na(data$minimumDepthInMeters | data$maximumDepthInMeters ),]
dim(data) 

# Making one column for depth 
data$depth <- NA

i <- 1
j <- 1

for(row in 1:nrow(data)){
  if(is.na(data$maximumDepthInMeters[row])){
    data$depth[row] = data$minimumDepthInMeters[row]
   }
  else if(is.na(data$minimumDepthInMeters[row])){
    data$depth[row] = data$maximumDepthInMeters[row]
    j <- j+1}
  else {
    difference = data$maximumDepthInMeters[row] - data$minimumDepthInMeters[row]
    if(difference <= 300){
      data$depth[row] = median(c(data$minimumDepthInMeters[row], data$maximumDepthInMeters[row])) 
      i <- i+1} else {data$depth[row] = NA}
  }
}

data <- data[!is.na(data$depth),]
data$depth <- as.numeric(data$depth)
hist(data$depth)



# Combine water volume columns #
data <- data[!is.na(data$volumeFiltered | data$volumeFiltered2),]

data$volumeFiltered <- as.numeric(data$volumeFiltered)
data$volumeFiltered2 <- as.numeric(data$volumeFiltered2)
data$volume <- rowSums(data[,c("volumeFiltered", "volumeFiltered2")], na.rm = TRUE)



# Filling missing diel values #
for(i in 1:nrow(data)){
  if(is.na(data$diel[i])){
    time = hour(dmy_hm(data$start_eventTime[i]))
    if(time >= 12){
      if(-12 <= data$solarPosition[i] & data$solarPosition[i] <= 12){
        data$diel[i] <- "dusk"
      } else if(data$solarPosition[i] < 0){
        data$diel[i] <- "night"
      } else if(data$solarPosition[i] > 0){
        data$diel[i] <- "day"
      }
    } else {
      if(-12 <= data$solarPosition[i] & data$solarPosition[i] <= 12){
        data$diel[i] <- "dawn"
      } else if(data$solarPosition[i] < 0){
        data$diel[i] <- "night"
      } else if(data$solarPosition[i] > 0){
        data$diel[i] <- "day"
      }
    }
  }
}



# Making numerical diel column #
data$diel_num <- NA 

for(i in 1:nrow(data)){
  if(data$diel[i] == "dawn"){
    data$diel_num[i] <- 1
  } else if(data$diel[i] == "day"){
    data$diel_num[i] <- 2
  } else if(data$diel[i] == "dusk"){
    data$diel_num[i] <- 3
  } else if(data$diel[i] == "night"){
    data$diel_num[i] <- 4
  } else {data$diel[i] <- NA}
}



# Editing variable names #

# Grouping cruises
data[grep("AMLR", data$cruiseCode),]$cruiseCode <- "AMLR"
data[grep("Ichtyo", data$cruiseCode),]$cruiseCode <- "Ichtyo"
data[grep("ichtyo", data$cruiseCode),]$cruiseCode <- "Ichtyo"

# Naming IYGPT nets by size
mean(as.numeric(data[data$cruiseCode == "HIPPIES" & data$netType == "IYGPT",]$netSurfaceArea), na.rm = T)

for(i in 1:nrow(data)){
  if(data$netType[i] == "MOHT opening/closing nets"){
    data$netType[i] <- "MOHT"
  } else if(data$netType[i] == "RMT 8"){
    data$netType[i] <- "RMT8"
  } else if(data$netType[i] == "RMT 1"){
    data$netType[i] <- "RMT1"
  } else if(data$netType[i] == "IYGPT"){
    print(data$netType[i])
    if(data$cruiseCode[i] == "Ichtyo"){
      print(data$cruiseCode[i])
      data$netType[i] <- "IYGPT (66)"
    } else if(data$cruiseCode[i] == "HIPPIES"){
      data$netType[i] <- "IYGPT (171)"
    }
  } else if(data$netType[i] == "IYGPT with MIDOC"){
    if(data$cruiseCode[i] == "SS1999_01"){
      data$netType[i] <- "IYGPT with MIDOC (57)"
    } else if(data$cruiseCode[i] == "KAXIS"){
      data$netType[i] <- "IYGPT with MIDOC (188)"
    }
  } else{next}
}

Cruise <- sort(unique(data$cruiseCode))
print(Cruise)

# Filling missing data 
data$meanTowSpeed[data$cruiseCode == "brokeeast"] <- 2
data$meanTowSpeed[data$cruiseCode == "SS1999_01"] <- (1.5+3)/2
data$meanTowSpeed[data$cruiseCode == "JB11"] <- 3

data$codendMesh[data$cruiseCode == "HIPPIES" & data$netType == "IYGPT (171)"] <- 10
data$codendMesh[data$cruiseCode == "HIPPIES" & data$netType == "RMT8"] <- 0.85
data$codendMesh[data$netType == "MOHT"] <- 1.95



##########################
# Volume anomaly removal #
##########################

colors <- c("normal" = "black", "outlier" = "red")

# AMLR
Cruise[1]
cruise <- data[data$cruiseCode == Cruise[1],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

for(i in 1:nrow(cruise)){
  if(cruise$z_score[i] < 3){
    cruise$color[i] <- "red"
  } else if(cruise$z_score[i] > 3){
    cruise$color[i] <- "red"
  } else {cruise$color[i] <- "grey"}
}

plot1 <- ggplot(data = cruise, aes(x = volume)) +
  geom_histogram(aes(fill = cut(z_score, c(-Inf, -3,3, Inf)))) +
  scale_fill_manual(values = c("red", "grey", "red")) +
  xlab(bquote("Volume "(m^3))) + 
  geom_vline(aes(xintercept = mean(volume)), col = "blue") +
  ggtitle(paste0("(A) ", Cruise[1], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "none")
plot1

sum(cruise$z_score > 3)
sum(cruise$z_score < -3)
min(cruise$volume)

cruise[cruise$z_score > 3,]$volume <- mean(cruise$volume)
cruise[cruise$z_score < -3,]$volume <- mean(cruise$volume)

new_data <- cruise


# brokeeast
Cruise[2]
cruise <- data[data$cruiseCode == Cruise[2],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)$sigma

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot2 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(B) ", Cruise[2], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot2

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

new_data <- bind_rows(new_data, cruise)


# brokewest
Cruise[3]
cruise <- data[data$cruiseCode == Cruise[3],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot3 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(C) ", Cruise[3], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot3

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

new_data <- bind_rows(new_data, cruise)


# HIPPIES
Cruise[4]
cruise <- data[data$cruiseCode == Cruise[4],]
unique(cruise$netType)
head(cruise)

# RMT8
cruise1 <- cruise[cruise$netType == "RMT8",]
head(cruise1)

cruise1$z_score <- (cruise1$volume-mean(cruise1$volume))/sd(cruise1$volume)

model <- lm(data = cruise1, volume~haulDuration)
summary(model)

cruise1$predict <- predict.lm(model, cruise1)
cruise1$max <- cruise1$predict + 3*summary(model)$sigma
cruise1$min <- cruise1$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise1)){
  if(cruise1$volume[i] < cruise1$min[i]){
    cruise1$color[i] <- "outlier"
  } else if(cruise1$volume[i] > cruise1$max[i]){
    cruise1$color[i] <- "outlier"
  } else {cruise1$color[i] <- "normal"}
}

plot4 <- ggplot(data = cruise1, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(D) ", Cruise[4], ": ", cruise1$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot4

which(cruise1$volume > cruise1$max)
which(cruise1$volume < cruise1$min)

new_data <- bind_rows(new_data, cruise1)


# IYGPT - 171m
cruise2 <- cruise[cruise$netType == "IYGPT (171)",]
head(cruise2)

cruise2$z_score <- (cruise2$volume-mean(cruise2$volume))/sd(cruise2$volume)

model <- lm(data = cruise2, volume~haulDuration)
summary(model)

cruise2$predict <- predict.lm(model, cruise2)
cruise2$max <- cruise2$predict + 3*summary(model)$sigma
cruise2$min <- cruise2$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise2)){
  if(cruise2$volume[i] < cruise2$min[i]){
    cruise2$color[i] <- "outlier"
  } else if(cruise2$volume[i] > cruise2$max[i]){
    cruise2$color[i] <- "outlier"
  } else {cruise2$color[i] <- "normal"}
}

plot5 <- ggplot(data = cruise2, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(E) ", Cruise[4], ": ", cruise2$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot5

which(cruise2$volume > cruise2$max)
which(cruise2$volume < cruise2$min)

new_data <- bind_rows(new_data, cruise2)


# Ichtyo
Cruise[5]
cruise <- data[data$cruiseCode == Cruise[5],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot6 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(F) ", Cruise[5], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot6

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

cruise[which(cruise$volume > cruise$max),]
cruise[which(cruise$volume < cruise$min),]

cruise[which(cruise$volume > cruise$max),]$volume <- cruise[which(cruise$volume > cruise$max),]$predict
cruise[which(cruise$volume < cruise$min),]$volume <- cruise[which(cruise$volume < cruise$min),]$predict

new_data <- bind_rows(new_data, cruise)


# JB11
Cruise[6]
cruise <- data[data$cruiseCode == Cruise[6],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot7 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(G) ", Cruise[6], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot7

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

new_data <- bind_rows(new_data, cruise)


# JR100
Cruise[7]
cruise <- data[data$cruiseCode == Cruise[7],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot8 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(H) ", Cruise[7], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot8

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

new_data <- bind_rows(new_data, cruise)


# JR161
Cruise[8]
cruise <- data[data$cruiseCode == Cruise[8],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot9 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(I) ", Cruise[8], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot9

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

new_data <- bind_rows(new_data, cruise)


# JR177
Cruise[9]
cruise <- data[data$cruiseCode == Cruise[9],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot10 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(J) ", Cruise[9], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot10

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

new_data <- bind_rows(new_data, cruise)


# JR200
Cruise[10]
cruise <- data[data$cruiseCode == Cruise[10],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot11 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(K) ", Cruise[10], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot11

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

new_data <- bind_rows(new_data, cruise)


# KAXIS
Cruise[11]
cruise <- data[data$cruiseCode == Cruise[11],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

cruise$max <- mean(cruise$volume) + 3*sd(cruise$volume)
cruise$min <- mean(cruise$volume) - 3*sd(cruise$volume)

for(i in 1:nrow(cruise)){
  if(cruise$z_score[i] < 3){
    cruise$color[i] <- "red"
  } else if(cruise$z_score[i] > 3){
    cruise$color[i] <- "red"
  } else {cruise$color[i] <- "grey"}
}

plot12 <- ggplot(data = cruise, aes(x = volume)) +
  geom_histogram(aes(fill = cut(z_score, c(-Inf, -3,3, Inf)))) +
  scale_fill_manual(values = c("grey", "red")) +
  xlab(bquote("Volume "(m^3))) + 
  geom_vline(aes(xintercept = mean(volume)), col = "blue") +
  ggtitle(paste0("(L) ", Cruise[11], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "none")
plot12

sum(cruise$z_score > 3)
sum(cruise$z_score < -3)

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

cruise[which(cruise$z_score > 3),]

cruise[cruise$z_score > 3,]$volume <- mean(cruise$volume)

new_data <- bind_rows(new_data, cruise)


# PS65
Cruise[12]
cruise <- data[data$cruiseCode == Cruise[12],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot13 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(N) ", Cruise[12], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot13

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

cruise[which(cruise$volume > cruise$max),]$volume <- cruise[which(cruise$volume > cruise$max),]$predict

new_data <- bind_rows(new_data, cruise)


# PS69_4
Cruise[13]
cruise <- data[data$cruiseCode == Cruise[13],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma
head(cruise)

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot14 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(O) ", Cruise[13], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot14

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

cruise[which(cruise$volume > cruise$max),]$volume <- cruise[which(cruise$volume > cruise$max),]$predict
cruise[which(cruise$volume < cruise$min),]$volume <- cruise[which(cruise$volume < cruise$min),]$predict

new_data <- bind_rows(new_data, cruise)


# PS69_6
Cruise[14]
cruise <- data[data$cruiseCode == Cruise[14],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot15 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(P) ", Cruise[14], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot15

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

cruise[which(cruise$volume < cruise$min),]$volume <- cruise[which(cruise$volume < cruise$min),]$predict

new_data <- bind_rows(new_data, cruise)


# SS1999_01
Cruise[15]
cruise <- data[data$cruiseCode == Cruise[15],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

cruise$max <- mean(cruise$volume) + 3*sd(cruise$volume)
cruise$min <- mean(cruise$volume) - 3*sd(cruise$volume)

for(i in 1:nrow(cruise)){
  if(cruise$z_score[i] < 3){
    cruise$color[i] <- "red"
  } else if(cruise$z_score[i] > 3){
    cruise$color[i] <- "red"
  } else {cruise$color[i] <- "grey"}
}

plot16 <- ggplot(data = cruise, aes(x = volume)) +
  geom_histogram(aes(fill = cut(z_score, c(-Inf, -3,3, Inf)))) +
  scale_fill_manual(values = c("red","grey")) +
  xlab(bquote("Volume "(m^3))) + 
  geom_vline(aes(xintercept = mean(volume)), col = "blue") +
  ggtitle(paste0("(Q) ", Cruise[15], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "none")
plot16

sum(cruise$z_score > 3)
sum(cruise$z_score < -3)

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

cruise[which(cruise$z_score < -3),]

cruise[cruise$z_score < -3,]$volume <- mean(cruise$volume)

new_data <- bind_rows(new_data, cruise)


# UM-18-08
Cruise[16]
cruise <- data[data$cruiseCode == Cruise[16],]
unique(cruise$netType)
head(cruise)

cruise$z_score <- (cruise$volume-mean(cruise$volume))/sd(cruise$volume)

model <- lm(data = cruise, volume~haulDuration)
summary(model)

cruise$predict <- predict.lm(model, cruise)
cruise$max <- cruise$predict + 3*summary(model)$sigma
cruise$min <- cruise$predict - 3*summary(model)$sigma
head(cruise)

for(i in 1:nrow(cruise)){
  if(cruise$volume[i] < cruise$min[i]){
    cruise$color[i] <- "outlier"
  } else if(cruise$volume[i] > cruise$max[i]){
    cruise$color[i] <- "outlier"
  } else {cruise$color[i] <- "normal"}
}

plot17 <- ggplot(data = cruise, aes(x = haulDuration, y = volume)) +
  geom_ribbon(aes(ymax = max, ymin = min), alpha = 0.1) +
  geom_point(aes(color = color)) +
  scale_color_manual(values = colors) +
  geom_smooth(method = "lm", se = F) +
  xlab("Haul duration (min)") + ylab(bquote("Volume "(m^3))) +
  ggtitle(paste0("(R) ", Cruise[16], ": ", cruise$netType[1])) +
  theme_classic() +
  scale_y_continuous(labels = scales::comma) +
  xlim(0,max(cruise$haulDuration)) +
  theme(legend.position = "none")
plot17

which(cruise$volume > cruise$max)
which(cruise$volume < cruise$min)

new_data <- bind_rows(new_data, cruise)


##############################
# Making net effect variable #
##############################

data_IKMT <- new_data[new_data$netType != "IKMT",]

data_IKMT$codendMesh <- as.numeric(data_IKMT$codendMesh)
data_IKMT$netSurfaceArea <- as.numeric(data_IKMT$netSurfaceArea)
data_IKMT$meanTowSpeed <- as.numeric(data_IKMT$meanTowSpeed)

correlation <- cbind(data_IKMT$meanTowSpeed, data_IKMT$netSurfaceArea, data_IKMT$codendMesh)
rcorr(correlation, type = "spearman")

pca_all <- prcomp(data = data_IKMT, ~ netSurfaceArea + codendMesh + meanTowSpeed, scale = T)
summary(pca_all)
plot(pca_all)
fviz_pca_biplot(pca_all)

axes_all <- predict(pca_all, newdata = data_IKMT)

pca <- data.frame(pca = axes_all[,"PC1"], eventID = data_IKMT$eventID)
head(pca)

data_pca <- merge(pca, new_data, by = "eventID", all.x = TRUE, all.y = TRUE)
head(data_pca)
dim(data_pca)


# Exporting data #

# Volume anomaly image 
Volume_plot <- plot1 + plot2 + plot3 + plot4 + plot5 +
  plot6 + plot7 + plot8 + plot9 + plot10 + plot11 + plot12 + plot13 +
  plot14 + plot15 + plot16 + plot17 + plot_layout(ncol = 6)
Volume_plot
# 3300W x 1500H (33 x 15)


# Post-quality filtration data 
write.csv(data_pca, "Myctobase/event_edited2.csv", row.names = F)


# Net summary info 
net_types <- new_data %>%
  group_by(netType, institutionCode) %>%
  summarise(net = paste0(round(mean(as.numeric(netSurfaceArea), na.rm = T), digits = 1),
                         " (", round(min(as.numeric(netSurfaceArea), na.rm = T), digits = 1),
                         " - ", round(max(as.numeric(netSurfaceArea), na.rm = T), digits = 1), ")"), 
            mesh = mean(as.numeric(codendMesh), na.rm = T), 
            speed = paste0(round(mean(as.numeric(meanTowSpeed), na.rm = T), digits = 3),
                           " (", round(min(as.numeric(meanTowSpeed), na.rm = T), digits = 1),
                           " - ", round(max(as.numeric(meanTowSpeed), na.rm = T), digits = 1), ")"), 
            lat = paste0(round(abs(max(start_latitude, na.rm = T)), digits = 0), "-", 
                         round(abs(min(start_latitude, na.rm = T)), digits = 0), "°S"),
            time = paste(format(min(as.Date(start_eventTime, format = "%d/%m/%Y"), na.rm = T), "%d/%m/%Y"), "-",
                         format(max(as.Date(start_eventTime, format = "%d/%m/%Y"), na.rm = T), "%d/%m/%Y")),
            count = n()) 
net_types
write.csv(net_types, "Net_summary.csv", row.names = F)
