###################################
### DVM analysis (light effect) ###
###################################

rm(list = ls())

require(dplyr)
require(geomtextpath)
require(ggplot2)
require(mgcv)
require(patchwork)
require(DescTools)



####################
## Data wrangling ##
####################

# Loading data #
data <- as.data.frame(read.csv("../Myctophid_DVM-main/Myctobase/event_edited2.csv", header = TRUE, stringsAsFactors = F))
str(data)
length(unique(data$eventID))

og_group <- as.data.frame(read.csv("Myctobase/groupOccurrence.csv", header = TRUE, stringsAsFactors = F))
str(og_group)
length(unique(og_group$eventID))



# Extracting data #
group <- og_group[c("eventID", "family", "scientificName", "organismQuantity", "individualCount")]
str(group)



# Choosing study species #

table <- table(group$scientificName) 
table

order <- order(table, decreasing = TRUE)
table[order]

species_choice <- c(1,3,5,7,8,10,11,12)
table[order[species_choice]]

species <- list()

for(i in species_choice){
  name <- names(table[order[i]])
  print(name)
  
  dataframe <- data.frame(group[(group$scientificName == name),])
  dataframe <- dataframe[rowSums(is.na(dataframe)) != ncol(dataframe),]
  
  species[[i]] <- dataframe
}



# Combining data #
for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  print(length(species[[i]]$eventID))
  new <- merge(data, species[[i]], by = "eventID", all.x = TRUE, all.y = FALSE)
  
  print(length(new$eventID)) 
  
  new$scientificName <- name
  
  species[[i]] <- new
  
}  


# Filling 0 counts #
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



# Distribution range function 
f <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


# Function to find intercept 
which.mins <- function(x, mins=3) {
  head(order(x), mins)
}


# Define knots 
knots <- list(diel_num = c(0.5,1.5,2.5,3.5,4.5))


# Result tables
Results_table <- data.frame(Species = NA,
                            median_night = NA,
                      median_day1 = NA,
                      median_day2 = NA,
                      median_day3 = NA,
                      median_day4 = NA,
                      Threshold1 = NA,
                      Threshold2 = NA,
                      Threshold3 = NA,
                      Threshold4 = NA,
                      proportion1 = NA,
                      proportion2 = NA,
                      proportion3 = NA,
                      proportion4 = NA)
Results_table



########################
# SO Light Attenuation #
########################

Depth <- 0:1000
Depth

Chl <- 3.5*exp(-((Depth - 30)^2)/(2*(22)^2))
Chl
which.max(Chl)
plot(Chl, Depth)

ifelse(Chl < 0.001, 
       K <- 0.0166 + 0.0325*exp(-0.014*(485-440)),
       K <- 0.0166 + 0.07242*(Chl^0.68955))
plot(K, Depth)
mean(K)

I0 <- 9.84*10^17
Light <- I0 * exp(-mean(K)*Depth)

Attenuation <- data.frame(Depth, K, Light)
head(Attenuation)

top <- Closest(Light, I0/100, which = T)
top
bottom <- Closest(Light, 10^9, which = T)
bottom
middle <- mean(c(174,780))
middle

light_attenuation <- ggplot(data = Attenuation, aes(x = Light, y = Depth, col = "blue")) +
  geom_hline(yintercept = bottom, col = "darkgreen") + 
  geom_hline(yintercept = top, col = "red") +
  geom_hline(yintercept = middle, col = "blue") +
  geom_point() + geom_line() +
  scale_y_reverse(breaks=seq(0,1000,200)) +
  theme_classic() + theme(legend.position = "none") +
  ylab("Depth (m)") + xlab(expression(Light~intensity~(photons~m^-2~s^-1))) +
  labs(title = "(A)") 

light_attenuation


########################
# Electrona antarctica #
########################

Eant <- species[[1]]
name <- Eant$scientificName[1]
head(Eant)

# Range selection #

Eant$presence <- NA

for(i in 1:nrow(Eant)){
  if(Eant$individualCount[i] > 0){
    Eant$presence[i] <- 1
  } else {Eant$presence[i] <- 0}
}

Eant$presence <- as.factor(Eant$presence)
presence_count <- Eant %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
presence_count 

count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
count_table

count <- vector()

for(i in 1:nrow(count_table)){
  if(count_table$`1`[i] > 0){
    count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
  }
}
count
count <- as.data.frame(count)

f(count$count)

Eant <- Eant[Eant$lat >= f(count$count)[[1]], ]


# Net removal
Eant <- Eant[Eant$netType != "IKMT",]
unique(Eant$netType)


# Raw mean abundance 
Eant_raw <- Eant %>%
  group_by(diel_num, depth, lat) %>%
  summarise(n = n(),
            mean = mean(CPUE))
head(Eant_raw)


# Modelling 
Eant_model <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Eant, select = TRUE,
                   family = "tw", knots = knots, method = "REML")
Eant_summary <- summary(Eant_model)
Eant_summary


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), 
                               pca = 0)
day_predict_data

day_predict <- predict.gam(Eant_model, day_predict_data, 
                           type = "link", se.fit = T)
day_predict$day_predict <- exp(Eant_model$family$linkinv(day_predict$fit)) -1

day_predict_data <- cbind(day_predict, day_predict_data)
head(day_predict_data)

day_abundance <- sum(day_predict_data$day_predict)
day_abundance



# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 3.5, depth = seq(1, 1000, by = 1),
                                 pca = 0)

night_predict <- predict.gam(Eant_model, night_predict_data, 
                           type = "link", se.fit = T)
night_predict$night_predict <- exp(Eant_model$family$linkinv(night_predict$fit)) -1

night_predict_data <- cbind(night_predict, night_predict_data)
head(night_predict_data)

night_abundance <- sum(night_predict_data$night_predict)
night_abundance


# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference

night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop$abun <- round(night_prop$y*10^5, digits = 0)

night_count <- vector()
for(i in 1:nrow(night_prop)){
  night_count <- c(night_count, rep(night_prop$x[i], night_prop$abun[i]))
}
median_night <- median(night_count)
median_night


# Inflation at 174m 
inflation_factor1 <- absolute_difference/top
inflation_factor1

inflation_factor1 <- absolute_difference/sum(day_predict_data[day_predict_data$depth < top,]$day_predict)

day_predict_data$inflated_day1 <- day_predict_data$day_predict
day_predict_data$inflated_day1[0:top] <- day_predict_data$inflated_day1[0:top] * inflation_factor1
day_predict_data$inflated_day1
sum(day_predict_data$inflated_day1)

day_prop <- data.frame(y = day_predict_data$inflated_day1, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day1 <- median(day_count)
median_day1

intercept1 <- (median_night + median_day1)/2
intercept1

proportion1 <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept1,]$inflated_day1) - 
                 sum(night_predict_data[night_predict_data$depth > intercept1,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion1


# Inflation at 477m 
inflation_factor2 <- absolute_difference/middle
inflation_factor2

day_predict_data$inflated_day2 <- day_predict_data$day_predict
day_predict_data$inflated_day2[0:middle] <- day_predict_data$inflated_day2[0:middle] + inflation_factor2
day_predict_data$inflated_day2

day_prop <- data.frame(y = day_predict_data$inflated_day2, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day2 <- median(day_count)
median_day2

intercept2 <- (median_night + median_day2)/2
intercept2

proportion2 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept2,]$inflated_day2) - 
                  sum(night_predict_data[night_predict_data$depth > intercept2,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion2


# Inflation at 780m 
inflation_factor3 <- absolute_difference/bottom
inflation_factor3

day_predict_data$inflated_day3 <- day_predict_data$day_predict
day_predict_data$inflated_day3[0:bottom] <- day_predict_data$inflated_day3[0:bottom] + inflation_factor3
day_predict_data$inflated_day3

day_prop <- data.frame(y = day_predict_data$inflated_day3, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day3 <- median(day_count)
median_day3

intercept3 <- (median_night + median_day3)/2
intercept3

proportion3 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept3,]$inflated_day3) - 
                  sum(night_predict_data[night_predict_data$depth > intercept3,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion3


# Original (equal inflation across all depth)
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
sum(day_predict_data$inflated_day)

day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day4 <- median(day_count)
median_day4

intercept4 <- (median_night + median_day4)/2
intercept4

proportion4 <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept4,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept4,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion4


# Plot
Eant_fit1 <- ggplot() +
  geom_vline(xintercept = intercept1, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept2, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept3, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept4, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day), col = "blue") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day1), col = "red") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day2), col = "pink") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day3, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(A)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
Eant_fit1

Results <- data.frame(Species = name,
                      median_night = median_night,
                      median_day1 = median_day1,
                      median_day2 = median_day2,
                      median_day3 = median_day3,
                      median_day4 = median_day4,
                      Threshold1 = intercept1,
                      Threshold2 = intercept2,
                      Threshold3 = intercept3,
                      Threshold4 = intercept4,
                      proportion1 = round(proportion1*100, digits = 1),
                      proportion2 = round(proportion2*100, digits = 1),
                      proportion3 = round(proportion3*100, digits = 1),
                      proportion4 = round(proportion4*100, digits = 1))
Results
Results_table <- rbind(Results_table, Results)
Results_table



############################
# Krefftichthys anderssoni #
############################

Kand <- species[[3]]
name <- Kand$scientificName[1]
head(Kand)

# Range selection #

Kand$presence <- NA

for(i in 1:nrow(Kand)){
  if(Kand$individualCount[i] > 0){
    Kand$presence[i] <- 1
  } else {Kand$presence[i] <- 0}
}

Kand$presence <- as.factor(Kand$presence)
presence_count <- Kand %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
presence_count 

count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
count_table

count <- vector()

for(i in 1:nrow(count_table)){
  if(count_table$`1`[i] > 0){
    count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
  }
}
count
count <- as.data.frame(count)

f(count$count)

Kand <- Kand[Kand$lat >= f(count$count)[[1]], ]


# Net removal 
Kand <- Kand[Kand$netType != "IKMT",]
unique(Kand$netType)


# Raw mean abundance 
Kand_raw <- Kand %>%
  group_by(diel_num, depth) %>%
  summarise(n = n(),
            mean = mean(CPUE))
head(Kand_raw)


# Modelling 
Kand$netType <- as.factor(Kand$netType)
Kand_model <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Kand, select = TRUE,
                   family = "tw", knots = knots, method = "REML")
Kand_summary <- summary(Kand_model)
Kand_summary


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 2, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Kand_model, day_predict_data, 
                           type = "link", se.fit = T)
day_predict$day_predict <- exp(Kand_model$family$linkinv(day_predict$fit)) -1

day_predict_data <- cbind(day_predict, day_predict_data)
head(day_predict_data)

day_abundance <- sum(day_predict_data$day_predict)
day_abundance


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 4, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Kand_model, night_predict_data, 
                           type = "link", se.fit = T)

night_predict$night_predict <- exp(Kand_model$family$linkinv(night_predict$fit)) -1
night_predict_data <- cbind(night_predict, night_predict_data)
head(night_predict_data)

night_abundance <- sum(night_predict_data$night_predict)
night_abundance

# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference

night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop$abun <- round(night_prop$y*10^5, digits = 0)

night_count <- vector()
for(i in 1:nrow(night_prop)){
  night_count <- c(night_count, rep(night_prop$x[i], night_prop$abun[i]))
}
median_night <- median(night_count)
median_night


# Inflation at 174m 
inflation_factor1 <- absolute_difference/top
inflation_factor1

day_predict_data$inflated_day1 <- day_predict_data$day_predict
day_predict_data$inflated_day1[0:top] <- day_predict_data$inflated_day1[0:top] + inflation_factor1
day_predict_data$inflated_day1

day_prop <- data.frame(y = day_predict_data$inflated_day1, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day1 <- median(day_count)
median_day1

intercept1 <- (median_night + median_day1)/2
intercept1

proportion1 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept1,]$inflated_day1) - 
                  sum(night_predict_data[night_predict_data$depth > intercept1,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion1


# Inflation at 477m 
inflation_factor2 <- absolute_difference/middle
inflation_factor2

day_predict_data$inflated_day2 <- day_predict_data$day_predict
day_predict_data$inflated_day2[0:middle] <- day_predict_data$inflated_day2[0:middle] + inflation_factor2
day_predict_data$inflated_day2

day_prop <- data.frame(y = day_predict_data$inflated_day2, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day2 <- median(day_count)
median_day2

intercept2 <- (median_night + median_day2)/2
intercept2

proportion2 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept2,]$inflated_day2) - 
                  sum(night_predict_data[night_predict_data$depth > intercept2,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion2


# Inflation at 780m 
inflation_factor3 <- absolute_difference/bottom
inflation_factor3

day_predict_data$inflated_day3 <- day_predict_data$day_predict
day_predict_data$inflated_day3[0:bottom] <- day_predict_data$inflated_day3[0:bottom] + inflation_factor3
day_predict_data$inflated_day3

day_prop <- data.frame(y = day_predict_data$inflated_day3, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day3 <- median(day_count)
median_day3

intercept3 <- (median_night + median_day3)/2
intercept3

proportion3 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept3,]$inflated_day3) - 
                  sum(night_predict_data[night_predict_data$depth > intercept3,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion3


# Original (equal inflation across all depth)
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
sum(day_predict_data$inflated_day)

day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day4 <- median(day_count)
median_day4

intercept4 <- (median_night + median_day4)/2
intercept4

proportion4 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept4,]$inflated_day) - 
                  sum(night_predict_data[night_predict_data$depth > intercept4,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion4


# Plot
Kand_fit1 <- ggplot() +
  geom_vline(xintercept = intercept1, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept2, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept3, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept4, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day), col = "blue") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day1), col = "red") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day2), col = "pink") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day3, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(B)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
Kand_fit1

Results <- data.frame(Species = name,
                      median_night = median_night,
                      median_day1 = median_day1,
                      median_day2 = median_day2,
                      median_day3 = median_day3,
                      median_day4 = median_day4,
                      Threshold1 = intercept1,
                      Threshold2 = intercept2,
                      Threshold3 = intercept3,
                      Threshold4 = intercept4,
                      proportion1 = round(proportion1*100, digits = 1),
                      proportion2 = round(proportion2*100, digits = 1),
                      proportion3 = round(proportion3*100, digits = 1),
                      proportion4 = round(proportion4*100, digits = 1))
Results
Results_table <- rbind(Results_table, Results)
Results_table



#########################
# Gymnoscopelus braueri #
#########################

Gbra <- species[[5]]
name <- Gbra$scientificName[1]
head(Gbra)


# Range selection #

Gbra$presence <- NA

for(i in 1:nrow(Gbra)){
  if(Gbra$individualCount[i] > 0){
    Gbra$presence[i] <- 1
  } else {Gbra$presence[i] <- 0}
}

Gbra$presence <- as.factor(Gbra$presence)
presence_count <- Gbra %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
presence_count 

count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
count_table

count <- vector()

for(i in 1:nrow(count_table)){
  if(count_table$`1`[i] > 0){
    count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
  }
}
count
count <- as.data.frame(count)

f(count$count)

Gbra <- Gbra[Gbra$lat >= f(count$count)[[1]], ]


# Net removal 
Gbra <- Gbra[Gbra$netType != "IKMT",]
unique(Gbra$netType)


# Raw mean abundance
Gbra_raw <- Gbra %>%
  group_by(diel_num, depth) %>%
  summarise(n = n(),
            mean = mean(CPUE))
head(Gbra_raw)


# Modelling
Gbra_model <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gbra, select = TRUE,
                   family = "tw", knots = knots, method = "REML")
Gbra_summary <- summary(Gbra_model)
Gbra_summary


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 2, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Gbra_model, day_predict_data, 
                           type = "link", se.fit = T)

day_predict$day_predict <- exp(Gbra_model$family$linkinv(day_predict$fit)) -1
day_predict_data <- cbind(day_predict, day_predict_data)
head(day_predict_data)

day_abundance <- sum(day_predict_data$day_predict)
day_abundance


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 4, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Gbra_model, night_predict_data, 
                           type = "link", se.fit = T)

night_predict$night_predict <- exp(Gbra_model$family$linkinv(night_predict$fit)) -1
night_predict_data <- cbind(night_predict, night_predict_data)
head(night_predict_data)

night_abundance <- sum(night_predict_data$night_predict)
night_abundance


# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference

night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop$abun <- round(night_prop$y*10^5, digits = 0)

night_count <- vector()
for(i in 1:nrow(night_prop)){
  night_count <- c(night_count, rep(night_prop$x[i], night_prop$abun[i]))
}
median_night <- median(night_count)
median_night


# Inflation at 174m 
inflation_factor1 <- absolute_difference/top
inflation_factor1

day_predict_data$inflated_day1 <- day_predict_data$day_predict
day_predict_data$inflated_day1[0:top] <- day_predict_data$inflated_day1[0:top] + inflation_factor1
day_predict_data$inflated_day1

day_prop <- data.frame(y = day_predict_data$inflated_day1, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day1 <- median(day_count)
median_day1

intercept1 <- (median_night + median_day1)/2
intercept1

proportion1 <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept1,]$inflated_day1) - 
                 sum(night_predict_data[night_predict_data$depth > intercept1,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion1


# Inflation at 477m 
inflation_factor2 <- absolute_difference/middle
inflation_factor2

day_predict_data$inflated_day2 <- day_predict_data$day_predict
day_predict_data$inflated_day2[0:middle] <- day_predict_data$inflated_day2[0:middle] + inflation_factor2
day_predict_data$inflated_day2

day_prop <- data.frame(y = day_predict_data$inflated_day2, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day2 <- median(day_count)
median_day2

intercept2 <- (median_night + median_day2)/2
intercept2

proportion2 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept2,]$inflated_day2) - 
                  sum(night_predict_data[night_predict_data$depth > intercept2,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion2


# Inflation at 780m 
inflation_factor3 <- absolute_difference/bottom
inflation_factor3

day_predict_data$inflated_day3 <- day_predict_data$day_predict
day_predict_data$inflated_day3[0:bottom] <- day_predict_data$inflated_day3[0:bottom] + inflation_factor3
day_predict_data$inflated_day3

day_prop <- data.frame(y = day_predict_data$inflated_day3, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day3 <- median(day_count)
median_day3

intercept3 <- (median_night + median_day3)/2
intercept3

proportion3 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept3,]$inflated_day3) - 
                  sum(night_predict_data[night_predict_data$depth > intercept3,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion3


# Original (equal inflation across all depth)
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
sum(day_predict_data$inflated_day)

day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day4 <- median(day_count)
median_day4

intercept4 <- (median_night + median_day4)/2
intercept4

proportion4 <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept4,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept4,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion4


# Plot
Gbra_fit1 <- ggplot() +
  geom_vline(xintercept = intercept1, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept2, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept3, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept4, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day), col = "blue") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day1), col = "red") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day2), col = "pink") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day3, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(A)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
Gbra_fit1

Results <- data.frame(Species = name,
                      median_night = median_night,
                      median_day1 = median_day1,
                      median_day2 = median_day2,
                      median_day3 = median_day3,
                      median_day4 = median_day4,
                      Threshold1 = intercept1,
                      Threshold2 = intercept2,
                      Threshold3 = intercept3,
                      Threshold4 = intercept4,
                      proportion1 = round(proportion1*100, digits = 1),
                      proportion2 = round(proportion2*100, digits = 1),
                      proportion3 = round(proportion3*100, digits = 1),
                      proportion4 = round(proportion4*100, digits = 1))
Results
Results_table <- rbind(Results_table, Results)
Results_table


#########################
# Protomyctophum bolini #
#########################

Pbol <- species[[7]]
name <- Pbol$scientificName[1]
head(Pbol)


# Range selection #

Pbol$presence <- NA

for(i in 1:nrow(Pbol)){
  if(Pbol$individualCount[i] > 0){
    Pbol$presence[i] <- 1
  } else {Pbol$presence[i] <- 0}
}

Pbol$presence <- as.factor(Pbol$presence)
presence_count <- Pbol %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
presence_count 

count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
count_table

count <- vector()

for(i in 1:nrow(count_table)){
  if(count_table$`1`[i] > 0){
    count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
  }
}
count
count <- as.data.frame(count)

f(count$count)

Pbol <- Pbol[Pbol$lat >= f(count$count)[[1]], ]


# Net removal 
Pbol <- Pbol[Pbol$netType != "IKMT",]
unique(Pbol$netType)


# Raw mean abundance 
Pbol_raw <- Pbol %>%
  group_by(diel_num, depth) %>%
  summarise(n = n(),
            mean = mean(CPUE))
head(Pbol_raw)


# Modelling
Pbol_model <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pbol, select = TRUE,
                   family = "tw", knots = knots, method = "REML")
Pbol_summary <- summary(Pbol_model)
Pbol_summary


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Pbol_model, day_predict_data, 
                           type = "link", se.fit = T)
day_predict$day_predict <- exp(Pbol_model$family$linkinv(day_predict$fit)) -1

day_predict_data <- cbind(day_predict, day_predict_data)
head(day_predict_data)

day_abundance <- sum(day_predict_data$day_predict)
day_abundance


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 3.5, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Pbol_model, night_predict_data, 
                           type = "link", se.fit = T)
night_predict$night_predict <- exp(Pbol_model$family$linkinv(night_predict$fit)) -1

night_predict_data <- cbind(night_predict, night_predict_data)
head(night_predict_data)

night_abundance <- sum(night_predict_data$night_predict)
night_abundance


# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference

night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop$abun <- round(night_prop$y*10^5, digits = 0)

night_count <- vector()
for(i in 1:nrow(night_prop)){
  night_count <- c(night_count, rep(night_prop$x[i], night_prop$abun[i]))
}
median_night <- median(night_count)
median_night


# Inflation at 174m 
inflation_factor1 <- absolute_difference/top
inflation_factor1

day_predict_data$inflated_day1 <- day_predict_data$day_predict
day_predict_data$inflated_day1[0:top] <- day_predict_data$inflated_day1[0:top] + inflation_factor1
day_predict_data$inflated_day1

day_prop <- data.frame(y = day_predict_data$inflated_day1, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day1 <- median(day_count)
median_day1

intercept1 <- (median_night + median_day1)/2
intercept1

proportion1 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept1,]$inflated_day1) - 
                  sum(night_predict_data[night_predict_data$depth > intercept1,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion1


# Inflation at 477m 
inflation_factor2 <- absolute_difference/middle
inflation_factor2

day_predict_data$inflated_day2 <- day_predict_data$day_predict
day_predict_data$inflated_day2[0:middle] <- day_predict_data$inflated_day2[0:middle] + inflation_factor2
day_predict_data$inflated_day2

day_prop <- data.frame(y = day_predict_data$inflated_day2, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day2 <- median(day_count)
median_day2

intercept2 <- (median_night + median_day2)/2
intercept2

proportion2 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept2,]$inflated_day2) - 
                  sum(night_predict_data[night_predict_data$depth > intercept2,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion2


# Inflation at 780m 
inflation_factor3 <- absolute_difference/bottom
inflation_factor3

day_predict_data$inflated_day3 <- day_predict_data$day_predict
day_predict_data$inflated_day3[0:bottom] <- day_predict_data$inflated_day3[0:bottom] + inflation_factor3
day_predict_data$inflated_day3

day_prop <- data.frame(y = day_predict_data$inflated_day3, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day3 <- median(day_count)
median_day3

intercept3 <- (median_night + median_day3)/2
intercept3

proportion3 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept3,]$inflated_day3) - 
                  sum(night_predict_data[night_predict_data$depth > intercept3,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion3


# Original (equal inflation across all depth)
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
sum(day_predict_data$inflated_day)

day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day4 <- median(day_count)
median_day4

intercept4 <- (median_night + median_day4)/2
intercept4

proportion4 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept4,]$inflated_day) - 
                  sum(night_predict_data[night_predict_data$depth > intercept4,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion4


# Plot
Pbol_fit1 <- ggplot() +
  geom_vline(xintercept = intercept1, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept2, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept3, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept4, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day), col = "blue") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day1), col = "red") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day2), col = "pink") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day3, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(A)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
Pbol_fit1

Results <- data.frame(Species = name,
                      median_night = median_night,
                      median_day1 = NA,
                      median_day2 = NA,
                      median_day3 = NA,
                      median_day4 = median_day4,
                      Threshold1 = NA,
                      Threshold2 = NA,
                      Threshold3 = NA,
                      Threshold4 = intercept4,
                      proportion1 = NA,
                      proportion2 = NA,
                      proportion3 = NA,
                      proportion4 = round(proportion4*100, digits = 1))
Results
Results_table <- rbind(Results_table, Results)
Results_table




##########################
# Gymnoscopelus nicholsi #
##########################

Gnic <- species[[8]]
name <- Gnic$scientificName[1]
head(Gnic)


# Range selection #

Gnic$presence <- NA

for(i in 1:nrow(Gnic)){
  if(Gnic$individualCount[i] > 0){
    Gnic$presence[i] <- 1
  } else {Gnic$presence[i] <- 0}
}

Gnic$presence <- as.factor(Gnic$presence)
presence_count <- Gnic %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
presence_count 

count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
count_table

count <- vector()

for(i in 1:nrow(count_table)){
  if(count_table$`1`[i] > 0){
    count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
  }
}
count
count <- as.data.frame(count)

f(count$count)

Gnic <- Gnic[Gnic$lat >= f(count$count)[[1]], ]


# Net removal 
Gnic <- Gnic[Gnic$netType != "IKMT",]
unique(Gnic$netType)


# Raw mean abundance 
Gnic_raw <- Gnic %>%
  group_by(diel_num, depth) %>%
  summarise(n = n(),
            mean = mean(CPUE))
head(Gnic_raw)


# Modelling 
Gnic_model <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gnic, select = TRUE,
                   family = "tw", knots = knots, method = "REML")
Gnic_summary <- summary(Gnic_model)
Gnic_summary


# Daytime abundance
day_predict_data <- data.frame(diel_num = 2, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Gnic_model, day_predict_data, 
                           type = "link", se.fit = T)
day_predict$day_predict <- exp(Gnic_model$family$linkinv(day_predict$fit)) -1

day_predict_data <- cbind(day_predict, day_predict_data)
head(day_predict_data)

day_abundance <- sum(day_predict_data$day_predict)
day_abundance


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 4, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Gnic_model, night_predict_data, 
                           type = "link", se.fit = T)
night_predict$night_predict <- exp(Gnic_model$family$linkinv(night_predict$fit)) -1

night_predict_data <- cbind(night_predict, night_predict_data)
head(night_predict_data)

night_abundance <- sum(night_predict_data$night_predict)
night_abundance


# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference

night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop$abun <- round(night_prop$y*10^5, digits = 0)

night_count <- vector()
for(i in 1:nrow(night_prop)){
  night_count <- c(night_count, rep(night_prop$x[i], night_prop$abun[i]))
}
median_night <- median(night_count)
median_night


# Inflation at 174m 
inflation_factor1 <- absolute_difference/top
inflation_factor1

day_predict_data$inflated_day1 <- day_predict_data$day_predict
day_predict_data$inflated_day1[0:top] <- day_predict_data$inflated_day1[0:top] + inflation_factor1
day_predict_data$inflated_day1

day_prop <- data.frame(y = day_predict_data$inflated_day1, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day1 <- median(day_count)
median_day1

intercept1 <- (median_night + median_day1)/2
intercept1

proportion1 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept1,]$inflated_day1) - 
                  sum(night_predict_data[night_predict_data$depth > intercept1,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion1


# Inflation at 477m 
inflation_factor2 <- absolute_difference/middle
inflation_factor2

day_predict_data$inflated_day2 <- day_predict_data$day_predict
day_predict_data$inflated_day2[0:middle] <- day_predict_data$inflated_day2[0:middle] + inflation_factor2
day_predict_data$inflated_day2

day_prop <- data.frame(y = day_predict_data$inflated_day2, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day2 <- median(day_count)
median_day2

intercept2 <- (median_night + median_day2)/2
intercept2

proportion2 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept2,]$inflated_day2) - 
                  sum(night_predict_data[night_predict_data$depth > intercept2,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion2


# Inflation at 780m 
inflation_factor3 <- absolute_difference/bottom
inflation_factor3

day_predict_data$inflated_day3 <- day_predict_data$day_predict
day_predict_data$inflated_day3[0:bottom] <- day_predict_data$inflated_day3[0:bottom] + inflation_factor3
day_predict_data$inflated_day3

day_prop <- data.frame(y = day_predict_data$inflated_day3, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day3 <- median(day_count)
median_day3

intercept3 <- (median_night + median_day3)/2
intercept3

proportion3 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept3,]$inflated_day3) - 
                  sum(night_predict_data[night_predict_data$depth > intercept3,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion3


# Original (equal inflation across all depth)
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
sum(day_predict_data$inflated_day)

day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day4 <- median(day_count)
median_day4

intercept4 <- (median_night + median_day4)/2
intercept4

proportion4 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept4,]$inflated_day) - 
                  sum(night_predict_data[night_predict_data$depth > intercept4,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion4


# Plot
Gnic_fit1 <- ggplot() +
  geom_vline(xintercept = intercept1, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept2, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept3, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept4, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day), col = "blue") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day1), col = "red") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day2), col = "pink") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day3, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(A)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
Gnic_fit1

Results <- data.frame(Species = name,
                      median_night = median_night,
                      median_day1 = median_day1,
                      median_day2 = median_day2,
                      median_day3 = median_day3,
                      median_day4 = median_day4,
                      Threshold1 = intercept1,
                      Threshold2 = intercept2,
                      Threshold3 = intercept3,
                      Threshold4 = intercept4,
                      proportion1 = round(proportion1*100, digits = 1),
                      proportion2 = round(proportion2*100, digits = 1),
                      proportion3 = round(proportion3*100, digits = 1),
                      proportion4 = round(proportion4*100, digits = 1))
Results
Results_table <- rbind(Results_table, Results)
Results_table



##########################
# Gymnoscopelus fraseri #
##########################

Gfra <- species[[10]]
name <- Gfra$scientificName[1]
head(Gfra)


# Range selection #

Gfra$presence <- NA

for(i in 1:nrow(Gfra)){
  if(Gfra$individualCount[i] > 0){
    Gfra$presence[i] <- 1
  } else {Gfra$presence[i] <- 0}
}

Gfra$presence <- as.factor(Gfra$presence)
presence_count <- Gfra %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
presence_count 

count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
count_table

count <- vector()

for(i in 1:nrow(count_table)){
  if(count_table$`1`[i] > 0){
    count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
  }
}
count
count <- as.data.frame(count)

f(count$count)

Gfra <- Gfra[Gfra$lat >= f(count$count)[[1]], ]


# Net removal 
Gfra <- Gfra[Gfra$netType != "IKMT",]
unique(Gfra$netType)


# Raw mean abundance 
Gfra_raw <- Gfra %>%
  group_by(diel_num, depth) %>%
  summarise(n = n(),
            mean = mean(CPUE))
head(Gfra_raw)


# Modelling 
Gfra_model <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gfra, select = TRUE,
                   family = "tw", knots = knots, method = "REML")
Gfra_summary <- summary(Gfra_model)
Gfra_summary


# Daytime abundance
day_predict_data <- data.frame(diel_num = 2, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Gfra_model, day_predict_data, 
                           type = "link", se.fit = T)
day_predict$day_predict <- exp(Gfra_model$family$linkinv(day_predict$fit)) -1

day_predict_data <- cbind(day_predict, day_predict_data)
head(day_predict_data)

day_abundance <- sum(day_predict_data$day_predict)
day_abundance


# Nighttime abundance
night_predict_data <- data.frame(diel_num = 4, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Gfra_model, night_predict_data, 
                           type = "link", se.fit = T)
night_predict$night_predict <- exp(Gfra_model$family$linkinv(night_predict$fit)) -1

night_predict_data <- cbind(night_predict, night_predict_data)
head(night_predict_data)

night_abundance <- sum(night_predict_data$night_predict)
night_abundance


# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference

night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop$abun <- round(night_prop$y*10^5, digits = 0)

night_count <- vector()
for(i in 1:nrow(night_prop)){
  night_count <- c(night_count, rep(night_prop$x[i], night_prop$abun[i]))
}
median_night <- median(night_count)
median_night


# Inflation at 174m 
inflation_factor1 <- absolute_difference/top
inflation_factor1

day_predict_data$inflated_day1 <- day_predict_data$day_predict
day_predict_data$inflated_day1[0:top] <- day_predict_data$inflated_day1[0:top] + inflation_factor1
day_predict_data$inflated_day1

day_prop <- data.frame(y = day_predict_data$inflated_day1, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day1 <- median(day_count)
median_day1

intercept1 <- (median_night + median_day1)/2
intercept1

proportion1 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept1,]$inflated_day1) - 
                  sum(night_predict_data[night_predict_data$depth > intercept1,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion1


# Inflation at 477m 
inflation_factor2 <- absolute_difference/middle
inflation_factor2

day_predict_data$inflated_day2 <- day_predict_data$day_predict
day_predict_data$inflated_day2[0:middle] <- day_predict_data$inflated_day2[0:middle] + inflation_factor2
day_predict_data$inflated_day2

day_prop <- data.frame(y = day_predict_data$inflated_day2, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day2 <- median(day_count)
median_day2

intercept2 <- (median_night + median_day2)/2
intercept2

proportion2 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept2,]$inflated_day2) - 
                  sum(night_predict_data[night_predict_data$depth > intercept2,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion2


# Inflation at 780m 
inflation_factor3 <- absolute_difference/bottom
inflation_factor3

day_predict_data$inflated_day3 <- day_predict_data$day_predict
day_predict_data$inflated_day3[0:bottom] <- day_predict_data$inflated_day3[0:bottom] + inflation_factor3
day_predict_data$inflated_day3

day_prop <- data.frame(y = day_predict_data$inflated_day3, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day3 <- median(day_count)
median_day3

intercept3 <- (median_night + median_day3)/2
intercept3

proportion3 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept3,]$inflated_day3) - 
                  sum(night_predict_data[night_predict_data$depth > intercept3,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion3


# Original (equal inflation across all depth)
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
sum(day_predict_data$inflated_day)

day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day4 <- median(day_count)
median_day4

intercept4 <- (median_night + median_day4)/2
intercept4

proportion4 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept4,]$inflated_day) - 
                  sum(night_predict_data[night_predict_data$depth > intercept4,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion4


# Plot
Gfra_fit1 <- ggplot() +
  geom_vline(xintercept = intercept1, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept2, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept3, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept4, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day), col = "blue") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day1), col = "red") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day2), col = "pink") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day3, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(A)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
Gfra_fit1

Results <- data.frame(Species = name,
                      median_night = median_night,
                      median_day1 = median_day1,
                      median_day2 = median_day2,
                      median_day3 = median_day3,
                      median_day4 = median_day4,
                      Threshold1 = intercept1,
                      Threshold2 = intercept2,
                      Threshold3 = intercept3,
                      Threshold4 = intercept4,
                      proportion1 = round(proportion1*100, digits = 1),
                      proportion2 = round(proportion2*100, digits = 1),
                      proportion3 = round(proportion3*100, digits = 1),
                      proportion4 = round(proportion4*100, digits = 1))
Results
Results_table <- rbind(Results_table, Results)
Results_table



###########################
# Protomyctophum tenisoni #
###########################

Pten <- species[[11]]
name <- Pten$scientificName[1]
head(Pten)


# Range selection #

Pten$presence <- NA

for(i in 1:nrow(Pten)){
  if(Pten$individualCount[i] > 0){
    Pten$presence[i] <- 1
  } else {Pten$presence[i] <- 0}
}

Pten$presence <- as.factor(Pten$presence)
presence_count <- Pten %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
presence_count 

count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
count_table

count <- vector()

for(i in 1:nrow(count_table)){
  if(count_table$`1`[i] > 0){
    count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
  }
}
count
count <- as.data.frame(count)

f(count$count)

Pten <- Pten[Pten$lat >= f(count$count)[[1]], ]


# Net removal 
Pten <- Pten[Pten$netType != "IKMT",]
unique(Pten$netType)


# Mean raw abundance 
Pten_raw <- Pten %>%
  group_by(diel_num, depth) %>%
  summarise(n = n(),
            mean = mean(CPUE))
head(Pten_raw)


# Modelling
Pten_model <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pten, select = TRUE,
                   family = "tw", knots = knots, method = "REML")
Pten_summary <- summary(Pten_model)
Pten_summary


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 2, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Pten_model, day_predict_data, 
                           type = "link", se.fit = T)
day_predict$day_predict <- exp(Pten_model$family$linkinv(day_predict$fit)) -1

day_predict_data <- cbind(day_predict, day_predict_data)
head(day_predict_data)

day_abundance <- sum(day_predict_data$day_predict)
day_abundance


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 4, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Pten_model, night_predict_data, 
                           type = "link", se.fit = T)
night_predict$night_predict <- exp(Pten_model$family$linkinv(night_predict$fit)) -1

night_predict_data <- cbind(night_predict, night_predict_data)
head(night_predict_data)

night_abundance <- sum(night_predict_data$night_predict)
night_abundance

# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference

night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop$abun <- round(night_prop$y*10^5, digits = 0)

night_count <- vector()
for(i in 1:nrow(night_prop)){
  night_count <- c(night_count, rep(night_prop$x[i], night_prop$abun[i]))
}
median_night <- median(night_count)
median_night


# Inflation at 174m 
inflation_factor1 <- absolute_difference/top
inflation_factor1

day_predict_data$inflated_day1 <- day_predict_data$day_predict
day_predict_data$inflated_day1[0:top] <- day_predict_data$inflated_day1[0:top] + inflation_factor1
day_predict_data$inflated_day1

day_prop <- data.frame(y = day_predict_data$inflated_day1, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day1 <- median(day_count)
median_day1

intercept1 <- (median_night + median_day1)/2
intercept1

proportion1 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept1,]$inflated_day1) - 
                  sum(night_predict_data[night_predict_data$depth > intercept1,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion1


# Inflation at 477m 
inflation_factor2 <- absolute_difference/middle
inflation_factor2

day_predict_data$inflated_day2 <- day_predict_data$day_predict
day_predict_data$inflated_day2[0:middle] <- day_predict_data$inflated_day2[0:middle] + inflation_factor2
day_predict_data$inflated_day2

day_prop <- data.frame(y = day_predict_data$inflated_day2, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day2 <- median(day_count)
median_day2

intercept2 <- (median_night + median_day2)/2
intercept2

proportion2 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept2,]$inflated_day2) - 
                  sum(night_predict_data[night_predict_data$depth > intercept2,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion2


# Inflation at 780m 
inflation_factor3 <- absolute_difference/bottom
inflation_factor3

day_predict_data$inflated_day3 <- day_predict_data$day_predict
day_predict_data$inflated_day3[0:bottom] <- day_predict_data$inflated_day3[0:bottom] + inflation_factor3
day_predict_data$inflated_day3

day_prop <- data.frame(y = day_predict_data$inflated_day3, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day3 <- median(day_count)
median_day3

intercept3 <- (median_night + median_day3)/2
intercept3

proportion3 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept3,]$inflated_day3) - 
                  sum(night_predict_data[night_predict_data$depth > intercept3,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion3


# Original (equal inflation across all depth)
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
sum(day_predict_data$inflated_day)

day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day4 <- median(day_count)
median_day4

intercept4 <- (median_night + median_day4)/2
intercept4

proportion4 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept4,]$inflated_day) - 
                  sum(night_predict_data[night_predict_data$depth > intercept4,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion4


# Plot
Pten_fit1 <- ggplot() +
  geom_vline(xintercept = intercept1, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept2, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept3, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept4, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day), col = "blue") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day1), col = "red") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day2), col = "pink") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day3, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(A)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
Pten_fit1

Results <- data.frame(Species = name,
                      median_night = median_night,
                      median_day1 = median_day1,
                      median_day2 = median_day2,
                      median_day3 = median_day3,
                      median_day4 = median_day4,
                      Threshold1 = intercept1,
                      Threshold2 = intercept2,
                      Threshold3 = intercept3,
                      Threshold4 = intercept4,
                      proportion1 = round(proportion1*100, digits = 1),
                      proportion2 = round(proportion2*100, digits = 1),
                      proportion3 = round(proportion3*100, digits = 1),
                      proportion4 = round(proportion4*100, digits = 1))
Results
Results_table <- rbind(Results_table, Results)
Results_table



########################
# Electrona carlsbergi #
########################

Ecar <- species[[12]]
name <- Ecar$scientificName[1]
head(Ecar)


# Range selection #

Ecar$presence <- NA

for(i in 1:nrow(Ecar)){
  if(Ecar$individualCount[i] > 0){
    Ecar$presence[i] <- 1
  } else {Ecar$presence[i] <- 0}
}

Ecar$presence <- as.factor(Ecar$presence)
presence_count <- Ecar %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
presence_count 

count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
count_table

count <- vector()

for(i in 1:nrow(count_table)){
  if(count_table$`1`[i] > 0){
    count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
  }
}
count
count <- as.data.frame(count)

f(count$count)

Ecar <- Ecar[Ecar$lat >= f(count$count)[[1]], ]


# Net removal 
Ecar <- Ecar[Ecar$netType != "IKMT",]
unique(Ecar$netType)


# Raw mean abundance 
Ecar_raw <- Ecar %>%
  group_by(diel_num, depth) %>%
  summarise(n = n(),
            mean = mean(CPUE))
head(Ecar_raw)


# Modelling 
Ecar_model <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Ecar, select = TRUE,
                   family = "tw", knots = knots, method = "REML")
Ecar_summary <- summary(Ecar_model)
Ecar_summary


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Ecar_model, day_predict_data, 
                           type = "link", se.fit = T)
day_predict$day_predict <- exp(Ecar_model$family$linkinv(day_predict$fit)) -1

day_predict_data <- cbind(day_predict, day_predict_data)
head(day_predict_data)

day_abundance <- sum(day_predict_data$day_predict)
day_abundance


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 3.5, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Ecar_model, night_predict_data, 
                           type = "link", se.fit = T)
night_predict$night_predict <- exp(Ecar_model$family$linkinv(night_predict$fit)) -1

night_predict_data <- cbind(night_predict, night_predict_data)
head(night_predict_data)

night_abundance <- sum(night_predict_data$night_predict)
night_abundance

# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference

night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop$abun <- round(night_prop$y*10^5, digits = 0)

night_count <- vector()
for(i in 1:nrow(night_prop)){
  night_count <- c(night_count, rep(night_prop$x[i], night_prop$abun[i]))
}
median_night <- median(night_count)
median_night


# Inflation at 174m 
inflation_factor1 <- absolute_difference/top
inflation_factor1

day_predict_data$inflated_day1 <- day_predict_data$day_predict
day_predict_data$inflated_day1[0:top] <- day_predict_data$inflated_day1[0:top] + inflation_factor1
day_predict_data$inflated_day1

day_prop <- data.frame(y = day_predict_data$inflated_day1, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day1 <- median(day_count)
median_day1

intercept1 <- (median_night + median_day1)/2
intercept1

proportion1 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept1,]$inflated_day1) - 
                  sum(night_predict_data[night_predict_data$depth > intercept1,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion1


# Inflation at 477m 
inflation_factor2 <- absolute_difference/middle
inflation_factor2

day_predict_data$inflated_day2 <- day_predict_data$day_predict
day_predict_data$inflated_day2[0:middle] <- day_predict_data$inflated_day2[0:middle] + inflation_factor2
day_predict_data$inflated_day2

day_prop <- data.frame(y = day_predict_data$inflated_day2, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day2 <- median(day_count)
median_day2

intercept2 <- (median_night + median_day2)/2
intercept2

proportion2 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept2,]$inflated_day2) - 
                  sum(night_predict_data[night_predict_data$depth > intercept2,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion2


# Inflation at 780m 
inflation_factor3 <- absolute_difference/bottom
inflation_factor3

day_predict_data$inflated_day3 <- day_predict_data$day_predict
day_predict_data$inflated_day3[0:bottom] <- day_predict_data$inflated_day3[0:bottom] + inflation_factor3
day_predict_data$inflated_day3

day_prop <- data.frame(y = day_predict_data$inflated_day3, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day3 <- median(day_count)
median_day3

intercept3 <- (median_night + median_day3)/2
intercept3

proportion3 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept3,]$inflated_day3) - 
                  sum(night_predict_data[night_predict_data$depth > intercept3,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion3


# Original (equal inflation across all depth)
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
sum(day_predict_data$inflated_day)

day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day4 <- median(day_count)
median_day4

intercept4 <- (median_night + median_day4)/2
intercept4

proportion4 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept4,]$inflated_day) - 
                  sum(night_predict_data[night_predict_data$depth > intercept4,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion4


# Plot
Ecar_fit1 <- ggplot() +
  geom_vline(xintercept = intercept1, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept2, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept3, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_vline(xintercept = intercept4, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day), col = "blue") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day1), col = "red") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day2), col = "pink") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day3, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(A)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
Ecar_fit1

Results <- data.frame(Species = name,
                      median_night = median_night,
                      median_day1 = median_day1,
                      median_day2 = median_day2,
                      median_day3 = median_day3,
                      median_day4 = median_day4,
                      Threshold1 = intercept1,
                      Threshold2 = intercept2,
                      Threshold3 = intercept3,
                      Threshold4 = intercept4,
                      proportion1 = round(proportion1*100, digits = 1),
                      proportion2 = round(proportion2*100, digits = 1),
                      proportion3 = round(proportion3*100, digits = 1),
                      proportion4 = round(proportion4*100, digits = 1))
Results
Results_table <- rbind(Results_table, Results)
Results_table



###############
# All species #
###############

# Choosing study species #
species <- group[which(group$scientificName %in% c("Electrona antarctica", "Electrona carlsbergi", 
                                                   "Gymnoscopelus braueri", 
                                                   "Gymnoscopelus fraseri", "Gymnoscopelus nicholsi", 
                                                   "Krefftichthys anderssoni", "Protomyctophum bolini",
                                                   "Protomyctophum tenisoni")),]
unique(species$scientificName)

species_sum <- species %>%
  group_by(eventID) %>%
  summarise(sum = sum(individualCount))
species_sum

new <- merge(data, species_sum, by = "eventID", all.x = TRUE, all.y = FALSE)

new["sum"][is.na(new["sum"])] <- 0

new$CPUE <- (new$sum/new$volume)*1000
new$logCPUE <- log(new$CPUE + 1)


# Modelling #

name <- "All species"

# Mean raw abundance 
raw <- new %>% group_by(lat, diel_num, depth) %>% 
  summarise(mean = mean(CPUE), sum = sum(sum), count = n())


# Remove nets
new <- new[new$lat >= -65,]
new_IKMT <- new[new$netType != "IKMT",]


# Modelling 
All_model1 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                    s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                  select = TRUE, data = new_IKMT,
                  family = "tw", knots = knots, method = "REML")

All_model_summary <- summary(All_model1)
All_model_summary


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(All_model1, day_predict_data, 
                           type = "link", se.fit = T)
day_predict$CI <- 1.96 * day_predict$se.fit

day_predict$day_predict <- exp(All_model1$family$linkinv(day_predict$fit)) -1
day_predict$upr <- exp(All_model1$family$linkinv(day_predict$fit + day_predict$CI)) -1
day_predict$lwr <- exp(All_model1$family$linkinv(day_predict$fit - day_predict$CI)) -1

day_predict_data <- cbind(day_predict, day_predict_data)
head(day_predict_data)

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_upr <- sum(day_predict_data$upr)
day_abundance_upr
day_abundance_lwr <- sum(day_predict_data$lwr)
day_abundance_lwr


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 3.5, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(All_model1, night_predict_data, 
                           type = "link", se.fit = T)
night_predict$CI <- 1.96 * night_predict$se.fit

night_predict$night_predict <- exp(All_model1$family$linkinv(night_predict$fit)) -1
night_predict$upr <- exp(All_model1$family$linkinv(night_predict$fit + night_predict$CI)) -1
night_predict$lwr <- exp(All_model1$family$linkinv(night_predict$fit - night_predict$CI)) -1

night_predict_data <- cbind(night_predict, night_predict_data)
head(night_predict_data)


night_abundance <- sum(night_predict_data$night_predict)
night_abundance
night_abundance_upr <- sum(night_predict_data$upr)
night_abundance_upr
night_abundance_lwr <- sum(night_predict_data$lwr)
night_abundance_lwr


# Day-night difference
absolute_difference <- night_abundance-day_abundance
absolute_difference

night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop$abun <- round(night_prop$y*10^5, digits = 0)

night_count <- vector()
for(i in 1:nrow(night_prop)){
  night_count <- c(night_count, rep(night_prop$x[i], night_prop$abun[i]))
}
median_night <- median(night_count)
median_night


# Inflation at 174m 
inflation_factor1 <- absolute_difference/top
inflation_factor1

day_predict_data$inflated_day1 <- day_predict_data$day_predict
day_predict_data$inflated_day1[0:top] <- day_predict_data$inflated_day1[0:top] + inflation_factor1
day_predict_data$inflated_day1

day_prop <- data.frame(y = day_predict_data$inflated_day1, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day1 <- median(day_count)
median_day1

intercept1 <- (median_night + median_day1)/2
intercept1

proportion1 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept1,]$inflated_day1) - 
                  sum(night_predict_data[night_predict_data$depth > intercept1,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion1


# Inflation at 477m 
inflation_factor2 <- absolute_difference/middle
inflation_factor2

day_predict_data$inflated_day2 <- day_predict_data$day_predict
day_predict_data$inflated_day2[0:middle] <- day_predict_data$inflated_day2[0:middle] + inflation_factor2
day_predict_data$inflated_day2

day_prop <- data.frame(y = day_predict_data$inflated_day2, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day2 <- median(day_count)
median_day2

intercept2 <- (median_night + median_day2)/2
intercept2

proportion2 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept2,]$inflated_day2) - 
                  sum(night_predict_data[night_predict_data$depth > intercept2,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion2


# Inflation at 780m 
inflation_factor3 <- absolute_difference/bottom
inflation_factor3

day_predict_data$inflated_day3 <- day_predict_data$day_predict
day_predict_data$inflated_day3[0:bottom] <- day_predict_data$inflated_day3[0:bottom] + inflation_factor3
day_predict_data$inflated_day3

day_prop <- data.frame(y = day_predict_data$inflated_day3, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day3 <- median(day_count)
median_day3

intercept3 <- (median_night + median_day3)/2
intercept3

proportion3 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept3,]$inflated_day3) - 
                  sum(night_predict_data[night_predict_data$depth > intercept3,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion3


# Original (equal inflation across all depth)
difference_factor <- night_abundance/day_abundance 
difference_factor

day_predict_data$inflated_day <- day_predict_data$day_predict * difference_factor
sum(day_predict_data$inflated_day)

day_prop <- data.frame(y = day_predict_data$inflated_day, x = day_predict_data$depth)

day_prop$abun <- round(day_prop$y*10^5, digits = 0)
head(day_prop)

day_count <- vector()
for(i in 1:nrow(day_prop)){
  day_count <- c(day_count, rep(day_prop$x[i], day_prop$abun[i]))
}
median_day4 <- median(day_count)
median_day4

intercept4 <- (median_night + median_day4)/2
intercept4

proportion4 <- (sum(night_predict_data$night_predict) - 
                  sum(day_predict_data[day_predict_data$depth <= intercept4,]$inflated_day) - 
                  sum(night_predict_data[night_predict_data$depth > intercept4,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion4


# Plot
All_fit1 <- ggplot() +
  geom_vline(xintercept = intercept1, alpha = 0.3, linewidth = 1, color = "red") +
  geom_vline(xintercept = intercept2, alpha = 0.3, linewidth = 1, color = "blue") +
  geom_vline(xintercept = intercept3, alpha = 0.3, linewidth = 1, color = "darkgreen") +
  geom_vline(xintercept = intercept4, alpha = 0.3, linewidth = 1, color = "orange") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.15) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.15) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 20),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 20),], 
             aes(x = depth, y = inflated_day1), col = "red") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 20),], 
             aes(x = depth, y = inflated_day2), col = "blue") +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 20),], 
             aes(x = depth, y = inflated_day3), col = "darkgreen") +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = "(B)") +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
All_fit1

Results <- data.frame(Species = name,
                      median_night = median_night,
                      median_day1 = median_day1,
                      median_day2 = median_day2,
                      median_day3 = median_day3,
                      median_day4 = median_day4,
                      Threshold1 = intercept1,
                      Threshold2 = intercept2,
                      Threshold3 = intercept3,
                      Threshold4 = intercept4,
                      proportion1 = round(proportion1*100, digits = 1),
                      proportion2 = round(proportion2*100, digits = 1),
                      proportion3 = round(proportion3*100, digits = 1),
                      proportion4 = round(proportion4*100, digits = 1))
Results
Results_table <- rbind(Results_table, Results)
Results_table




####################
# Exporting result # 
####################

# Model prediction
Results_table
Final_Results_table <- Results_table[2:10,]
Final_Results_table

mean(Results_table$proportion1/Results_table$proportion4, na.rm = T)
mean(Results_table$proportion2/Results_table$proportion4, na.rm = T)
mean(Results_table$proportion3/Results_table$proportion4, na.rm = T)

Results_table

write.csv(Results_table, "Light_inflation.csv", row.names = F)


# Light inflation
light_attenuation + All_fit1
# 1200W x 600H
