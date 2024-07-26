####################
### DVM analysis ###
####################

rm(list = ls())

require(dplyr)
require(geomtextpath)
require(ggplot2)
require(mgcv)
require(patchwork)



####################
## Data wrangling ##
####################

# Loading data #
data <- as.data.frame(read.csv("Myctobase/event_edited2.csv", header = TRUE, stringsAsFactors = F))
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
                      Peak_day = NA, Peak_night = NA,
                      centre_day = NA, centre_night = NA,
                      day_abundance = NA, night_abundance = NA,
                      Threshold = NA,
                      proportion = NA)
Results_table

GAM_result <- data.frame(Species = NA, Interaction = NA, depth = NA, time = NA, net = NA,
                         R_squared = NA, N = NA, lat = NA)
GAM_result



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
                   family = "tw", knots = knots, method = "ML")
Eant_summary <- summary(Eant_model)
Eant_summary
plot(Eant_model, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)


result <- data.frame(Species = name, 
                         Interaction = round(Eant_summary$s.table[1,4], digits = 5), 
                         depth = round(Eant_summary$s.table[2,4], digits = 5),
                         time = round(Eant_summary$s.table[3,4], digits = 5),
                         net = round(Eant_summary$p.pv[[2]], digits =5),
                         R_squared = round(Eant_summary$dev.expl, digits = 3), 
                         N = Eant_summary$n,
                         lat = abs(f(count$count)[[1]]))
result
GAM_result <- rbind(GAM_result, result)
GAM_result


# Plotting
Eant_predict_data <- data.frame(diel_num = rep(seq(0.5, 4.5, length = 100), each = 100), 
                                depth = rep(seq(0, 1000, length = 100), 100),
                                volume = 1000, pca = 0)
Eant_predict1 <- predict.gam(Eant_model, Eant_predict_data, 
                             type = "response")
Eant_predict_data1 <- cbind(Eant_predict1, Eant_predict_data)
head(Eant_predict_data1)
Eant_predict_data1$Eant_predict <- exp(Eant_predict_data1$Eant_predict1) -1

Eant_plot1 <- ggplot() +
  geom_tile(data = Eant_predict_data1, 
            aes(x = diel_num, y = depth, fill = Eant_predict)) +
  scale_fill_gradientn(colors = hcl.colors(20, "viridis")) +
  geom_point(data = Eant_raw[Eant_raw$mean > 0 & Eant_raw$depth < 1000,], 
             aes(x = diel_num, y = depth, size = mean, alpha = 0.5,col = "white"),
             shape = 1) +
  geom_textcontour(data = Eant_predict_data1, 
                   aes(x = diel_num, y = depth, z = Eant_predict,
                       linecolor = "white", textcolor = "white",
                       label = after_stat(round(level, digits =  2))),
                   bins = 8, straight = TRUE) +
  geom_rug(data = Eant[Eant$depth <1000,], 
           inherit.aes = FALSE,
           aes(y=depth, col = "white")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,10,5),
        plot.title = element_text(size = 12, margin = margin(0,0,0.5,0))) +
  scale_y_reverse(expand = c(0, 0), breaks=seq(0,1000,200)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks=c(1,2,3,4),
                     labels = c("Dawn", "Day", "Dusk", "Night")) +
  ggtitle(expression(italic('(A) Electrona antarctica')))
Eant_plot1


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), 
                               pca = 0)
day_predict_data

day_predict <- predict.gam(Eant_model, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_predict_data$upr <- day_predict_data$day_predict + (2 * day_predict_data$SE)
day_predict_data$lwr <- day_predict_data$day_predict - (2 * day_predict_data$SE)
day_predict_data

peak_day <- day_predict_data$depth[which.max(day_predict_data$day_predict)]
peak_day

central_day <- sum(day_predict_data$day_predict * day_predict_data$depth)/
  sum(day_predict_data$day_predict)
central_day

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_SE <- sum(day_predict_data$SE)
day_abundance_SE


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 3.5, depth = seq(1, 1000, by = 1),
                                 pca = 0)
night_predict_data

night_predict <- predict.gam(Eant_model, night_predict_data, 
                             type = "response", se.fit = T)
night_predict_data <- cbind(night_predict, night_predict_data)
night_predict_data
night_predict_data$night_predict <- exp(night_predict_data$fit) -1
night_predict_data$SE <- exp(night_predict_data$se.fit) -1

peak_night <- night_predict_data$depth[which.max(night_predict_data$night_predict)]
peak_night
  
central_night <- sum(night_predict_data$night_predict * night_predict_data$depth)/
  sum(night_predict_data$night_predict)
central_night

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

intercept <- day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 1)]
intercept


# Plot
day_predict_data$lwr[day_predict_data$lwr < 0] <- 0
night_predict_data$lwr[night_predict_data$lwr < 0] <- 0

Eant_fit <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Eant[Eant$depth < 1000 & 
                           Eant$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Eant[Eant$depth < 1000 & 
                           Eant$diel_num == 4,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(A)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() 
Eant_fit

Eant_fit1 <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Eant[Eant$depth < 1000 & 
                           Eant$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Eant[Eant$depth < 1000 & 
                           Eant$diel_num == 4,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
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


# Proportion of population doing DVM
proportion <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion


Results <- data.frame(Species = name,
                      Peak_day = peak_day, Peak_night = peak_night,
                      centre_day = round(central_day, digits = 0), 
                      centre_night = round(central_night, digits = 0),
                      day_abundance = paste0(round(day_abundance/1000, digits = 3), 
                                             " (", round(day_abundance_SE/1000, digits = 3), ")"), 
                      night_abundance = paste0(round(night_abundance/1000, digits = 3), 
                                               " (", round(night_abundance_SE/1000, digits = 3), ")"),
                      Threshold = intercept,
                      proportion = round(proportion*100, digits = 1))
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
                   family = "tw", knots = knots, method = "ML")
Kand_summary <- summary(Kand_model)
Kand_summary
plot(Kand_model, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)

result <- data.frame(Species = name, 
                     Interaction = round(Kand_summary$s.table[1,4], digits = 5), 
                     depth = round(Kand_summary$s.table[2,4], digits = 5),
                     time = round(Kand_summary$s.table[3,4], digits = 5),
                     net = round(Kand_summary$p.pv[[2]], digits =5),
                     R_squared = round(Kand_summary$dev.expl, digits = 3), 
                     N = Kand_summary$n,
                     lat = abs(f(count$count)[[1]]))
result
GAM_result <- rbind(GAM_result, result)
GAM_result


# Plotting
Kand_predict_data <- data.frame(diel_num = rep(seq(0.5, 4.5, length = 100), each = 100), 
                                depth = rep(seq(0, 1000, length = 100), 100),
                                volume = 1000, pca = 0)
Kand_predict1 <- predict.gam(Kand_model, Kand_predict_data, 
                             type = "response", exclude = "s(netType)", newdata.guaranteed=TRUE)
Kand_predict_data1 <- cbind(Kand_predict1, Kand_predict_data)
head(Kand_predict_data1)
Kand_predict_data1$Kand_predict <- exp(Kand_predict_data1$Kand_predict1) -1

Kand_plot1 <- ggplot() +
  geom_tile(data = Kand_predict_data1, 
            aes(x = diel_num, y = depth, fill = Kand_predict)) +
  scale_fill_gradientn(colors = hcl.colors(20, "viridis")) +
  geom_point(data = Kand_raw[Kand_raw$mean > 0 & Kand_raw$depth < 1000,], 
             aes(x = diel_num, y = depth, size = mean, alpha = 0.5,col = "white"),
             shape = 1) +
  geom_textcontour(data = Kand_predict_data1, 
                   aes(x = diel_num, y = depth, z = Kand_predict,
                       linecolor = "white", textcolor = "white",
                       label = after_stat(round(level, digits =  2))),
                   bins = 8, straight = TRUE) +
  geom_rug(data = Kand[Kand$depth <1000,], 
           inherit.aes = FALSE,
           aes(y=depth, col = "white")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,10,5),
        plot.title = element_text(size = 12, margin = margin(0,0,0.5,0))) +
  scale_y_reverse(expand = c(0, 0), breaks=seq(0,1000,200)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks=c(1,2,3,4),
                     labels = c("Dawn", "Day", "Dusk", "Night")) +
  ggtitle(expression(italic('(B) Krefftichthys anderssoni')))
Kand_plot1


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 2, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Kand_model, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_predict_data$upr <- day_predict_data$day_predict + (2 * day_predict_data$SE)
day_predict_data$lwr <- day_predict_data$day_predict - (2 * day_predict_data$SE)
day_predict_data

peak_day <- day_predict_data$depth[which.max(day_predict_data$day_predict)]
peak_day

central_day <- sum(day_predict_data$day_predict * day_predict_data$depth)/
  sum(day_predict_data$day_predict)
central_day

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_SE <- sum(day_predict_data$SE)
day_abundance_SE


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 4, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Kand_model, night_predict_data, 
                             type = "response", se.fit = T)
night_predict_data <- cbind(night_predict, night_predict_data)
night_predict_data
night_predict_data$night_predict <- exp(night_predict_data$fit) -1
night_predict_data$SE <- exp(night_predict_data$se.fit) -1

night_predict_data$upr <- night_predict_data$night_predict + (2 * night_predict_data$SE)
night_predict_data$lwr <- night_predict_data$night_predict - (2 * night_predict_data$SE)
night_predict_data

peak_night <- night_predict_data$depth[which.max(night_predict_data$night_predict)]
peak_night

central_night <- sum(night_predict_data$night_predict * night_predict_data$depth)/
  sum(night_predict_data$night_predict)
central_night

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

day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 2)]
intercept <- day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 1)]
intercept


# Plot
day_predict_data$lwr[day_predict_data$lwr < 0] <- 0
night_predict_data$lwr[night_predict_data$lwr < 0] <- 0

Kand_fit <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Kand[Kand$depth < 1000 & 
                           Kand$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Kand[Kand$depth < 1000 & 
                           Kand$diel_num == 4 ,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(B)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = c(0.95, 1.1), legend.justification = c("right", "top")) +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() 
Kand_fit

Kand_fit1 <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Kand[Kand$depth < 1000 & 
                           Kand$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Kand[Kand$depth < 1000 & 
                           Kand$diel_num == 4,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(B)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = c(0.95, 1.1), legend.justification = c("right", "top")) +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() +  coord_flip() 
Kand_fit1


# Proportion of population doing DVM
proportion <- (sum(night_predict_data$night_predict) - 
  sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
  sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion


Results <- data.frame(Species = name,
                      Peak_day = peak_day, Peak_night = peak_night,
                      centre_day = round(central_day, digits = 0), 
                      centre_night = round(central_night, digits = 0),
                      day_abundance = paste0(round(day_abundance/1000, digits = 3), 
                                             " (", round(day_abundance_SE/1000, digits = 3), ")"), 
                      night_abundance = paste0(round(night_abundance/1000, digits = 3), 
                                               " (", round(night_abundance_SE/1000, digits = 3), ")"),
                      Threshold = intercept,
                      proportion = round(proportion*100, digits = 1))
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
                   family = "tw", knots = knots, method = "ML")
Gbra_summary <- summary(Gbra_model)
Gbra_summary
plot(Gbra_model, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)


result <- data.frame(Species = name, 
                     Interaction = round(Gbra_summary$s.table[1,4], digits = 5), 
                     depth = round(Gbra_summary$s.table[2,4], digits = 5),
                     time = round(Gbra_summary$s.table[3,4], digits = 5),
                     net = round(Gbra_summary$p.pv[[2]], digits =5),
                     R_squared = round(Gbra_summary$dev.expl, digits = 3), 
                     N = Gbra_summary$n,
                     lat = abs(f(count$count)[[1]]))
result
GAM_result <- rbind(GAM_result, result)
GAM_result


# Plotting
Gbra_predict_data <- data.frame(diel_num = rep(seq(0.5, 4.5, length = 100), each = 100), 
                                depth = rep(seq(0, 1000, length = 100), 100),
                                volume = 1000, pca = 0)
Gbra_predict1 <- predict.gam(Gbra_model, Gbra_predict_data, 
                             type = "response")
Gbra_predict_data1 <- cbind(Gbra_predict1, Gbra_predict_data)
head(Gbra_predict_data1)
Gbra_predict_data1$Gbra_predict <- exp(Gbra_predict_data1$Gbra_predict1) -1

Gbra_plot1 <- ggplot() +
  geom_tile(data = Gbra_predict_data1, 
            aes(x = diel_num, y = depth, fill = Gbra_predict)) +
  scale_fill_gradientn(colors = hcl.colors(20, "viridis")) +
  geom_point(data = Gbra_raw[Gbra_raw$mean > 0 & Gbra_raw$depth < 1000,], 
             aes(x = diel_num, y = depth, size = mean, alpha = 0.5,col = "white"),
             shape = 1) +
  geom_textcontour(data = Gbra_predict_data1, 
                   aes(x = diel_num, y = depth, z = Gbra_predict,
                       linecolor = "white", textcolor = "white",
                       label = after_stat(round(level, digits =  2))),
                   bins = 8, straight = TRUE) +
  geom_rug(data = Gbra[Gbra$depth <1000,], 
           inherit.aes = FALSE,
           aes(y=depth, col = "white")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,10,5),
        plot.title = element_text(size = 12, margin = margin(0,0,0.5,0))) +
  scale_y_reverse(expand = c(0, 0), breaks=seq(0,1000,200)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks=c(1,2,3,4),
                     labels = c("Dawn", "Day", "Dusk", "Night")) +
  ggtitle(expression(italic('(C) Gymnoscopelus braueri')))
Gbra_plot1


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 2, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Gbra_model, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_predict_data$upr <- day_predict_data$day_predict + (2 * day_predict_data$SE)
day_predict_data$lwr <- day_predict_data$day_predict - (2 * day_predict_data$SE)
day_predict_data

peak_day <- day_predict_data$depth[which.max(day_predict_data$day_predict)]
peak_day

central_day <- sum(day_predict_data$day_predict * day_predict_data$depth)/
  sum(day_predict_data$day_predict)
central_day

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_SE <- sum(day_predict_data$SE)
day_abundance_SE


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 4, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Gbra_model, night_predict_data, 
                             type = "response", se.fit = T)
night_predict_data <- cbind(night_predict, night_predict_data)
night_predict_data
night_predict_data$night_predict <- exp(night_predict_data$fit) -1
night_predict_data$SE <- exp(night_predict_data$se.fit) -1

night_predict_data$upr <- night_predict_data$night_predict + (2 * night_predict_data$SE)
night_predict_data$lwr <- night_predict_data$night_predict - (2 * night_predict_data$SE)
night_predict_data

peak_night <- night_predict_data$depth[which.max(night_predict_data$night_predict)]
peak_night

central_night <- sum(night_predict_data$night_predict * night_predict_data$depth)/
  sum(night_predict_data$night_predict)
central_night

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

intercept <- day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 1)]
intercept


# Plot
day_predict_data$lwr[day_predict_data$lwr < 0] <- 0
night_predict_data$lwr[night_predict_data$lwr < 0] <- 0

Gbra_fit <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Gbra[Gbra$depth < 1000 & 
                           Gbra$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Gbra[Gbra$depth < 1000 & 
                           Gbra$diel_num == 4 ,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(C)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() 
Gbra_fit

Gbra_fit1 <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Gbra[Gbra$depth < 1000 & 
                           Gbra$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Gbra[Gbra$depth < 1000 & 
                           Gbra$diel_num == 4,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(C)"~italic(.(name)))) +
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


# Proportion of population doing DVM
proportion <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion


Results <- data.frame(Species = name,
                      Peak_day = peak_day, Peak_night = peak_night,
                      centre_day = round(central_day, digits = 0), 
                      centre_night = round(central_night, digits = 0),
                      day_abundance = paste0(round(day_abundance/1000, digits = 3), 
                                             " (", round(day_abundance_SE/1000, digits = 3), ")"), 
                      night_abundance = paste0(round(night_abundance/1000, digits = 3), 
                                               " (", round(night_abundance_SE/1000, digits = 3), ")"),
                      Threshold = intercept,
                      proportion = round(proportion*100, digits = 1))
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
                   family = "tw", knots = knots, method = "ML")
Pbol_summary <- summary(Pbol_model)
Pbol_summary
plot(Pbol_model, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)


result <- data.frame(Species = name, 
                     Interaction = round(Pbol_summary$s.table[1,4], digits = 5), 
                     depth = round(Pbol_summary$s.table[2,4], digits = 5),
                     time = round(Pbol_summary$s.table[3,4], digits = 5),
                     net = round(Pbol_summary$p.pv[[2]], digits =5),
                     R_squared = round(Pbol_summary$dev.expl, digits = 3), 
                     N = Pbol_summary$n,
                     lat = abs(f(count$count)[[1]]))
result
GAM_result <- rbind(GAM_result, result)
GAM_result


# Plotting
Pbol_predict_data <- data.frame(diel_num = rep(seq(0.5, 4.5, length = 100), each = 100), 
                                depth = rep(seq(0, 1000, length = 100), 100),
                                volume = 1000, pca = 0)
Pbol_predict1 <- predict.gam(Pbol_model, Pbol_predict_data, 
                             type = "response")
Pbol_predict_data1 <- cbind(Pbol_predict1, Pbol_predict_data)
head(Pbol_predict_data1)
Pbol_predict_data1$Pbol_predict <- exp(Pbol_predict_data1$Pbol_predict1) -1

Pbol_plot1 <- ggplot() +
  geom_tile(data = Pbol_predict_data1, 
            aes(x = diel_num, y = depth, fill = Pbol_predict)) +
  scale_fill_gradientn(colors = hcl.colors(20, "viridis")) +
  geom_point(data = Pbol_raw[Pbol_raw$mean > 0 & Pbol_raw$depth < 1000,], 
             aes(x = diel_num, y = depth, size = mean, alpha = 0.5,col = "white"),
             shape = 1) +
  geom_textcontour(data = Pbol_predict_data1, 
                   aes(x = diel_num, y = depth, z = Pbol_predict,
                       linecolor = "white", textcolor = "white",
                       label = after_stat(round(level, digits =  2))),
                   bins = 8, straight = TRUE) +
  geom_rug(data = Pbol[Pbol$depth <1000,], 
           inherit.aes = FALSE,
           aes(y=depth, col = "white")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,10,5),
        plot.title = element_text(size = 12, margin = margin(0,0,0.5,0))) +
  scale_y_reverse(expand = c(0, 0), breaks=seq(0,1000,200)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks=c(1,2,3,4),
                     labels = c("Dawn", "Day", "Dusk", "Night")) +
  ggtitle(expression(italic('(D) Protomyctophum bolini')))
Pbol_plot1


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Pbol_model, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_predict_data$upr <- day_predict_data$day_predict + (2 * day_predict_data$SE)
day_predict_data$lwr <- day_predict_data$day_predict - (2 * day_predict_data$SE)
day_predict_data

peak_day <- day_predict_data$depth[which.max(day_predict_data$day_predict)]
peak_day

central_day <- sum(day_predict_data$day_predict * day_predict_data$depth)/
  sum(day_predict_data$day_predict)
central_day

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_SE <- sum(day_predict_data$SE)
day_abundance_SE


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 3.5, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Pbol_model, night_predict_data, 
                             type = "response", se.fit = T)
night_predict_data <- cbind(night_predict, night_predict_data)
night_predict_data
night_predict_data$night_predict <- exp(night_predict_data$fit) -1
night_predict_data$SE <- exp(night_predict_data$se.fit) -1

night_predict_data$upr <- night_predict_data$night_predict + (2 * night_predict_data$SE)
night_predict_data$lwr <- night_predict_data$night_predict - (2 * night_predict_data$SE)
night_predict_data

peak_night <- night_predict_data$depth[which.max(night_predict_data$night_predict)]
peak_night

central_night <- sum(night_predict_data$night_predict * night_predict_data$depth)/
  sum(night_predict_data$night_predict)
central_night

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

day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 3)]
intercept <- day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 3)][1]
intercept


# Plot
day_predict_data$lwr[day_predict_data$lwr < 0] <- 0
night_predict_data$lwr[night_predict_data$lwr < 0] <- 0

Pbol_fit <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Pbol[Pbol$depth < 1000 & 
                           Pbol$diel_num == 2,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Pbol[Pbol$depth < 1000 & 
                           Pbol$diel_num == 4 ,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 30),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(D)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() 
Pbol_fit


Pbol_fit1 <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Pbol[Pbol$depth < 1000 & 
                           Pbol$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Pbol[Pbol$depth < 1000 & 
                           Pbol$diel_num == 4,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(D)"~italic(.(name)))) +
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


# Proportion of population doing DVM
proportion <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion


Results <- data.frame(Species = name,
                      Peak_day = peak_day, Peak_night = peak_night,
                      centre_day = round(central_day, digits = 0), 
                      centre_night = round(central_night, digits = 0),
                      day_abundance = paste0(round(day_abundance/1000, digits = 3), 
                                             " (", round(day_abundance_SE/1000, digits = 3), ")"), 
                      night_abundance = paste0(round(night_abundance/1000, digits = 3), 
                                               " (", round(night_abundance_SE/1000, digits = 3), ")"),
                      Threshold = intercept,
                      proportion = round(proportion*100, digits = 1))
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
                   family = "tw", knots = knots, method = "ML")
Gnic_summary <- summary(Gnic_model)
Gnic_summary
plot(Gnic_model, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)

result <- data.frame(Species = name, 
                     Interaction = round(Gnic_summary$s.table[1,4], digits = 5), 
                     depth = round(Gnic_summary$s.table[2,4], digits = 5),
                     time = round(Gnic_summary$s.table[3,4], digits = 5),
                     net = round(Gnic_summary$p.pv[[2]], digits =5),
                     R_squared = round(Gnic_summary$dev.expl, digits = 3), 
                     N = Gnic_summary$n,
                     lat = abs(f(count$count)[[1]]))
result
GAM_result <- rbind(GAM_result, result)
GAM_result


# Plotting
Gnic_predict_data <- data.frame(diel_num = rep(seq(0.5, 4.5, length = 100), each = 100), 
                                depth = rep(seq(0, 1000, length = 100), 100),
                                volume = 1000, pca = 0)
Gnic_predict1 <- predict.gam(Gnic_model, Gnic_predict_data, 
                             type = "response")
Gnic_predict_data1 <- cbind(Gnic_predict1, Gnic_predict_data)
head(Gnic_predict_data1)
Gnic_predict_data1$Gnic_predict <- exp(Gnic_predict_data1$Gnic_predict1) -1

Gnic_plot1 <- ggplot() +
  geom_tile(data = Gnic_predict_data1, 
            aes(x = diel_num, y = depth, fill = Gnic_predict)) +
  scale_fill_gradientn(colors = hcl.colors(20, "viridis")) +
  geom_point(data = Gnic_raw[Gnic_raw$mean > 0 & Gnic_raw$depth < 1000,], 
             aes(x = diel_num, y = depth, size = mean, alpha = 0.5,col = "white"),
             shape = 1) +
  geom_textcontour(data = Gnic_predict_data1, 
                   aes(x = diel_num, y = depth, z = Gnic_predict,
                       linecolor = "white", textcolor = "white",
                       label = after_stat(round(level, digits =  2))),
                   bins = 10, straight = TRUE) +
  geom_rug(data = Gnic[Gnic$depth <1000,], 
           inherit.aes = FALSE,
           aes(y=depth, col = "white")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,10,5),
        plot.title = element_text(size = 12, margin = margin(0,0,0.5,0))) +
  scale_y_reverse(expand = c(0, 0), breaks=seq(0,1000,200)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks=c(1,2,3,4),
                     labels = c("Dawn", "Day", "Dusk", "Night")) +
  ggtitle(expression(italic('(E) Gymnoscopelus nicholsi')))
Gnic_plot1


# Daytime abundance
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Gnic_model, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_predict_data$upr <- day_predict_data$day_predict + (2 * day_predict_data$SE)
day_predict_data$lwr <- day_predict_data$day_predict - (2 * day_predict_data$SE)
day_predict_data

peak_day <- day_predict_data$depth[which.max(day_predict_data$day_predict)]
peak_day

central_day <- sum(day_predict_data$day_predict * day_predict_data$depth)/
  sum(day_predict_data$day_predict)
central_day

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_SE <- sum(day_predict_data$SE)
day_abundance_SE


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 3.5, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Gnic_model, night_predict_data, 
                             type = "response", se.fit = T)
night_predict_data <- cbind(night_predict, night_predict_data)
night_predict_data
night_predict_data$night_predict <- exp(night_predict_data$fit) -1
night_predict_data$SE <- exp(night_predict_data$se.fit) -1

night_predict_data$upr <- night_predict_data$night_predict + (2 * night_predict_data$SE)
night_predict_data$lwr <- night_predict_data$night_predict - (2 * night_predict_data$SE)
night_predict_data

peak_night <- night_predict_data$depth[which.max(night_predict_data$night_predict)]
peak_night

central_night <- sum(night_predict_data$night_predict * night_predict_data$depth)/
  sum(night_predict_data$night_predict)
central_night

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
day_prop <- day_prop[day_prop$x < 400,]
night_prop <- data.frame(y = night_predict_data$night_predict, x = night_predict_data$depth)
night_prop <- night_prop[night_prop$x < 400,]

day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 1)]
intercept <- day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 1)]
intercept


# Plot
day_predict_data$lwr[day_predict_data$lwr < 0] <- 0
night_predict_data$lwr[night_predict_data$lwr < 0] <- 0

Gnic_fit <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Gnic[Gnic$depth < 1000 & 
                           Gnic$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Gnic[Gnic$depth < 1000 & 
                           Gnic$diel_num == 4 ,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(E)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() 
Gnic_fit


Gnic_fit1 <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Gnic[Gnic$depth < 1000 & 
                           Gnic$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Gnic[Gnic$depth < 1000 & 
                           Gnic$diel_num == 4,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(E)"~italic(.(name)))) +
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


# Proportion of population doing DVM
proportion <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion


Results <- data.frame(Species = name,
                      Peak_day = peak_day, Peak_night = peak_night,
                      centre_day = round(central_day, digits = 0), 
                      centre_night = round(central_night, digits = 0),
                      day_abundance = paste0(round(day_abundance/1000, digits = 3), 
                                             " (", round(day_abundance_SE/1000, digits = 3), ")"), 
                      night_abundance = paste0(round(night_abundance/1000, digits = 3), 
                                               " (", round(night_abundance_SE/1000, digits = 3), ")"),
                      Threshold = intercept,
                      proportion = round(proportion*100, digits = 1))
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
                   family = "tw", knots = knots, method = "ML")
Gfra_summary <- summary(Gfra_model)
Gfra_summary
plot(Gfra_model, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)

result <- data.frame(Species = name, 
                     Interaction = round(Gfra_summary$s.table[1,4], digits = 5), 
                     depth = round(Gfra_summary$s.table[2,4], digits = 5),
                     time = round(Gfra_summary$s.table[3,4], digits = 5),
                     net = round(Gfra_summary$p.pv[[2]], digits =5),
                     R_squared = round(Gfra_summary$dev.expl, digits = 3), 
                     N = Gfra_summary$n,
                     lat = abs(f(count$count)[[1]]))
result
GAM_result <- rbind(GAM_result, result)
GAM_result


# Plotting
Gfra_predict_data <- data.frame(diel_num = rep(seq(0.5, 4.5, length = 100), each = 100), 
                                depth = rep(seq(0, 1000, length = 100), 100),
                                volume = 1000, pca = 0)
Gfra_predict1 <- predict.gam(Gfra_model, Gfra_predict_data, 
                             type = "response")
Gfra_predict_data1 <- cbind(Gfra_predict1, Gfra_predict_data)
head(Gfra_predict_data1)
Gfra_predict_data1$Gfra_predict <- exp(Gfra_predict_data1$Gfra_predict1) -1

Gfra_plot1 <- ggplot() +
  geom_tile(data = Gfra_predict_data1, 
            aes(x = diel_num, y = depth, fill = Gfra_predict)) +
  scale_fill_gradientn(colors = hcl.colors(20, "viridis")) +
  geom_point(data = Gfra_raw[Gfra_raw$mean > 0 & Gfra_raw$depth < 1000,], 
             aes(x = diel_num, y = depth, size = mean, alpha = 0.5,col = "white"),
             shape = 1) +
  geom_textcontour(data = Gfra_predict_data1, 
                   aes(x = diel_num, y = depth, z = Gfra_predict,
                       linecolor = "white", textcolor = "white",
                       label = after_stat(round(level, digits =  2))),
                   bins = 10, straight = TRUE) +
  geom_rug(data = Gfra[Gfra$depth <1000,], 
           inherit.aes = FALSE,
           aes(y=depth, col = "white")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,10,5),
        plot.title = element_text(size = 12, margin = margin(0,0,0.5,0))) +
  scale_y_reverse(expand = c(0, 0), breaks=seq(0,1000,200)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks=c(1,2,3,4),
                     labels = c("Dawn", "Day", "Dusk", "Night")) +
  ggtitle(expression(italic('(F) Gymnoscopelus fraseri')))
Gfra_plot1


# Daytime abundance
day_predict_data <- data.frame(diel_num = 2, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Gfra_model, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_predict_data$upr <- day_predict_data$day_predict + (2 * day_predict_data$SE)
day_predict_data$lwr <- day_predict_data$day_predict - (2 * day_predict_data$SE)
day_predict_data

peak_day <- day_predict_data$depth[which.max(day_predict_data$day_predict)]
peak_day

central_day <- sum(day_predict_data$day_predict * day_predict_data$depth)/
  sum(day_predict_data$day_predict)
central_day

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_SE <- sum(day_predict_data$SE)
day_abundance_SE


# Nighttime abundance
night_predict_data <- data.frame(diel_num = 4, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Gfra_model, night_predict_data, 
                             type = "response", se.fit = T)
night_predict_data <- cbind(night_predict, night_predict_data)
night_predict_data
night_predict_data$night_predict <- exp(night_predict_data$fit) -1
night_predict_data$SE <- exp(night_predict_data$se.fit) -1

night_predict_data$upr <- night_predict_data$night_predict + (2 * night_predict_data$SE)
night_predict_data$lwr <- night_predict_data$night_predict - (2 * night_predict_data$SE)
night_predict_data

peak_night <- night_predict_data$depth[which.max(night_predict_data$night_predict)]
peak_night

central_night <- sum(night_predict_data$night_predict * night_predict_data$depth)/
  sum(night_predict_data$night_predict)
central_night

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

day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 2)]
intercept <- day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 2)][2]
intercept


# Plot
day_predict_data$lwr[day_predict_data$lwr < 0] <- 0
night_predict_data$lwr[night_predict_data$lwr < 0] <- 0

Gfra_fit <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Gfra[Gfra$depth < 1000 & 
                           Gfra$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Gfra[Gfra$depth < 1000 & 
                           Gfra$diel_num == 4 ,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 30),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(F)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() 
Gfra_fit

Gfra_fit1 <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Gfra[Gfra$depth < 1000 & 
                           Gfra$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Gfra[Gfra$depth < 1000 & 
                           Gfra$diel_num == 4,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(F)"~italic(.(name)))) +
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


# Proportion of population doing DVM
proportion <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion


Results <- data.frame(Species = name,
                      Peak_day = peak_day, Peak_night = peak_night,
                      centre_day = round(central_day, digits = 0), 
                      centre_night = round(central_night, digits = 0),
                      day_abundance = paste0(round(day_abundance/1000, digits = 3), 
                                             " (", round(day_abundance_SE/1000, digits = 3), ")"), 
                      night_abundance = paste0(round(night_abundance/1000, digits = 3), 
                                               " (", round(night_abundance_SE/1000, digits = 3), ")"),
                      Threshold = intercept,
                      proportion = round(proportion*100, digits = 1))
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
                   family = "tw", knots = knots, method = "ML")
Pten_summary <- summary(Pten_model)
Pten_summary
plot(Pten_model, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)

result <- data.frame(Species = name, 
                     Interaction = round(Pten_summary$s.table[1,4], digits = 5), 
                     depth = round(Pten_summary$s.table[2,4], digits = 5),
                     time = round(Pten_summary$s.table[3,4], digits = 5),
                     net = round(Pten_summary$p.pv[[2]], digits =5),
                     R_squared = round(Pten_summary$dev.expl, digits = 3), 
                     N = Pten_summary$n,
                     lat = abs(f(count$count)[[1]]))
result
GAM_result <- rbind(GAM_result, result)
GAM_result


# Plotting
Pten_predict_data <- data.frame(diel_num = rep(seq(0.5, 4.5, length = 100), each = 100), 
                                depth = rep(seq(0, 1000, length = 100), 100),
                                volume = 1000, pca = 0)
Pten_predict1 <- predict.gam(Pten_model, Pten_predict_data, 
                             type = "response")
Pten_predict_data1 <- cbind(Pten_predict1, Pten_predict_data)
head(Pten_predict_data1)
Pten_predict_data1$Pten_predict <- exp(Pten_predict_data1$Pten_predict1) -1

Pten_plot1 <- ggplot() +
  geom_tile(data = Pten_predict_data1, 
            aes(x = diel_num, y = depth, fill = Pten_predict)) +
  scale_fill_gradientn(colors = hcl.colors(20, "viridis")) +
  geom_point(data = Pten_raw[Pten_raw$mean > 0 & Pten_raw$depth < 1000,], 
             aes(x = diel_num, y = depth, size = mean, alpha = 0.5,col = "white"),
             shape = 1) +
  geom_textcontour(data = Pten_predict_data1, 
                   aes(x = diel_num, y = depth, z = Pten_predict,
                       linecolor = "white", textcolor = "white",
                       label = after_stat(round(level, digits =  2))),
                   bins = 10, straight = TRUE) +
  geom_rug(data = Pten[Pten$depth <1000,], 
           inherit.aes = FALSE,
           aes(y=depth, col = "white")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,10,5),
        plot.title = element_text(size = 12, margin = margin(0,0,0.5,0))) +
  scale_y_reverse(expand = c(0, 0), breaks=seq(0,1000,200)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks=c(1,2,3,4),
                     labels = c("Dawn", "Day", "Dusk", "Night")) +
  ggtitle(expression(italic('(G) Protomyctophum tenisoni')))
Pten_plot1


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 2, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Pten_model, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_predict_data$upr <- day_predict_data$day_predict + (2 * day_predict_data$SE)
day_predict_data$lwr <- day_predict_data$day_predict - (2 * day_predict_data$SE)
day_predict_data

peak_day <- day_predict_data$depth[which.max(day_predict_data$day_predict)]
peak_day

central_day <- sum(day_predict_data$day_predict * day_predict_data$depth)/
  sum(day_predict_data$day_predict)
central_day

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_SE <- sum(day_predict_data$SE)
day_abundance_SE


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 4, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Pten_model, night_predict_data, 
                             type = "response", se.fit = T)
night_predict_data <- cbind(night_predict, night_predict_data)
night_predict_data
night_predict_data$night_predict <- exp(night_predict_data$fit) -1
night_predict_data$SE <- exp(night_predict_data$se.fit) -1

night_predict_data$upr <- night_predict_data$night_predict + (2 * night_predict_data$SE)
night_predict_data$lwr <- night_predict_data$night_predict - (2 * night_predict_data$SE)
night_predict_data

peak_night <- night_predict_data$depth[which.max(night_predict_data$night_predict)]
peak_night

central_night <- sum(night_predict_data$night_predict * night_predict_data$depth)/
  sum(night_predict_data$night_predict)
central_night

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

day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 10)]
intercept <- day_prop$x[which.mins(abs(day_prop$y - night_prop$y), 10)][3]
intercept


# Plot
day_predict_data$lwr[day_predict_data$lwr < 0] <- 0
night_predict_data$lwr[night_predict_data$lwr < 0] <- 0

Pten_fit <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Pten[Pten$depth < 1000 & 
                           Pten$diel_num == 2,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Pten[Pten$depth < 1000 & 
                           Pten$diel_num == 4 ,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(G)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() 
Pten_fit

Pten_fit1 <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Pten[Pten$depth < 1000 & 
                           Pten$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Pten[Pten$depth < 1000 & 
                           Pten$diel_num == 4,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(G)"~italic(.(name)))) +
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


# Proportion of population doing DVM
proportion <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion


Results <- data.frame(Species = name,
                      Peak_day = peak_day, Peak_night = peak_night,
                      centre_day = round(central_day, digits = 0), 
                      centre_night = round(central_night, digits = 0),
                      day_abundance = paste0(round(day_abundance/1000, digits = 3), 
                                             " (", round(day_abundance_SE/1000, digits = 3), ")"), 
                      night_abundance = paste0(round(night_abundance/1000, digits = 3), 
                                               " (", round(night_abundance_SE/1000, digits = 3), ")"),
                      Threshold = intercept,
                      proportion = round(proportion*100, digits = 1))
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
                   family = "tw", knots = knots, method = "ML")
Ecar_summary <- summary(Ecar_model)
Ecar_summary
plot(Ecar_model, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)


result <- data.frame(Species = name, 
                     Interaction = round(Ecar_summary$s.table[1,4], digits = 5), 
                     depth = round(Ecar_summary$s.table[2,4], digits = 5),
                     time = round(Ecar_summary$s.table[3,4], digits = 5),
                     net = round(Ecar_summary$p.pv[[2]], digits =5),
                     R_squared = round(Ecar_summary$dev.expl, digits = 3), 
                     N = Ecar_summary$n,
                     lat = abs(f(count$count)[[1]]))
result
GAM_result <- rbind(GAM_result, result)
GAM_result


# Plotting
Ecar_predict_data <- data.frame(diel_num = rep(seq(0.5, 4.5, length = 100), each = 100), 
                                depth = rep(seq(0, 1000, length = 100), 100),
                                volume = 1000, pca = 0)
Ecar_predict1 <- predict.gam(Ecar_model, Ecar_predict_data, 
                             type = "response")
Ecar_predict_data1 <- cbind(Ecar_predict1, Ecar_predict_data)
head(Ecar_predict_data1)
Ecar_predict_data1$Ecar_predict <- exp(Ecar_predict_data1$Ecar_predict1) -1

Ecar_plot1 <- ggplot() +
  geom_tile(data = Ecar_predict_data1, 
            aes(x = diel_num, y = depth, fill = Ecar_predict)) +
  scale_fill_gradientn(colors = hcl.colors(20, "viridis")) +
  geom_point(data = Ecar_raw[Ecar_raw$mean > 0 & Ecar_raw$depth < 1000,], 
             aes(x = diel_num, y = depth, size = mean, alpha = 0.5,col = "white"),
             shape = 1) +
  geom_textcontour(data = Ecar_predict_data1, 
                   aes(x = diel_num, y = depth, z = Ecar_predict,
                       linecolor = "white", textcolor = "white",
                       label = after_stat(round(level, digits =  2))),
                   bins = 8, straight = TRUE) +
  geom_rug(data = Ecar[Ecar$depth <1000,], 
           inherit.aes = FALSE,
           aes(y=depth, col = "white")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,10,5),
        plot.title = element_text(size = 12, margin = margin(0,0,0.5,0))) +
  scale_y_reverse(expand = c(0, 0), breaks=seq(0,1000,200)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks=c(1,2,3,4),
                     labels = c("Dawn", "Day", "Dusk", "Night")) +
  ggtitle(expression(italic('(H) Electrona carlsbergi')))
Ecar_plot1


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(Ecar_model, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_predict_data$upr <- day_predict_data$day_predict + (2 * day_predict_data$SE)
day_predict_data$lwr <- day_predict_data$day_predict - (2 * day_predict_data$SE)
day_predict_data

peak_day <- day_predict_data$depth[which.max(day_predict_data$day_predict)]
peak_day

central_day <- sum(day_predict_data$day_predict * day_predict_data$depth)/
  sum(day_predict_data$day_predict)
central_day

day_abundance <- sum(day_predict_data$day_predict)
day_abundance
day_abundance_SE <- sum(day_predict_data$SE)
day_abundance_SE


# Nighttime abundance 
night_predict_data <- data.frame(diel_num = 3.5, depth = seq(1, 1000, by = 1), pca = 0)
night_predict_data

night_predict <- predict.gam(Ecar_model, night_predict_data, 
                             type = "response", se.fit = T)
night_predict_data <- cbind(night_predict, night_predict_data)
night_predict_data
night_predict_data$night_predict <- exp(night_predict_data$fit) -1
night_predict_data$SE <- exp(night_predict_data$se.fit) -1

night_predict_data$upr <- night_predict_data$night_predict + (2 * night_predict_data$SE)
night_predict_data$lwr <- night_predict_data$night_predict - (2 * night_predict_data$SE)
night_predict_data

peak_night <- night_predict_data$depth[which.max(night_predict_data$night_predict)]
peak_night

central_night <- sum(night_predict_data$night_predict * night_predict_data$depth)/
  sum(night_predict_data$night_predict)
central_night

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

day_prop$x[which.mins(abs(day_prop[day_prop$x < 400,]$y - night_prop[night_prop$x < 400,]$y), 10)]
intercept <- day_prop$x[which.mins(abs(day_prop[day_prop$x < 400,]$y - night_prop[night_prop$x < 400,]$y), 10)][7]
intercept


# Plot
day_predict_data$lwr[day_predict_data$lwr < 0] <- 0
night_predict_data$lwr[night_predict_data$lwr < 0] <- 0

Ecar_fit <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Ecar[Ecar$depth < 1000 & 
                           Ecar$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Ecar[Ecar$depth < 1000 & 
                           Ecar$diel_num == 4 ,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 15),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(H)"~italic(.(name)))) +
  labs(y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = "none") +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() 
Ecar_fit

Ecar_fit1 <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = Ecar[Ecar$depth < 1000 & 
                           Ecar$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = Ecar[Ecar$depth < 1000 & 
                           Ecar$diel_num == 4,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = bquote("(H)"~italic(.(name)))) +
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


# Proportion of population doing DVM
proportion <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion


Results <- data.frame(Species = name,
                      Peak_day = peak_day, Peak_night = peak_night,
                      centre_day = round(central_day, digits = 0), 
                      centre_night = round(central_night, digits = 0),
                      day_abundance = paste0(round(day_abundance/1000, digits = 3), 
                                             " (", round(day_abundance_SE/1000, digits = 3), ")"), 
                      night_abundance = paste0(round(night_abundance/1000, digits = 3), 
                                               " (", round(night_abundance_SE/1000, digits = 3), ")"),
                      Threshold = intercept,
                      proportion = round(proportion*100, digits = 1))
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
                  family = "tw", knots = knots, method = "ML")

All_model_summary <- summary(All_model1)
All_model_summary
plot(All_model1, scheme = 2, scale = 0, pages = 1, all.terms = TRUE)

result <- data.frame(Species = name, 
                     Interaction = round(All_model_summary$s.table[1,4], digits = 5), 
                     depth = round(All_model_summary$s.table[2,4], digits = 5),
                     time = round(All_model_summary$s.table[3,4], digits = 5),
                     net = round(All_model_summary$p.pv[[2]], digits =5),
                     R_squared = round(All_model_summary$dev.expl, digits = 3), 
                     N = All_model_summary$n,
                     lat = 65)
result
GAM_result <- rbind(GAM_result, result)
GAM_result


# Plotting
All_predict_data <- data.frame(diel_num = rep(seq(0.5, 4.5, length = 100), each = 100), 
                               depth = rep(seq(0, 1000, length = 100), 100),
                               volume = 1000, pca = 0)

All_predict <- predict.gam(All_model1, All_predict_data, 
                           type = "response")
All_predict_data <- cbind(All_predict, All_predict_data)
All_predict_data$All_predict <- exp(All_predict_data$All_predict) -1

All_plot1 <- ggplot() +
  geom_tile(data = All_predict_data, 
            aes(x = diel_num, y = depth, fill = All_predict)) +
  scale_fill_gradientn(colors = hcl.colors(20, "viridis")) +
  geom_point(data = raw, 
             aes(x = diel_num, y = depth, size = mean, alpha = 0.5,col = "white"),
             shape = 1) +
  geom_textcontour(data = All_predict_data, 
                   aes(x = diel_num, y = depth, z = All_predict,
                       linecolor = "white", textcolor = "white",
                       label = after_stat(round(level, digits =  2))),
                   bins = 8, straight = TRUE) +
  geom_rug(data = new_IKMT, 
           inherit.aes = FALSE,
           aes(y=depth, col = "white")) +
  theme_classic() + labs(y = "Depth (m)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.margin = margin(1,20,10,5)) +
  scale_y_reverse(expand = c(0, 0), breaks=seq(0,1000,200)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks=c(1,2,3,4),
                     labels = c("Dawn", "Day", "Dusk", "Night")) +
  labs(title ="(A)") 
All_plot1


# Daytime abundance 
day_predict_data <- data.frame(diel_num = 1.5, depth = seq(1, 1000, by = 1), pca = 0)
day_predict_data

day_predict <- predict.gam(All_model1, day_predict_data, 
                           type = "response", se.fit = T)
day_predict_data <- cbind(day_predict, day_predict_data)
day_predict_data
day_predict_data$day_predict <- exp(day_predict_data$fit) -1
day_predict_data$SE <- exp(day_predict_data$se.fit) -1

day_predict_data$upr <- day_predict_data$day_predict + (2 * day_predict_data$SE)
day_predict_data$lwr <- day_predict_data$day_predict - (2 * day_predict_data$SE)
day_predict_data

peak_day <- day_predict_data$depth[which.max(day_predict_data$day_predict)]
peak_day

central_day <- sum(day_predict_data$day_predict * day_predict_data$depth)/
  sum(day_predict_data$day_predict)
central_day

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

peak_night <- night_predict_data$depth[which.max(night_predict_data$night_predict)]
peak_night

central_night <- sum(night_predict_data$night_predict * night_predict_data$depth)/
  sum(night_predict_data$night_predict)
central_night

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

intercept <- day_prop$x[which.mins(abs(day_prop$y - night_prop$y))][1]
intercept


# Plot
day_predict_data$lwr[day_predict_data$lwr < 0] <- 0
night_predict_data$lwr[night_predict_data$lwr < 0] <- 0

All_fit <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = new_IKMT[new_IKMT$depth < 1000 & 
                           new_IKMT$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = new_IKMT[new_IKMT$depth < 1000 & 
                           new_IKMT$diel_num == 4 ,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,1000,100)) + 
  labs(title = "(B)", y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = c(0.95, 1), legend.justification = c("right", "top")) +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt()  
All_fit

All_fit1 <- ggplot() +
  geom_vline(xintercept = intercept, alpha = 0.3, linewidth = 1, color = "#0066cc") +
  geom_ribbon(data = day_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Day"), alpha = 0.3) +
  geom_ribbon(data = night_predict_data, aes(x = depth, ymin=lwr, ymax=upr, fill = "Night"), alpha = 0.3) +
  geom_point(data = new_IKMT[new_IKMT$depth < 1000 & 
                               new_IKMT$diel_num == 2 ,], 
             aes(x = depth, y = CPUE, col = "Day", shape = "Day"), alpha = 0.5) +
  geom_point(data = new_IKMT[new_IKMT$depth < 1000 & 
                               new_IKMT$diel_num == 4 ,], 
             aes(x = depth, y = CPUE, col = "Night", shape = "Night"), alpha = 0.5) +
  geom_point(data = day_predict_data, aes(x = depth, y = day_predict, col = "Day")) +
  geom_point(data = day_predict_data[seq(1, nrow(day_predict_data), 50),], 
             aes(x = depth, y = inflated_day, col = "Day")) +
  geom_point(data = night_predict_data, aes(x = depth, y = night_predict, col = "Night"),
             show.legend = T) +
  theme_classic() +
  scale_x_reverse(breaks = seq(0,1000,100)) + 
  labs(title = "(B)", y = expression(Abundance~(ind.~per~'1000'~m^3)), x = "Depth (m)") +
  theme(legend.position = c(0.95, 1), legend.justification = c("right", "top")) +
  scale_fill_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                    breaks=c("Day","Night")) +
  scale_color_manual(name = "", values=c(Day="#E69F00", Night="#56B4E9"),
                     breaks=c("Day","Night")) +
  scale_shape_manual(name = "", values=c(Day=2, Night=1),
                     breaks=c("Day","Night")) +
  scale_y_sqrt() + coord_flip()
All_fit1


# Proportion of population doing DVM
proportion <- (sum(night_predict_data$night_predict) - 
                 sum(day_predict_data[day_predict_data$depth <= intercept,]$inflated_day) - 
                 sum(night_predict_data[night_predict_data$depth > intercept,]$night_predict))/
  sum(night_predict_data$night_predict)
proportion

name <- "All species"
Results <- data.frame(Species = name,
                      Peak_day = peak_day, Peak_night = peak_night,
                      centre_day = round(central_day, digits = 0), 
                      centre_night = round(central_night, digits = 0),
                      day_abundance = paste0(round(day_abundance/1000, digits = 3), 
                                             " (", round(day_abundance_SE/1000, digits = 3), ")"), 
                      night_abundance = paste0(round(night_abundance/1000, digits = 3), 
                                               " (", round(night_abundance_SE/1000, digits = 3), ")"),
                      Threshold = intercept,
                      proportion = round(proportion*100, digits = 1))

Results
Results_table <- rbind(Results_table, Results)
Results_table



####################
# Exporting result # 
####################

# Model prediction
Results_table <- na.omit(Results_table)
Results_table

mean(Results_table$proportion)

write.csv(Results_table, "DVM_pattern.csv", row.names = F)

# GAM result
GAM_result <- na.omit(GAM_result)
GAM_result

write.csv(GAM_result, "GAM_result.csv", row.names = F)


# DVM pattern
DVM_plot <- Eant_plot1 + Kand_plot1 + Gbra_plot1 + Pbol_plot1 +
  Gnic_plot1 + Gfra_plot1 + Pten_plot1 + Ecar_plot1 + plot_layout(ncol = 2) 
DVM_plot

wrap_elements(DVM_plot) +
  labs(tag = "Depth (m)") +
  theme(
    plot.tag = element_text(size = rel(1), angle = 90),
    plot.tag.position = "left"
  )
# 800w x 1200H (8x12)

# DVM pattern (day/night only)
Eant_fit +
  Kand_fit + Gbra_fit + Pbol_fit +
  Gnic_fit + Gfra_fit + Pten_fit +
  Ecar_fit +  plot_layout(ncol = 2, axis_titles = "collect")
# 800W x 1400H (8 x 14)

# DVM pattern (day/night only)
Eant_fit1 +
  Kand_fit1 + Gbra_fit1 + Pbol_fit1 +
  Gnic_fit1 + Gfra_fit1 + Pten_fit1 +
  Ecar_fit1 +  plot_layout(ncol = 2, axis_titles = "collect")
# 800W x 1400H (8 x 14)


# All species model 
All_plot1 + All_fit
# 1200W x 600H

# All species model 
All_plot1 + All_fit1
# 1200W x 600H
