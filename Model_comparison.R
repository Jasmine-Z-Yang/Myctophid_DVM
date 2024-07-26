########################
### Model comparison ###
########################

rm(list = ls())

require(dplyr)
require(mgcv)



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



# Grouping cruises
data$cruiseCode[grep("AMLR", data$cruiseCode)] <- "AMLR"
data$cruiseCode[grep("Ichtyo", data$cruiseCode)] <- "Ichtyo"
data$cruiseCode[grep("ichtyo", data$cruiseCode)] <- "Ichtyo"
data$cruiseCode[grep("JR", data$cruiseCode)] <- "BAS"
data$cruiseCode[grep("JB", data$cruiseCode)] <- "BAS"



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


knots <- list(diel_num = c(0.5,1.5,2.5,3.5,4.5))


# Result table
Comparison <- data.frame(Species = vector(), Distribution = vector(), Transformation = vector(),
                         ML = vector(), Variance = vector(), AIC = vector())
Comparison


par(mfrow = c(4,4))
par(mar=c(2,2,2,1))


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


# Tweedie distribution #

#Untransformed 
Eant_model1 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Eant, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Eant_model1)

Result <- data.frame(Species = Eant$scientificName[1], 
                     Distribution = summary(Eant_model1)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Eant_model1)$sp.criterion, 
                     Variance = summary(Eant_model1)$dev.exp, 
                     AIC = Eant_model1$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

num <- 1
Eant_qq1 <- qq.gam(Eant_model1, 
                   main = paste0("(", as.roman(num), ") ", Eant$scientificName[1], ": Tweedie"),
                   adj = 0)
num <- num + 1


#Transformed 
Eant_model2 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Eant, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Eant_model2)

Result <- data.frame(Species = Eant$scientificName[1], 
                     Distribution = summary(Eant_model2)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Eant_model2)$sp.criterion, 
                     Variance = summary(Eant_model2)$dev.exp, 
                     AIC = Eant_model2$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Eant_qq2 <- qq.gam(Eant_model2, 
                   main = paste0("(", as.roman(num), ") ", Eant$scientificName[1], ": Tweedie (log)"),
                   adj = 0)
num <- num + 1


# Negative binomial distribution #

#Untransformed 
Eant_model3 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5)  + pca, 
                   data = Eant, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Eant_model3)

Result <- data.frame(Species = Eant$scientificName[1], 
                     Distribution = summary(Eant_model3)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Eant_model3)$sp.criterion, 
                     Variance = summary(Eant_model3)$dev.exp, 
                     AIC = Eant_model3$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Eant_qq3 <- qq.gam(Eant_model3, 
                   main = paste0("(", as.roman(num), ") ", Eant$scientificName[1], ": Negative binomial"),
                   adj = 0)
num <- num + 1


#Transformed 
Eant_model4 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5)  + pca, 
                   data = Eant, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Eant_model4)

Result <- data.frame(Species = Eant$scientificName[1], 
                     Distribution = summary(Eant_model4)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Eant_model4)$sp.criterion, 
                     Variance = summary(Eant_model4)$dev.exp, 
                     AIC = Eant_model4$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Eant_qq4 <- qq.gam(Eant_model4, 
                   main = paste0("(", as.roman(num), ") ", Eant$scientificName[1], ": Negative binomial (log)"),
                   adj = 0)
num <- num + 1



############################
# Krefftichthys anderssoni #
############################

Kand <- species[[3]]
sum(Kand$individualCount)
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


# Tweedie distribution #

#Untransformed 
Kand_model1 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Kand, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Kand_model1)
Result <- data.frame(Species = Kand$scientificName[1], 
                     Distribution = summary(Kand_model1)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Kand_model1)$sp.criterion, 
                     Variance = summary(Kand_model1)$dev.exp, 
                     AIC = Kand_model1$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Kand_qq1 <- qq.gam(Kand_model1, 
                   main = paste0("(", as.roman(num), ") ", Kand$scientificName[1], ": Tweedie"),
                   adj = 0)
num <- num + 1


#Transformed 
Kand_model2 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Kand, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Kand_model2)

Result <- data.frame(Species = Kand$scientificName[1], 
                     Distribution = summary(Kand_model2)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Kand_model2)$sp.criterion, 
                     Variance = summary(Kand_model2)$dev.exp, 
                     AIC = Kand_model2$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Kand_qq2 <- qq.gam(Kand_model2, 
                   main = paste0("(", as.roman(num), ") ", Kand$scientificName[1], ": Tweedie (log)"),
                   adj = 0)
num <- num + 1


# Negative binomial distribution #

#Untransformed 
Kand_model3 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5)  + pca, 
                   data = Kand, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Kand_model3)

Result <- data.frame(Species = Kand$scientificName[1], 
                     Distribution = summary(Kand_model3)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Kand_model3)$sp.criterion, 
                     Variance = summary(Kand_model3)$dev.exp, 
                     AIC = Kand_model3$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Kand_qq3 <- qq.gam(Kand_model3, 
                   main = paste0("(", as.roman(num), ") ", Kand$scientificName[1], ": Negative binomial"),
                   adj = 0)
num <- num + 1


#Transformed 
Kand_model4 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Kand, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Kand_model4)
Result <- data.frame(Species = Kand$scientificName[1], 
                     Distribution = summary(Kand_model4)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Kand_model4)$sp.criterion, 
                     Variance = summary(Kand_model4)$dev.exp, 
                     AIC = Kand_model4$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Kand_qq4 <- qq.gam(Kand_model4, 
                   main = paste0("(", as.roman(num), ") ", Kand$scientificName[1], ": Negative binomial (log)"),
                   adj = 0)
num <- num + 1



#########################
# Gymnoscopelus braueri #
#########################

Gbra <- species[[5]]
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



# Tweedie distribution #

#Untransformed 
Gbra_model1 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gbra, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Gbra_model1)

Result <- data.frame(Species = Gbra$scientificName[1], 
                     Distribution = summary(Gbra_model1)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Gbra_model1)$sp.criterion, 
                     Variance = summary(Gbra_model1)$dev.exp, 
                     AIC = Gbra_model1$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gbra_qq1 <- qq.gam(Gbra_model1, 
                   main = paste0("(", as.roman(num), ") ", Gbra$scientificName[1], ": Tweedie"),
                   adj = 0)
num <- num + 1


#Transformed 
Gbra_model2 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gbra, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Gbra_model2)

Result <- data.frame(Species = Gbra$scientificName[1], 
                     Distribution = summary(Gbra_model2)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Gbra_model2)$sp.criterion, 
                     Variance = summary(Gbra_model2)$dev.exp, 
                     AIC = Gbra_model2$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gbra_qq2 <- qq.gam(Gbra_model2, 
                   main = paste0("(", as.roman(num), ") ", Gbra$scientificName[1], ": Tweedie (log)"),
                   adj = 0)
num <- num + 1



# Negative binomial distribution #

#Untransformed 
Gbra_model3 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gbra, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Gbra_model3)
Result <- data.frame(Species = Gbra$scientificName[1], 
                     Distribution = summary(Gbra_model3)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Gbra_model3)$sp.criterion, 
                     Variance = summary(Gbra_model3)$dev.exp, 
                     AIC = Gbra_model3$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gbra_qq3 <- qq.gam(Gbra_model3, 
                   main = paste0("(", as.roman(num), ") ", Gbra$scientificName[1], ": Negative binomial"),
                   adj = 0)
num <- num + 1


#Transformed 
Gbra_model4 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gbra, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Gbra_model4)

Result <- data.frame(Species = Gbra$scientificName[1], 
                     Distribution = summary(Gbra_model4)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Gbra_model4)$sp.criterion, 
                     Variance = summary(Gbra_model4)$dev.exp, 
                     AIC = Gbra_model4$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gbra_qq4 <- qq.gam(Gbra_model4, 
                   main = paste0("(", as.roman(num), ") ", Gbra$scientificName[1], ": Negative binomial (log)"),
                   adj = 0)
num <- num + 1



#########################
# Protomyctophum bolini #
#########################

Pbol <- species[[7]]
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



# Tweedie distribution #

#Untransformed 
Pbol_model1 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pbol, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Pbol_model1)

Result <- data.frame(Species = Pbol$scientificName[1], 
                     Distribution = summary(Pbol_model1)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Pbol_model1)$sp.criterion, 
                     Variance = summary(Pbol_model1)$dev.exp, 
                     AIC = Pbol_model1$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Pbol_qq1 <- qq.gam(Pbol_model1, 
                   main = paste0("(", as.roman(num), ") ", Pbol$scientificName[1], ": Tweedie"),
                   adj = 0)
num <- num + 1


#Transformed 
Pbol_model2 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pbol, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Pbol_model2)

Result <- data.frame(Species = Pbol$scientificName[1], 
                     Distribution = summary(Pbol_model2)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Pbol_model2)$sp.criterion, 
                     Variance = summary(Pbol_model2)$dev.exp, 
                     AIC = Pbol_model2$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Pbol_qq2 <- qq.gam(Pbol_model2, 
                   main = paste0("(", as.roman(num), ") ", Pbol$scientificName[1], ": Tweedie (log)"),
                   adj = 0)
num <- num + 1



# Negative binomial distribution #

#Untransformed 
Pbol_model3 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pbol, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Pbol_model3)

Result <- data.frame(Species = Pbol$scientificName[1], 
                     Distribution = summary(Pbol_model3)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Pbol_model3)$sp.criterion, 
                     Variance = summary(Pbol_model3)$dev.exp, 
                     AIC = Pbol_model3$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Pbol_qq3 <- qq.gam(Pbol_model3, 
                   main = paste0("(", as.roman(num), ") ", Pbol$scientificName[1], ": Negative binomial"),
                   adj = 0)
num <- num + 1


#Transformed 
Pbol_model4 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pbol, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Pbol_model4)

Result <- data.frame(Species = Pbol$scientificName[1], 
                     Distribution = summary(Pbol_model4)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Pbol_model4)$sp.criterion, 
                     Variance = summary(Pbol_model4)$dev.exp, 
                     AIC = Pbol_model4$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Pbol_qq4 <- qq.gam(Pbol_model4, 
                   main = paste0("(", as.roman(num), ") ", Pbol$scientificName[1], ": Negative binomial (log)"),
                   adj = 0)
num <- num + 1



##########################
# Gymnoscopelus nicholsi #
##########################

Gnic <- species[[8]]
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

Gnic_Sub <- Gnic[Gnic$lat >= f(count$count)[[1]], ]


# Net removal 
Gnic <- Gnic[Gnic$netType != "IKMT",]
unique(Gnic$netType)



# Tweedie distribution #

#Untransformed 
Gnic_model1 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gnic, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Gnic_model1)

Result <- data.frame(Species = Gnic$scientificName[1], 
                     Distribution = summary(Gnic_model1)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Gnic_model1)$sp.criterion, 
                     Variance = summary(Gnic_model1)$dev.exp, 
                     AIC = Gnic_model1$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gnic_qq1 <- qq.gam(Gnic_model1, 
                   main = paste0("(", as.roman(num), ") ", Gnic$scientificName[1], ": Tweedie"),
                   adj = 0)
num <- num + 1


#Transformed 
Gnic_model2 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gnic, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Gnic_model2)

Result <- data.frame(Species = Gnic$scientificName[1], 
                     Distribution = summary(Gnic_model2)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Gnic_model2)$sp.criterion, 
                     Variance = summary(Gnic_model2)$dev.exp, 
                     AIC = Gnic_model2$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gnic_qq2 <- qq.gam(Gnic_model2, 
                   main = paste0("(", as.roman(num), ") ", Gnic$scientificName[1], ": Tweedie (log)"),
                   adj = 0)
num <- num + 1



# Negative binomial distribution #

#Untransformed 
Gnic_model3 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gnic, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Gnic_model3)

Result <- data.frame(Species = Gnic$scientificName[1], 
                     Distribution = summary(Gnic_model3)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Gnic_model3)$sp.criterion, 
                     Variance = summary(Gnic_model3)$dev.exp, 
                     AIC = Gnic_model3$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gnic_qq3 <- qq.gam(Gnic_model3, 
                   main = paste0("(", as.roman(num), ") ", Gnic$scientificName[1], ": Negative binomial"),
                   adj = 0)
num <- num + 1


#Transformed 
Gnic_model4 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gnic, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Gnic_model4)

Result <- data.frame(Species = Gnic$scientificName[1], 
                     Distribution = summary(Gnic_model4)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Gnic_model4)$sp.criterion, 
                     Variance = summary(Gnic_model4)$dev.exp, 
                     AIC = Gnic_model4$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gnic_qq4 <- qq.gam(Gnic_model4, 
                   main = paste0("(", as.roman(num), ") ", Gnic$scientificName[1], ": Negative binomial (log)"),
                   adj = 0)
num <- num + 1



##########################
# Gymnoscopelus fraseri #
##########################

Gfra <- species[[10]]
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


# Tweedie distribution #

#Untransformed 
Gfra_model1 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gfra, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Gfra_model1)

Result <- data.frame(Species = Gfra$scientificName[1], 
                     Distribution = summary(Gfra_model1)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Gfra_model1)$sp.criterion, 
                     Variance = summary(Gfra_model1)$dev.exp, 
                     AIC = Gfra_model1$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gfra_qq1 <- qq.gam(Gfra_model1, 
                   main = paste0("(", as.roman(num), ") ", Gfra$scientificName[1], ": Tweedie"),
                   adj = 0)
num <- num + 1


#Transformed 
Gfra_model2 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gfra, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Gfra_model2)

Result <- data.frame(Species = Gfra$scientificName[1], 
                     Distribution = summary(Gfra_model2)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Gfra_model2)$sp.criterion, 
                     Variance = summary(Gfra_model2)$dev.exp, 
                     AIC = Gfra_model2$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gfra_qq2 <- qq.gam(Gfra_model2, 
                   main = paste0("(", as.roman(num), ") ", Gfra$scientificName[1], ": Tweedie (log)"),
                   adj = 0)
num <- num + 1


# Negative binomial distribution #

#Untransformed 
Gfra_model3 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gfra, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Gfra_model3)

Result <- data.frame(Species = Gfra$scientificName[1], 
                     Distribution = summary(Gfra_model3)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Gfra_model3)$sp.criterion, 
                     Variance = summary(Gfra_model3)$dev.exp, 
                     AIC = Gfra_model3$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gfra_qq3 <- qq.gam(Gfra_model3, 
                   main = paste0("(", as.roman(num), ") ", Gfra$scientificName[1], ": Negative binomial"),
                   adj = 0)
num <- num + 1


#Transformed 
Gfra_model4 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Gfra, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Gfra_model4)

Result <- data.frame(Species = Gfra$scientificName[1], 
                     Distribution = summary(Gfra_model4)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Gfra_model4)$sp.criterion, 
                     Variance = summary(Gfra_model4)$dev.exp, 
                     AIC = Gfra_model4$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Gfra_qq4 <- qq.gam(Gfra_model4, 
                   main = paste0("(", as.roman(num), ") ", Gfra$scientificName[1], ": Negative binomial (log)"),
                   adj = 0)
num <- num + 1




###########################
# Protomyctophum tenisoni #
###########################

Pten <- species[[11]]
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


# Tweedie distribution #

#Untransformed 
Pten_model1 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pten, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Pten_model1)

Result <- data.frame(Species = Pten$scientificName[1], 
                     Distribution = summary(Pten_model1)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Pten_model1)$sp.criterion, 
                     Variance = summary(Pten_model1)$dev.exp, 
                     AIC = Pten_model1$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Pten_qq1 <- qq.gam(Pten_model1, 
                   main = paste0("(", as.roman(num), ") ", Pten$scientificName[1], ": Tweedie"),
                   adj = 0)
num <- num + 1


#Transformed 
Pten_model2 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pten, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Pten_model2)

Result <- data.frame(Species = Pten$scientificName[1], 
                     Distribution = summary(Pten_model2)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Pten_model2)$sp.criterion, 
                     Variance = summary(Pten_model2)$dev.exp, 
                     AIC = Pten_model2$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Pten_qq2 <- qq.gam(Pten_model2, 
                   main = paste0("(", as.roman(num), ") ", Pten$scientificName[1], ": Tweedie (log)"),
                   adj = 0)
num <- num + 1


# Negative binomial distribution #

#Untransformed 
Pten_model3 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pten, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Pten_model3)

Result <- data.frame(Species = Pten$scientificName[1], 
                     Distribution = summary(Pten_model3)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Pten_model3)$sp.criterion, 
                     Variance = summary(Pten_model3)$dev.exp, 
                     AIC = Pten_model3$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Pten_qq3 <- qq.gam(Pten_model3, 
                   main = paste0("(", as.roman(num), ") ", Pten$scientificName[1], ": Negative binomial"),
                   adj = 0)
num <- num + 1


#Transformed 
Pten_model4 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Pten, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Pten_model4)

Result <- data.frame(Species = Pten$scientificName[1], 
                     Distribution = summary(Pten_model4)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Pten_model4)$sp.criterion, 
                     Variance = summary(Pten_model4)$dev.exp, 
                     AIC = Pten_model4$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Pten_qq4 <- qq.gam(Pten_model4, 
                   main = paste0("(", as.roman(num), ") ", Pten$scientificName[1], ": Negative binomial (log)"),
                   adj = 0)
num <- num + 1



########################
# Electrona carlsbergi #
########################

Ecar <- species[[12]]
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


Ecar_raw <- Ecar %>%
  group_by(diel_num, depth) %>%
  summarise(n = n(),
            mean = mean(CPUE))
head(Ecar_raw)


# Tweedie distribution #

#Untransformed 
Ecar_model1 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Ecar, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Ecar_model1)
Result <- data.frame(Species = Ecar$scientificName[1], 
                     Distribution = summary(Ecar_model1)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Ecar_model1)$sp.criterion, 
                     Variance = summary(Ecar_model1)$dev.exp, 
                     AIC = Ecar_model1$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Ecar_qq1 <- qq.gam(Ecar_model1, 
                   main = paste0("(", as.roman(num), ") ", Ecar$scientificName[1], ": Tweedie"),
                   adj = 0)
num <- num + 1


#Transformed 
Ecar_model2 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Ecar, select = TRUE,
                   family = "tw", knots = knots, method = "ML")
summary(Ecar_model2)

Result <- data.frame(Species = Ecar$scientificName[1], 
                     Distribution = summary(Ecar_model2)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Ecar_model2)$sp.criterion, 
                     Variance = summary(Ecar_model2)$dev.exp, 
                     AIC = Ecar_model2$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Ecar_qq2 <- qq.gam(Ecar_model2, 
                   main = paste0("(", as.roman(num), ") ", Ecar$scientificName[1], ": Tweedie (log)"),
                   adj = 0)
num <- num + 1



# Negative binomial distribution #

#Untransformed 
Ecar_model3 <- gam(CPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Ecar, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Ecar_model3)

Result <- data.frame(Species = Ecar$scientificName[1], 
                     Distribution = summary(Ecar_model3)$family[[1]], 
                     Transformation = "No",
                     ML = summary(Ecar_model3)$sp.criterion, 
                     Variance = summary(Ecar_model3)$dev.exp, 
                     AIC = Ecar_model3$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Ecar_qq3 <- qq.gam(Ecar_model3, 
                   main = paste0("(", as.roman(num), ") ", Ecar$scientificName[1], ": Negative binomial"),
                   adj = 0)
num <- num + 1


#Transformed 
Ecar_model4 <- gam(logCPUE ~ ti(diel_num, depth, k =  c(5,5), bs = c("cc", "tp")) + 
                     s(diel_num, k = 5, bs = "cc") + s(depth, k = 5) + pca, 
                   data = Ecar, select = TRUE,
                   family = "nb", knots = knots, method = "ML")
summary(Ecar_model4)

Result <- data.frame(Species = Ecar$scientificName[1], 
                     Distribution = summary(Ecar_model4)$family[[1]], 
                     Transformation = "Log",
                     ML = summary(Ecar_model4)$sp.criterion, 
                     Variance = summary(Ecar_model4)$dev.exp, 
                     AIC = Ecar_model4$aic)
Comparison <- rbind(Comparison, Result)
tail(Comparison)

Ecar_qq4 <- qq.gam(Ecar_model4, 
                   main = paste0("(", as.roman(num), ") ", Ecar$scientificName[1], ": Negative binomial (log)"),
                   adj = 0)
num <- num + 1



##################
# Export results # 
##################

# 2000W x 1000H for figures
par(mfrow = c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))


# Model comparison parameters 
write.csv(Comparison, "Model_comparison.csv", row.names = F)




