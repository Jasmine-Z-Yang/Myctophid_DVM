###################################
### Distribution raincloud plot ###
###################################

rm(list = ls())

require(ggplot2)
require(ggfx)
require(dplyr)
require(patchwork)
require(orsifronts)



## Ocean fronts boundaries ##
mean(as.data.frame(coordinates(orsifronts@lines[[1]]@Lines$'1'))$V2)
mean(as.data.frame(coordinates(orsifronts@lines[[2]]@Lines$'1'))$V2)
mean(as.data.frame(coordinates(orsifronts@lines[[3]]@Lines$'1'))$V2)
mean(as.data.frame(coordinates(orsifronts@lines[[4]]@Lines$'1'))$V2)
mean(as.data.frame(coordinates(orsifronts@lines[[5]]@Lines$'1'))$V2)



####################
## Data wrangling ##
####################

# Loading data #
event <- as.data.frame(read.csv("Myctobase/event_edited2.csv", header = TRUE, stringsAsFactors = F))

group <- as.data.frame(read.csv("Myctobase/groupOccurrence.csv", header = TRUE, stringsAsFactors = F))



# Choosing study species #

table <- table(group$scientificName) 
table

order <- order(table, decreasing = TRUE)
table[order[1:20]]

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
  new <- merge(event, species[[i]], by = "eventID", all.x = TRUE, all.y = FALSE)
  
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
  data$CPUE <- (data$individualCount/data$volume)*1000
  data$logCPUE <- log(data$CPUE + 1)
  
  species[[i]] <- data
  
}



####################
## Mean abundance ##
####################

# Species mean abundance per latitude 

all_catch <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
                      c("lat", "min", "median", "max", "mean", "species"))
all_catch

name_order <- vector()
for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  name_order <- c(name_order, name)
  
  data <- species[[i]]

  data <- data[, c("CPUE", "lat")]
  
  data <- data %>%
    group_by(lat) %>%
    summarise(min = quantile(CPUE, probs = 0.1),
              median = quantile(CPUE, probs = 0.5),
              max = quantile(CPUE, probs = 0.9),
              mean = mean(CPUE))
  data$species <- name
  
  print(head(data))
  all_catch <- rbind(all_catch, data)
  
}

View(all_catch)

# Making plot 


max(all_catch$mean)

library(rcartocolor)
display_carto_pal(9, "Safe")
palette <- carto_pal(9, "Safe")

choice <- 1
for(i in unique(all_catch$species)){
  print(i)
  
  header <- LETTERS[choice]
  title <- bquote("("*.(header)*")"~italic(.(i)))
  
  assign(paste0("CPUE_",i), ggplot(all_catch[all_catch$species == i,], 
                                   aes(x = lat, y = mean, group = lat)) + 
           with_blur(geom_rect(data = all_catch, 
                     aes(xmin = -Inf, 
                         xmax = -65,
                         ymin = -Inf,
                         ymax = Inf),
                     alpha = 1, fill = "#E7E7E7"), 
                     sigma = 5) + 
           with_blur(geom_rect(data = all_catch,
                     aes(xmin = -55, 
                         xmax = -51,
                         ymin = -Inf,
                         ymax = Inf),
                     alpha = 1, fill = "#E7E7E7"), 
                     sigma = 5) +
           geom_bar(stat = "identity", fill = palette[choice]) +
           scale_x_continuous(breaks = seq(-70, -48, by = 2), limits = c(-70.5, -47.5), 
                              expand = c(0.01, 0.01)) +
           scale_y_continuous(limits = c(0,max(all_catch[all_catch$species == i,]$mean)),
                              expand = c(0.01,0)) +
           labs(x = "Latitude (°)", y = expression("Mean abundance per 10000m"^-3*"")) +
           theme_classic() +
           theme(legend.position = "none",
                 strip.text.y = element_text(angle=0, face = "italic"),
                 strip.background = element_blank(),
                 axis.line.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),
                 plot.margin = margin(0.1,0,0,0, unit = "points"),
                 plot.title = element_text(size = 10)) +
           coord_cartesian(clip = 'off') +
           ggtitle(title) )
  print(get(paste0("CPUE_",i)))
  
  choice <- choice + 1
}


ocean_text1 <- data.frame(
  label = c("Subpolar Zone"),
  species   = c("Electrona antarctica"),
  x     = c(-68),
  y     = c(Inf)
)

ocean_text2 <- data.frame(
  label = c("Antarctic\nZone","Polar Front\nZone","Subantarctic\nZone"),
  species   = c(rep("Electrona antarctica",3)),
  x     = c(-63,-53,-49),
  y     = c(rep(Inf, 3))
)

ocean_text3 <- data.frame(
  label = c("sACCF","PF","SAF"),
  species   = c(rep("Electrona antarctica",3)),
  x     = c(-65,-55,-51),
  y     = c(rep(Inf, 3))
)


i <- unique(all_catch$species)[1]
choice <- 1
header <- LETTERS[choice]
title <- bquote("("*.(header)*")"~italic(.(i)))
assign(paste0("CPUE_",i), ggplot(all_catch[all_catch$species == i,], 
                                 aes(x = lat, y = mean, group = lat)) + 
         with_blur(geom_rect(data = all_catch, 
                             aes(xmin = -Inf, 
                                 xmax = -65,
                                 ymin = -Inf,
                                 ymax = Inf),
                             alpha = 1, fill = "#E7E7E7"), 
                   sigma = 5) + 
         with_blur(geom_rect(data = all_catch,
                             aes(xmin = -55, 
                                 xmax = -51,
                                 ymin = -Inf,
                                 ymax = Inf),
                             alpha = 1, fill = "#E7E7E7"), 
                   sigma = 5) +
         geom_text(data = ocean_text1, inherit.aes=FALSE,
                   mapping = aes(x = x, y = y, label = label),
                   fontface = "bold",
                   col = "gray", size = 4, vjust = 2, position = position_stack()) +
         geom_text(data = ocean_text2, inherit.aes=FALSE,
                   mapping = aes(x = x, y = y, label = label),
                   fontface = "bold",
                   col = "gray", size = 4, vjust = 1.2, position = position_stack()) +
         geom_text(data = ocean_text3, inherit.aes=FALSE,
                   mapping = aes(x = x, y = y, label = label),
                   col = "black", size = 3, vjust = -1, position = position_stack()) +
         geom_bar(stat = "identity", fill = palette[choice]) +
         scale_x_continuous(breaks = seq(-70, -48, by = 2), limits = c(-70.5, -47.5), 
                            expand = c(0.01, 0.01)) +
         scale_y_continuous(limits = c(0,max(all_catch[all_catch$species == i,]$mean)),
                            expand = c(0.01,0)) +
         labs(x = "Latitude (°)", y = expression("Mean abundance per 10000m"^-3*"")) +
         theme_classic() +
         theme(legend.position = "none",
               strip.text.y = element_text(angle=0, face = "italic"),
               strip.background = element_blank(),
               axis.line.x = element_blank(),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               plot.margin = margin(0.1,0,0,0, unit = "points"),
               plot.title = element_text(size = 10)) +
         coord_cartesian(clip = 'off') +
         ggtitle(title) )



####################
## Presence count ##
####################

# Defining quantiles of boxplot #

f <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


all_count <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("species", "count"))
head(all_count)

name_order <- vector()


# Making dataframe of presence/absence

for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  name_order <- c(name_order, name)
  
  new <- species[[i]]
  
  
  # Making column of presence/absence 
  new$presence <- NA
  
  for(i in 1:nrow(new)){
    if(new$individualCount[i] > 0){
      new$presence[i] <- 1
    } else {new$presence[i] <- 0}
  }
  
  
  # Distribution of presence 
  new$presence <- as.factor(new$presence)
  
  presence_count <- new %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
  
  count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
  
  count <- vector()
  
  for(i in 1:nrow(count_table)){
    if(count_table$`1`[i] > 0){
      count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
    }
  }
  
  print(f(count))
  
  ind_count <- data.frame(species = rep(name, length(count)), count = count)
  all_count <- rbind(all_count, ind_count)
  
}

name_order <- NA
for(i in species_choice){
  
  name <- species[[i]]$scientificName[1]
  print(name)
  
  name_order <- c(name_order, name)
  
  new <- species[[i]]
  
  
  # Making column of presence/absence 
  new$presence <- NA
  
  for(i in 1:nrow(new)){
    if(new$individualCount[i] > 0){
      new$presence[i] <- 1
    } else {new$presence[i] <- 0}
  }
  
  
  # Distribution of presence 
  new$presence <- as.factor(new$presence)
  
  presence_count <- new %>% count(lat, presence, .drop = FALSE) # Making data frame of presence distribution 
  
  total <- presence_count %>% group_by(lat) %>% summarise(total = sum(n))
  total$presence <- presence_count[presence_count$presence == 1,]$n/total$total
  print(ggplot(data = total, aes(x = lat, y = presence)) + geom_bar(stat = "identity",))
  
  count_table <- as.data.frame(presence_count %>% tidyr::spread(presence, n))
  
  count <- vector()
  
  for(i in 1:nrow(count_table)){
    if(count_table$`1`[i] > 0){
      count <- c(count, rep(count_table$lat[i], count_table$`1`[i]))
    }
  }
  
  print(f(count))
  
  ind_count <- data.frame(species = rep(name, length(count)), count = count)
  all_count <- rbind(all_count, ind_count)
  
}

head(all_count)
name_order


# Making plot 

all_count$species <- factor(all_count$species, levels = name_order)

choice <- 1
for(i in unique(all_count$species)){
  print(i)
  
  assign(paste0("distribution_",i), ggplot(all_count[all_count$species == i,], aes(x=NA, y=count, fill = species)) + 
    stat_summary(fun.data = f, geom="boxplot") +
      scale_fill_manual(values = c(palette[choice])) +
    coord_flip() + 
    scale_y_continuous(breaks = seq(-70, -48, by = 2), limits = c(-70.5, -47.5), 
                       expand = c(0.01, 0.01)) +
    labs(y = "Latitude (°)") +
    theme_classic() +
    theme(legend.position = "none",
          axis.line.y = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(0,0,5,0, unit = "points"),
          plot.title = element_text(size = 10),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA)) +
      scale_x_discrete(expand = c(0.5,0.5)))
  
  print(get(paste0("distribution_",i)))
  
  choice <- choice + 1
}


pdf("Species distribution.pdf", height = 12, width = 8)

`CPUE_Electrona antarctica` + plot_spacer() + `distribution_Electrona antarctica` + 
`CPUE_Krefftichthys anderssoni` + plot_spacer() + `distribution_Krefftichthys anderssoni` + 
`CPUE_Gymnoscopelus braueri` + plot_spacer() + `distribution_Gymnoscopelus braueri` + 
`CPUE_Protomyctophum bolini` + plot_spacer() + `distribution_Protomyctophum bolini` + 
`CPUE_Gymnoscopelus nicholsi` + plot_spacer() + `distribution_Gymnoscopelus nicholsi` + 
`CPUE_Gymnoscopelus fraseri` + plot_spacer() + `distribution_Gymnoscopelus fraseri` + 
`CPUE_Protomyctophum tenisoni` + plot_spacer() + `distribution_Protomyctophum tenisoni` + 
`CPUE_Electrona carlsbergi` + plot_spacer() + `distribution_Electrona carlsbergi` + 
  plot_layout(ncol = 1, heights = rep(c(4, -1.1, 1),8))

dev.off()

