###########################
### Fish size variation ###
###########################

rm(list = ls())

require(ggplot2)
require(patchwork)



# Loading data #
data <- as.data.frame(read.csv("Myctobase/event_edited2.csv", header = TRUE, stringsAsFactors = F))
str(data)
length(unique(data$eventID))

group <- as.data.frame(read.csv("Myctobase/groupOccurrence.csv", header = TRUE, stringsAsFactors = F))
str(group)
length(unique(group$eventID))

og_ind <- as.data.frame(read.csv("Myctobase/individualOccurrence.csv", header = TRUE, stringsAsFactors = F))
str(og_ind)
length(unique(og_ind$eventID))


sp_list <- c("Electrona antarctica", "Electrona carlsbergi", 
             "Gymnoscopelus braueri", 
             "Gymnoscopelus fraseri", "Gymnoscopelus nicholsi", 
             "Krefftichthys anderssoni", "Protomyctophum bolini",
             "Protomyctophum tenisoni")
sp_list



###################################
# Histogram of size caught vs net #
###################################

library(rcartocolor)
display_carto_pal(9, "Safe")
palette <- carto_pal(6, "Safe")

choice <- 1
for(i in sp_list){
  print(paste("Species:", i))
  
  sp_size <- og_ind[og_ind$scientificName == i, ] 
  sp_size <- sp_size[!is.na(sp_size$standard_length),]
  
  
  sp <-  merge(data, sp_size, by = "eventID", all.x = FALSE, all.y = TRUE)
  print(dim(sp))
  
  
  # Filling data from literature 
  sp$standard_length[sp$standard_length == "SC 1"] <- 10
  sp$standard_length[sp$standard_length == "SC 2"] <- 30
  sp$standard_length[sp$standard_length == "SC 3"] <- 50
  sp <- sp[!is.na(as.numeric(sp$standard_length)),]
  sp$standard_length <- as.numeric(sp$standard_length)
  
  print(unique(sp$netType))
  sp$netType <- factor(sp$netType, 
                          levels = c("IYGPT (171)", 
                                     "IYGPT with MIDOC (188)",
                                     "MOHT",
                                     "RMT25",
                                     "RMT8"))
  
  print(paste("Mean:", mean(sp$standard_length)))
  print(paste("Median:", median(sp$standard_length)))
  print(paste("Range:", range(sp$standard_length)))
  
  # Size variation by net 
  header <- LETTERS[choice]
  title <- bquote("("*.(header)*")"~italic(.(i)))
  assign(paste0("size_",i), ggplot(data = sp[!is.na(sp$netType),], aes(x = standard_length, 
                                                                       fill = netType)) +
           geom_histogram(alpha = 0.5, position = "identity") +
           theme_classic() +
           xlim(0,165) +
           ggtitle(title) +
           ylab("Count") + xlab("Standard length (mm)") +
           labs(fill = "Net type") +
           scale_fill_manual(values = c("IYGPT (171)" = palette[1], 
                                          "IYGPT with MIDOC (188)" = palette[4],
                                          "MOHT" = palette[5],
                                          "RMT25" = palette[2],
                                          "RMT8" = palette[3]), drop = FALSE))
  print(get(paste0("size_",i)))
  choice <- choice + 1
}
  

hist <- `size_Electrona antarctica` + `size_Krefftichthys anderssoni` + `size_Gymnoscopelus braueri` + 
  `size_Protomyctophum bolini` + `size_Gymnoscopelus nicholsi` + `size_Gymnoscopelus fraseri` +
  `size_Protomyctophum tenisoni` + `size_Electrona carlsbergi` + plot_layout(ncol = 2) & 
  theme(legend.position = "bottom")
hist + plot_layout(guides = "collect", axis_titles = "collect") 
# 10W x 14H



#####################################################
### Effect of dpeth/lat on size of samples caught ###
#####################################################

library(rcartocolor)
display_carto_pal(6, "Safe")
palette <- carto_pal(6, "Safe")

max(data$lat)
min(data$lat)
max(data$depth)
min(data$depth)


choice <- 1
for(i in sp_list){
  print(paste("Species:", i))
  
  sp_size <- og_ind[og_ind$scientificName == i, ] 
  sp_size <- sp_size[!is.na(sp_size$standard_length),]
  
  
  sp <-  merge(data, sp_size, by = "eventID", all.x = FALSE, all.y = TRUE)
  print(dim(sp))
  
  # Filling data from literature
  sp$standard_length[sp$standard_length == "SC 1"] <- 10
  sp$standard_length[sp$standard_length == "SC 2"] <- 30
  sp$standard_length[sp$standard_length == "SC 3"] <- 50
  sp <- sp[!is.na(as.numeric(sp$standard_length)),]
  sp$standard_length <- as.numeric(sp$standard_length)
  
  print(unique(sp$netType))
  sp$netType <- factor(sp$netType, 
                       levels = c("IYGPT (171)", 
                                  "IYGPT with MIDOC (188)",
                                  "MOHT",
                                  "RMT25",
                                  "RMT8"))
  
  print(paste("Mean:", mean(sp$standard_length)))
  print(paste("Median:", median(sp$standard_length)))
  print(paste("Range:", range(sp$standard_length)))
  
  # Size variation by net 
  header <- LETTERS[choice]
  title <- bquote("("*.(header)*")"~italic(.(i)))
  
  assign(paste0("sizebylatnet_",i), ggplot(data = sp[!is.na(sp$lat) & !is.na(sp$netType),], 
                                           aes(x = lat, y = standard_length, col = netType,
                                               group = interaction(lat, netType))) +
           geom_point(position=position_jitter(h=0.1), alpha = 0.4) +
           theme_classic() +
           xlim(-70,-48)+
           labs(col = "Net type", fill = "Net type",
                y = "Standard length (mm)", x = "Latitude (°)") +
           ggtitle(title) +
    scale_color_manual(values = c("IYGPT (171)" = palette[1], 
                                 "IYGPT with MIDOC (188)" = palette[4],
                                 "MOHT" = palette[5],
                                 "RMT25" = palette[2],
                                 "RMT8" = palette[3]), drop = FALSE))
  print(get(paste0("sizebylatnet_",i)))
  
  
  assign(paste0("sizebydepth_",i), ggplot(data = sp[!is.na(sp$depth) & !is.na(sp$netType),], 
                                           aes(x = depth, y = standard_length, col = netType,
                                               group = interaction(depth, netType))) +
           geom_point(position=position_jitter(h=0.1), alpha = 0.4) +
           theme_classic() +
           scale_x_continuous(breaks=seq(0,1000,200)) +
           labs(col = "Net type", fill = "Net type",
                y = "Standard length (mm)", x = "Depth (m)") +
           ggtitle(title) +
           scale_color_manual(values = c("IYGPT (171)" = palette[1], 
                                        "IYGPT with MIDOC (188)" = palette[4],
                                        "MOHT" = palette[5],
                                        "RMT25" = palette[2],
                                        "RMT8" = palette[3]), drop = FALSE))
  print(get(paste0("sizebydepth_",i)))
  
  
  assign(paste0("sizebydepth2_",i), ggplot(data = sp[!is.na(sp$depth) & !is.na(sp$netType),], 
                                          aes(y = depth, x = standard_length, col = netType,
                                              group = interaction(depth, netType))) +
           geom_point(position=position_jitter(h=0.1), alpha = 0.4) +
           theme_classic() +
           scale_y_reverse(expand = c(0.1, 0.1), breaks=seq(0,1000,200)) +
           labs(col = "Net type", fill = "Net type", y = "Depth (m)", x = "Standard length (mm)") +
           ggtitle(title) +
           scale_color_manual(values = c("IYGPT (171)" = palette[1], 
                                        "IYGPT with MIDOC (188)" = palette[4],
                                        "MOHT" = palette[5],
                                        "RMT25" = palette[2],
                                        "RMT8" = palette[3]), drop = FALSE))
  print(get(paste0("sizebydepth2_",i)))

  choice <- choice + 1
}


# Size vs latitude plot 
lat_plot <- `sizebylatnet_Electrona antarctica` + `sizebylatnet_Krefftichthys anderssoni` +
  `sizebylatnet_Gymnoscopelus braueri` + `sizebylatnet_Protomyctophum bolini` +
  `sizebylatnet_Gymnoscopelus nicholsi` + `sizebylatnet_Gymnoscopelus fraseri` +
  `sizebylatnet_Protomyctophum tenisoni` + `sizebylatnet_Electrona carlsbergi`+ plot_layout(ncol = 2) & 
  theme(legend.position = "bottom")

lat_plot + plot_layout(guides = "collect", axis_titles = "collect") 
# 10W x 14H


# Size vs depth plot (depth on y-axis)
depth_plot <- `sizebydepth_Electrona antarctica` + `sizebydepth_Krefftichthys anderssoni` +
  `sizebydepth_Gymnoscopelus braueri` + `sizebydepth_Protomyctophum bolini` +
  `sizebydepth_Gymnoscopelus nicholsi` + `sizebydepth_Gymnoscopelus fraseri` +
  `sizebydepth_Protomyctophum tenisoni` + `sizebydepth_Electrona carlsbergi` +
  plot_layout(ncol = 2) & theme(legend.position = "bottom")

depth_plot + plot_layout(guides = "collect", axis_titles = "collect") 
# 10W x 14H


# Size vs depth plot (depth on x-axis)
depth_plot2 <- `sizebydepth2_Electrona antarctica` + `sizebydepth2_Krefftichthys anderssoni` +
  `sizebydepth2_Gymnoscopelus braueri` + `sizebydepth2_Protomyctophum bolini` +
  `sizebydepth2_Gymnoscopelus nicholsi` + `sizebydepth2_Gymnoscopelus fraseri` +
  `sizebydepth2_Protomyctophum tenisoni` + `sizebydepth2_Electrona carlsbergi` +
  plot_layout(ncol = 2) & theme(legend.position = "bottom")

depth_plot2 + plot_layout(guides = "collect", axis_titles = "collect") 



