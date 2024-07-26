#######################
### Abundance count ###
#######################

rm(list = ls())

require(dplyr)

# Loading data #
data <- as.data.frame(read.csv("Myctobase/event_edited2.csv", header = TRUE, stringsAsFactors = F))
str(data)
length(unique(data$eventID))

og_group <- as.data.frame(read.csv("Myctobase/groupOccurrence.csv", header = TRUE, stringsAsFactors = F))
str(og_group)
length(unique(og_group$eventID))

group <- og_group[c("eventID", "family", "scientificName", "organismQuantity", "individualCount")]
str(group)



# Selecting Myctophidae family identified to species level #
Myctophidae <- group[group$family == "Myctophidae",]

Myctophidae_species <- Myctophidae[sapply(strsplit(Myctophidae$scientificName, " "), length) == 2,]


# Total catch and number of hauls present before data quality control #
abundance_before <- Myctophidae_species %>%
  group_by(scientificName) %>%
  summarise(haul.before = n(),
            total.before = sum(individualCount)) 
abundance_before

abundance_before <- na.omit(abundance_before)
abundance_before


# Total catch and number of hauls present after data quality control #
Myctophidae.cropped <- merge(Myctophidae_species, data, by = "eventID", all.x = F, all.y = F)

abundance_after <- Myctophidae.cropped %>%
  group_by(scientificName) %>%
  summarise(haul.after = n(),
            total.after = sum(individualCount),
            lat = paste0(abs(max(lat)), "-", abs(min(lat)), "°S"))
abundance_after

abundance_after <- na.omit(abundance_after)
abundance_after


# Export table 
final_abundance <- merge(abundance_before, abundance_after, by = "scientificName")
final_abundance <- final_abundance %>%
  arrange(desc(haul.after))
final_abundance

write.csv(final_abundance, "Abundance.csv", row.names = F)
