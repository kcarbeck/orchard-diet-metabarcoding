##Script modified from Eliot Miller. Attempt to streamline how we validate species lists from metabarcoding data

library(rgbif)
library(ggplot2)
library(data.table)
library(dplyr)

setwd("/Users/jlw395/Desktop")
# the download code is hashed out below so you don't accidentally run it again.
# it took a data frame confusingly titled "species", with a column within it labeled,
# "species". basically, I just passed it a vector of species scientific names.

#download GBIF occurrences for all these species
usr <- "jenwalsh123" 
pwd <- "rucwor-hijSov-3xizmo"
email <- "jlw395@cornell.edu" 

#download GBIF occurrences for all these species
#Match species names to taxon keys from GBIF
gbif_taxon_keys <- 
Species %>% 
pull("species") %>% # use fewer names if you want to just test
name_backbone_checklist()  %>% # match to backbone
filter(matchType %in% c("EXACT","FUZZY")) %>% # get matched names
pull(usageKey) # get the gbif taxonkeys


#occ_download
occ_download(pred_in("taxonKey", gbif_taxon_keys), format = "SIMPLE_CSV", user=usr,pwd=pwd,email=email)

occ_download_wait('0009877-250415084134356')

system.time(dat <- as.data.frame(fread("0009877-250415084134356.csv")))

#cull any without coordinates or coordinateUncertaintyInMeters > 20,000
dat2 <- dat[complete.cases(dat[,c("decimalLatitude","decimalLongitude")]),]

#make the risky assumption that coordinate uncertainty NA = 0
dat2$coordinateUncertaintyInMeters[is.na(dat2$coordinateUncertaintyInMeters)] <- 0
dat2 <- dat2[dat2$coordinateUncertaintyInMeters <= 20000,]

#cut any coordinates outside of the US
dat2 <- dat2[dat2$decimalLongitude >= -170 & dat2$decimalLongitude <= -50,]
dat2 <- dat2[dat2$decimalLatitude >= 15 & dat2$decimalLatitude <= 80,]

test <- subset(dat2, species=="Trachelipus rathkii")

pdf("gbif_22April2025T.NYS.pdf")
ggplot(NYS, aes(decimalLongitude, decimalLatitude)) +
  borders(database="world", regions=c("canada","mexico","usa"), fill="gray70",
          xlim = c(-170, -50), ylim = c(15, 80)) +
  geom_point(alpha=0.02, size=0.9, pch=16) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_quickmap() +
  theme(legend.position="none")
dev.off()

US <- subset(dat2, countryCode == "US")
NYS <- dat2[grep("New York", dat2$stateProvince), ]

write.table(NYS, "NYS.txt", sep=",")

##Trying to just look at NYS
#cut any coordinates outside of the US
dat3 <- dat2[dat2$decimalLongitude >= -71.51 & dat2$decimalLongitude <= -79.77,]
dat3 <- dat2[dat2$decimalLatitude >= 40.51 & dat2$decimalLatitude <= 45.11,]



