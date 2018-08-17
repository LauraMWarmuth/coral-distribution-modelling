# MaxEnt modelling for all Scleractinia species, per region, per species

library(here)
library(rgbif)
# library(maxnet)
library(sp)
library(raster)
library(maptools)
library(rgdal)
library(dismo)
library(SDMTools)
library(SSDM)
library(sdmpredictors)
library(leaflet)

# gc() # for checking memory space

# Exploring the marine datasets
datasets <- list_datasets(terrestrial = FALSE, marine = TRUE)

# Exploring the marine layers 
layers <- list_layers(datasets)
here()
# Download specific layers to the current directory
# distance to shore maybe as well?

# BO2_salinitymean_ss Sea surface salinity (mean) instead normal salinity?

env.variables <- load_layers(c("BO_sstmean", "BO_salinity", "BO_dissox", "BO_bathymean"), datadir = "/Bio-ORACLE")
BO_bathymean <- load_layers("BO_bathymean", datadir = "/Bio-ORACLE")
names(BO_bathymean)
plot(env.variables)

###

# SHAPEFILE FOR DEPTH FROM BIO ORACLE OR NOAA

# Get bathymetry data from NOAA tropics
library(marmap)
NOAAbathy.tropics <- getNOAA.bathy(lon1 = -180,lon2 = 180, lat1 = -23.5, lat2 = 23.5, resolution = 4, keep = TRUE)

# Get bathymetry data from NOAA coral triangle
NOAAbathy.coral.triangle <- getNOAA.bathy(lon1 = 11,lon2 = 171, lat1 = -15, lat2 = 17, resolution = 4, keep = TRUE) #antimeridian = TRUE)
# longitude 11 not coral triangle???
NOAAbathy.coral.triangle.new <- getNOAA.bathy(lon1 = 95,lon2 = 171, lat1 = -15, lat2 = 17, resolution = 4, keep = TRUE, antimeridian = FALSE)


NOAAbathy.coral.triangle.new.df <- as.xyz(NOAAbathy.coral.triangle.new)
attach(NOAAbathy.coral.triangle.new.df)

# Subset shallow reefs
NOAAbathy.coral.triangle.shallow.new <- subset(NOAAbathy.coral.triangle.new.df, V3 == 0:30)
dim(NOAAbathy.coral.triangle.shallow.new)
head(NOAAbathy.coral.triangle.new.shallow)
colnames(NOAAbathy.coral.triangle.shallow.new) <- c("decimalLongitude","decimalLatitude", "Depth")
write.csv(NOAAbathy.coral.triangle.shallow.new, file = "NOAAbathy.coral.triangle.shallow30m.new.csv")
head(NOAAbathy.coral.triangle.shallow.new)
# Create shapefile from csv
library(maptools)
library(rgdal)
library(sp)
WGScoor <-  NOAAbathy.coral.triangle.shallow.new
coordinates(WGScoor)=~decimalLongitude+decimalLatitude
proj4string(WGScoor)<- CRS("+proj=longlat +datum=WGS84")
# could already save as shape file
# but project out of WGS84 (why?)
LLcoor <- spTransform(WGScoor, CRS("+proj=longlat"))
raster::shapefile(LLcoor, "NOAAbathy.coral.triangle.shallow30m.new.shp")

# Read shapefile.shp from current wd
library(rgdal)
NOAAbathy.coral.triangle.shape <- readOGR(dsn = ".", layer = "NOAAbathy.coral.triangle.shallow30m.new")
head(NOAAbathy.coral.triangle.shape)

env.variables.coral.triangle.southeast.asia.masked <- mask(env.variables.coral.triangle.southeast.asia, NOAAbathy.coral.triangle.shape)


# Exploring the available future/paleo marine layers
future_layers <- list_layers_future(terrestrial=FALSE)
# Available scenarios
future_scenario <- unique(future_layers$scenario)
future_year <- unique(future_layers$year)

# A1B: rapic economic growth 1.4 - 6.4 °C, balanced emphasis on all energy sources
# A2: regionally oriented economic development 2.0 -5.4 °C
# B1: global environmental sustainability 1.1 - 2.9 °C

env.variables.future <- load_layers(c("BO_A1B_2100_sstmean", "BO_A1B_2200_sstmean", "BO_A2_2100_sstmean", "BO_B1_2100_sstmean", "BO_B1_2200_sstmean",
                                      "BO_A1B_2100_salinity", "BO_A1B_2200_salinity", "BO_A2_2100_salinity", "BO_B1_2100_salinity", "BO_B1_2200_salinity"), datadir = "/Bio-ORACLE")

env.variables.future.salinitymean.ss <- load_layers(c("BO2_RCP26_2050_salinitymean_ss", "BO2_RCP45_2050_salinitymean_ss", "BO2_RCP60_2050_salinitymean_ss", "BO2_RCP85_2050_salinitymean_ss", "BO2_RCP26_2100_salinitymean_ss", "BO2_RCP45_2100_salinitymean_ss", "BO2_RCP60_2100_salinitymean_ss", "BO2_RCP85_2100_salinitymean_ss"), datadir = "/Bio-ORACLE")

env.variables.future.A1B <-load_layers(c("BO_A1B_2100_sstmean", "BO_A1B_2100_salinity"), datadir = "/Bio-ORACLE")
env.variables.future.B1 <- load_layers(c("BO_B1_2100_sstmean", "BO_B1_2100_salinity"), datadir = "/Bio-ORACLE")

env.variables.future.A1B.2100.sstmean <- load_layers("BO_A1B_2100_sstmean", datadir = "/Bio-ORACLE")
env.variables.future.B1.2100.sstmean <- load_layers("BO_B1_2100_sstmean", datadir = "/Bio-ORACLE")

env.variables.future.A1B.2220.sstmean <- load_layers("BO_A1B_2200_sstmean", datadir = "/Bio-ORACLE")
env.variables.future.B1.2200.sstmean <- load_layers("BO_B1_2200_sstmean", datadir = "/Bio-ORACLE")

# RGBIF Multiple species occurrences, example Scleractinia in tropics
# Tropics or warm water corals: latitude 23.5 or 30?
scleractinia.tropics.all <- occ_search(scientificName = "Scleractinia", decimalLatitude = '-23.5, 23.5',
                                       decimalLongitude = '-180, 180', limit = 200000)
scleractinia.tropics.all.data <- scleractinia.tropics.all[["data"]]
write.csv(scleractinia.tropics.all.data, file = "scleractinia.tropics.all.csv")

# Create subsets for coordinates and countries
scleractinia.tropics.all.coordinates <- scleractinia.tropics.all.data[,3:4]
write.csv(scleractinia.tropics.all.coordinates, file = "scleractinia.tropics.all.coordinates.csv")

scleractinia.tropics.all.coord.country <- scleractinia.tropics.all.data[,c(3,4,48)]
write.csv(scleractinia.tropics.all.coord.country, file = "scleractinia.tropics.all.coord.country.csv")

scleractinia.tropics.all.country <- scleractinia.tropics.all.data[,48]
write.csv(scleractinia.tropics.all.country, file = "scleractinia.tropics.all.country.csv")

# Data cleaning (SDM Dismo Vignette Hijmans, 2017)
# Remove missing coordinates
scleractinia.tropics.all.coordinates.clean <- subset(scleractinia.tropics.all.coordinates, !is.na(decimalLongitude) & !is.na(decimalLatitude))
dim(scleractinia.tropics.all.coordinates)
dim(scleractinia.tropics.all.coordinates.clean)

# Removing duplicates (95 duplicate rows across all columns)
# distinct() in dplyr to keep only unique/distinct rows from a data frame
library("tibble") # better to manage large datasets
scleractinia.tropics.all.data.df <- as_data_frame(scleractinia.tropics.all.data)
head(scleractinia.tropics.all.data.df)
library("dplyr")
scleractinia.tropics.all.data.df.unique <- distinct(scleractinia.tropics.all.data.df)

# Subsets without replicates
scleractinia.tropics.all.coordinates.unique <- scleractinia.tropics.all.data.df.unique[,3:4]
write.csv(scleractinia.tropics.all.coordinates.unique, file = "scleractinia.tropics.all.coordinates.unique.csv")

scleractinia.tropics.all.coord.country.unique <- scleractinia.tropics.all.data.df.unique[,c(3,4,48)]
write.csv(scleractinia.tropics.all.coord.country.unique, file = "scleractinia.tropics.all.coord.country.unique.csv")

scleractinia.tropics.all.country.unique <- scleractinia.tropics.all.data.df.unique[,48]
write.csv(scleractinia.tropics.all.country.unique, file = "scleractinia.tropics.all.country.unique.csv")


# Removing missing genera
attach(scleractinia.tropics.all.data.df.unique)
scleractinia.tropics.all.data.df.unique.new <- subset(scleractinia.tropics.all.data.df.unique, !is.na(genus))
dim(scleractinia.tropics.all.data.df.unique)
# [1] 197161    180
dim(scleractinia.tropics.all.data.df.unique.new)
# [1] 194092    180

# Map Scleractinia occurrences
data(wrld_simpl)
plot(wrld_simpl, 
     xlim = c(-180, 180),
     ylim = c(-23.5, 23.5), # ylim -23.5, 23.5 doesnt cut map?
     axes = TRUE, 
     col = "grey95")
points(x = scleractinia.tropics.all.coordinates.unique$decimalLongitude,
       y = scleractinia.tropics.all.coordinates.unique$decimalLatitude, 
       col = "red",
       pch = 16,
       cex = 0.3)
plot(wrld_simpl, add = TRUE, border = "grey5")

# Restore the box around the map
box()

# Remove points on land, "click" function or row numbers as labels to delete rows in dataset?
# NOT NECESSARY TO DO SINCE MARINE ENV VARIABLES ARE ONLY OCEAN BASED ANYWAYS!!!


scleractinia.tropics.all.coordinates.unique.df <- as.data.frame(scleractinia.tropics.all.coordinates.unique)
spatialpoints <- SpatialPoints(scleractinia.tropics.all.coordinates.unique, proj4string=CRS(as.character(NA)), bbox = NULL)
head(spatialpoints)
library("raster")
click(spatialpoints, id=TRUE, xy= FALSE, show = TRUE)

# Split data into subsets for different regions
unique.countries <- unique(scleractinia.tropics.all.country.unique, incomparables = FALSE)
write.csv(unique.countries, file = "unique.countries.csv")

attach(scleractinia.tropics.all.coord.country.unique)
western.pacific <- c("American Samoa", "Australia", "Brunei Darussalam", "China", "Cook Islands", "Fiji", "French Polynesia", "Guam", "Hong Kong", "Japan", "Kiribati", "Lao People's Democratic Republic", "Marshall Islands", "Micronesia, Federated State of", "New Caledonia", "New Zealand", "Niue", "Northern Mariana Islands", "Palau", "Samoa", "Singapore", "Taiwan", "Tonga", "Tuvalu", "United States Minor Outlying Islands", "Vanuatu", "Viet Nam", "Wallis and Futuna")
eastern.pacific <- c("Chile", "Colombia(decimal)", "Costa Rica(decimal)", "Ecuador", "El Salvador", "Honduras(decimal)", "Mexico(decimalLongitude)", "Nicaragua(decimal)", "Panama(decimal)", "Peru")
caribbean.gulfofmexico <- c("Anguilla", "Antigua and Barbuda", "Aruba", "Bahamas", "Barbados", "Belize", "Bonaire, Sint Eustatius and Saba", "Brazil", "Cayman Islands", "Colombia(decimal)", "Costa Rica(decimal)", "Cuba", "Curaçao", "Dominica", "Dominican Republic", "French Guiana", "Grenada", "Guadeloupe", "Guyana", "Haiti", "Honduras(decimal)", "Jamaica", "Martinique", "Mexico(decimal)", "Montserrat", "Nicaragua(decimal)", "Panama(decimal)", "Puerto Rico", "Saint Barthélemy", "Saint Kitts and Nevis", "Saint Lucia", "Saint Martin (French part)", "Saint Vincent and the Grenadines", "Sint Maarten (Dutch part)", "Suriname", "Trinidad and Tobago", "Turks and Caicos Islands", "Venezuela, Bolivarian Republic of", "Virgin Islands, British", "Virigin Islands, U.S.")
western.indian.ocean <- c("British Indian Ocean Territory", "Comoros", "Dijbouti", "Egypt", "Eritrea", "Ethiopia", "India", "Kenya", "Madagascar", "Maldives", "Mauritania", "Mauritius", "Mayotte", "Mozambique", "Oman", "Réunion", "Saudi Arabia(decimal)", "Seychelles", "Somalia", "Sudan(decimal)", "Tanzania, United Republic of", "Yemen(decimal)", "Zimbabwe")
eastern.indian.ocean <- c("Bangladesh", "Christmas Island", "Cocos (Keeling) Islands", "Myanmar", "Sri Lanka", "Thailand")
coral.triangle.southeast.asia <- c("Indonesia", "Malaysia", "Papua New Guinea", "Philippines", "Solomon Islands", "Timor-Leste") 
eastern.atlantic <- c("Angola", "Benin", "Cabo Verde", "Cameroon", "Côte d'Ivoire", "Equatorial Guinea", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Liberia", "Morocco", "Namibia", "Nigeria", "Saint Helena, Ascension and Tristan da Cunha", "Sao Tome and Principe", "Senegal", "Sierra Leone")
#red.sea ?

western.pacific.coord <- subset(scleractinia.tropics.all.coord.country.unique, country == western.pacific)
western.pacific.coord <- western.pacific.coord[,1:2]

eastern.pacific.coord <- subset(scleractinia.tropics.all.coord.country.unique, country == eastern.pacific)
eastern.pacific.coord <- eastern.pacific.coord[,1:2]

caribbean.gulfofmexico.coord <- subset(scleractinia.tropics.all.coord.country.unique, country == caribbean.gulfofmexico)
caribbean.gulfofmexico.coord <- caribbean.gulfofmexico.coord[,1:2]

western.indian.ocean.coord <- subset(scleractinia.tropics.all.coord.country.unique, country == western.indian.ocean)
western.indian.ocean.coord <- western.indian.ocean.coord[,1:2]

eastern.indian.ocean.coord <- subset(scleractinia.tropics.all.coord.country.unique, country == eastern.indian.ocean)
eastern.indian.ocean.coord <- eastern.indian.ocean.coord[,1:2]

coral.triangle.southeast.asia.coord <- subset(scleractinia.tropics.all.coord.country.unique, country == coral.triangle.southeast.asia)
coral.triangle.southeast.asia.coord <- coral.triangle.southeast.asia.coord[,1:2]

eastern.atlantic.coord <- subset(scleractinia.tropics.all.coord.country.unique, country == eastern.atlantic)
eastern.atlantic.coord <- eastern.atlantic.coord[,1:2]

# Warning message:
#In country == eastern.atlantic :
 # longer object length is not a multiple of shorter object length

# Map CR occurrences: Atlantic and Pacific?
head(coord.costa.rica)
coord.costa.rica <- subset(scleractinia.tropics.all.coord.country.unique, country == "Costa Rica")
coord.costa.rica <- subset(coord.costa.rica[,1:2])
dim(coord.costa.rica)

data(wrld_simpl)
plot(wrld_simpl, 
     xlim = c(-180, 180),
     ylim = c(-23.5, 23.5), # ylim -23.5, 23.5 doesnt cut map?
     axes = TRUE, 
     col = "grey95")
points(x = coord.costa.rica$decimalLongitude,
       y = coord.costa.rica$decimalLatitude, 
       col = "red",
       pch = 16,
       cex = 0.3)
plot(wrld_simpl, add = TRUE, border = "grey5")
# occurrences on both coasts - how to split?

# Crop map for coral triangle region/ apply model only for this region!!

## Crop rasters from script Occurrence_data.R
# (xmin, xmax, ymin, ymax) from 
#xmin=min(x.ocean)-5
#xmax=max(x.ocean)+5
#ymin=min(y.ocean)-5
#ymax=max(y.ocean)+5
#limits <- c(xmin, xmax, ymin, ymax)
#rasters.crop <- crop(rasters,limits)

# Determine geographic extent of coral triangle occurrences
max.lat <- ceiling(max(coral.triangle.southeast.asia.coord.df.new$decimalLatitude))
min.lat <- floor(min(coral.triangle.southeast.asia.coord.df.new$decimalLatitude))
max.lon <- ceiling(max(coral.triangle.southeast.asia.coord.df.new$decimalLongitude))
min.lon <- floor(min(coral.triangle.southeast.asia.coord.df.new$decimalLongitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))
max.lon
# Crop bio oracle data to geographic extent of coral triangle occurrences
env.variables.coral.triangle.southeast.asia <- crop(x = env.variables, y = geographic.extent)
env.variables.coral.triangle.southeast.asia.future.A1B <- crop(x = env.variables.future.A1B, y = geographic.extent)
env.variables.coral.triangle.southeast.asia.future.B1 <- crop(x = env.variables.future.B1, y = geographic.extent)

# Maxent
library(rJava)
maxent()
coral.triangle.southeast.asia.coord.df <- as.data.frame(coral.triangle.southeast.asia.coord)
coral.triangle.southeast.asia.coord.df.new <- coral.triangle.southeast.asia.coord.df[c("decimalLongitude", "decimalLatitude")]
coral.triangle.southeast.asia.coord.df.new <- coral.triangle.southeast.asia.coord.df.new[-373,]
class(coral.triangle.southeast.asia.coord.df.new)
class(env.variables.coral.triangle.southeast.asia.masked)
xm_masked30m <- maxent(env.variables.coral.triangle.southeast.asia.masked, coral.triangle.southeast.asia.coord.df.new, path="C:/Users/Laura/Desktop/MSc Dissertation/output") #nbg = 200000?
# Error in .local(x, p, ...) : 
#  more than half of the presence points have NA predictor values


# NAs get dropped or influence model??
# then for current xm incl bathymetry variable
# for future: take ss temp and ss salinity (check scenarios) -> should be automatically ss??
# alternatively cut model at certain distance from coast

# DECIDED INCLUDE BATHY AS VARIABLE

xm
plot(xm)

# To validate the model, we use the appropriately named function.
e <- evaluate(coral.triangle.southeast.asia.coord.df.new, xm, env.variables.coral.triangle.southeast.asia)
#Error in (function (classes, fdef, mtable)  : 
 #           unable to find an inherited method for function ‘evaluate’ for signature ‘"data.frame"’
pred_xm <- predict(xm_new.coord, env.variables.coral.triangle.southeast.asia) # generate the predictions
pred_xm
png(filename = "xm_new.coord.png", width = 480, height = 480)
plot(pred_xm)
dev.off()
# Future Scenario A1B balanced economic growth
futureA1B_xm <- maxent(env.variables.coral.triangle.southeast.asia.future.A1B, coral.triangle.southeast.asia.coord.df.new, path="C:/Users/Laura/Desktop/MSc Dissertation/output") #nbg = 200000
#Warning message:
# In .local(x, p, ...) :
#120 (14.72%) of the presence points have NA predictor values
futureA1B_xm
plot(futureA1B_xm)
# To validate the model, we use the appropriately named function.

e <- evaluate(coral.triangle.southeast.asia.coord.df.new, futureA1B_xm, env.variables.coral.triangle.southeast.asia.future.A1B)
#Error in (function (classes, fdef, mtable)  : 
#           unable to find an inherited method for function ‘evaluate’ for signature ‘"data.frame"’
pred_futureA1B_xm <- predict(futureA1B_xm, env.variables.coral.triangle.southeast.asia.future.A1B) # generate the predictions
pred_futureA1B_xm
png(filename = "pred_futureA1B_xm_new.coord.png", width = 480, height = 480)
plot(pred_futureA1B_xm)
dev.off()
# Future Scenario B1 global sustainability
futureB1_xm <- maxent(env.variables.coral.triangle.southeast.asia.future.B1, coral.triangle.southeast.asia.coord.df.new, path="C:/Users/Laura/Desktop/MSc Dissertation/output") #nbg = 200000, path = here
#Warning message:
# In .local(x, p, ...) :
#120 (14.72%) of the presence points have NA predictor values
futureB1_xm
plot(futureB1_xm)
# To validate the model, we use the appropriately named function.

e <- evaluate(coral.triangle.southeast.asia.coord.df.new, futureB1_xm, env.variables.coral.triangle.southeast.asia.future.B1)
#Error in (function (classes, fdef, mtable)  : 
#           unable to find an inherited method for function ‘evaluate’ for signature ‘"data.frame"’
pred_futureB1_xm <- predict(futureB1_xm, env.variables.coral.triangle.southeast.asia.future.B1) # generate the predictions
pred_futureB1_xm
png(filename = "pred_futureB1_xm_new.coord.png", width = 480, height = 480)
plot(pred_futureB1_xm)
dev.off()
# Maxent model coral triangle occurrences for each species as a loop
# Maybe genus makes more sense, since not enough presence points?
# evaluate and predict function within loop?

# Subset including species name for coral triangle
scleractinia.tropics.all.species.coord.country <- scleractinia.tropics.all.data.df.unique[,c(1,3,4,48)]
scleractinia.tropics.all.species.coord.country.df <- as.data.frame(scleractinia.tropics.all.species.coord.country)

coral.triangle.southeast.asia.species.coord <- subset(scleractinia.tropics.all.species.coord.country.df, country == coral.triangle.southeast.asia)
coral.triangle.southeast.asia.species.coord <- coral.triangle.southeast.asia.species.coord[,1:3]
coral.triangle.southeast.asia.species.coord <- coral.triangle.southeast.asia.species.coord[c("name", "decimalLongitude", "decimalLatitude")]
write.csv(coral.triangle.southeast.asia.species.coord, file = "coral.triangle.southeast.asia.species.coord.csv")

# Subset including genus name for coral triangle (new dataset without missing genera)
scleractinia.tropics.all.genera.coord.country <- scleractinia.tropics.all.data.df.unique.new[,c(3,4,30,48)]
scleractinia.tropics.all.genera.coord.country.df <- as.data.frame(scleractinia.tropics.all.genera.coord.country)

coral.triangle.southeast.asia.genera.coord <- subset(scleractinia.tropics.all.genera.coord.country.df, country == coral.triangle.southeast.asia)
coral.triangle.southeast.asia.genera.coord <- coral.triangle.southeast.asia.genera.coord[,1:3]
coral.triangle.southeast.asia.genera.coord <- coral.triangle.southeast.asia.genera.coord[c("genus", "decimalLongitude", "decimalLatitude")]
dim(coral.triangle.southeast.asia.genera.coord)
# Get genera list for coral triangle
coral.triangle.southeast.asia.genera <- coral.triangle.southeast.asia.genera.coord[,1]

unique.genera <- unique(coral.triangle.southeast.asia.genera, incomparables = FALSE)
write.csv(unique.genera, file = "unique.genera.csv")
# 173 genera
# many genera <10 observations???

# Split whole coral triangle dataset into list with object per genus
attach(coral.triangle.southeast.asia.genera.coord)
genuslist <- split(coral.triangle.southeast.asia.genera.coord, genus)
attach(genuslist)

genuslist.only.coord <- lapply(genuslist, function(x) x[(names(x) %in% c("decimalLongitude", "decimalLatitude"))] )
attach(genuslist.only.coord)

# Skip errors, put trycatch inside of loop
for (i in 1:length(genuslist.only.coord)){
  genus.df <- as.data.frame(genuslist.only.coord[[i]])
  tryCatch({mx_output <- maxent(env.variables.coral.triangle.southeast.asia, genus.df, path=paste("C:/Users/Laura/Desktop/MSc Dissertation/output/", i, sep=""))
            png(filename = paste(i, "mx_output.png", sep = "_"), width = 480, height = 480)
            plot(mx_output)
            dev.off()
            mx_predict <- predict(mx_output, env.variables.coral.triangle.southeast.asia)
            png(filename = paste(i, "mx_predict.png", sep = "_"), width = 480, height = 480)
            plot(mx_predict)
            dev.off()}, # first function
            error=function(err){print(err); print(genus.df)}) # print error
}

# RESULT: 123 genera, rest lost cause >50% NA predictor values

# Layer prediction maps for final maps with sp richness

# Trial with one file to get binary vector
species_samplePredictions.new <- species_samplePredictions[,-c(3:5)]
head(species_samplePredictions.new)
attach(species_samplePredictions.new)
# Make predictions values numeric, >0.8
species_samplePredictions.new$Cloglog.prediction.binary <- as.numeric(species_samplePredictions.new$Cloglog.prediction > 0.8)
head(species_samplePredictions.new)

# Loop over folders to get binary vectors
setwd("C:/Users/Laura/Desktop/MSc Dissertation/output")
prediction.files <- list.files("C:/Users/Laura/Desktop/MSc Dissertation/output", recursive = TRUE, pattern = "species_samplePredictions.csv")
# 124 files (actually just 123 folders??)
head(prediction.files)
# Clever way to get rid of certain files in various folders
junk <- dir(path="C:/Users/Laura/Desktop/MSc Dissertation/output", recursive = TRUE, pattern = "binary")
?dir
file.remove(junk)
?file.remove

?dir
binary.files.together <- dir(path = "C:/Users/Laura/Desktop/MSc Dissertation/output", recursive = TRUE, pattern = "binary")
file.copy(from=binary.files.together, to=new.folder, overwrite = FALSE, recursive = TRUE, copy.mode = TRUE)
OutPath <- "C:/Users/Laura/Desktop/MSc Dissertation/output/binary"

for (i in 1:lengthprediction.files) {
  # read data
  species_samplePredictions <- read.csv(fileName, 
                          header=TRUE,
                          sep= ",")
  species_samplePredictions <- species_samplePredictions[,-c(3:5)]
  species_samplePredictions$Cloglog.prediction.binary <- as.numeric(species_samplePredictions$Cloglog.prediction > 0.8)
  write.csv(species_samplePredictions, file = "C:/Users/Laura/Desktop/MSc Dissertation/output/binary.csv", sep = "")
}

binary.files <- list.files(current.folder, recursive = TRUE, pattern = "binary.csv")
current.folder <- "C:/Users/Laura/Desktop/MSc Dissertation/output"
new.folder <- "C:/Users/Laura/Desktop/MSc Dissertation/binary"
file.copy(from=binary.files, to=new.folder, overwrite = FALSE, recursive=TRUE, copy.mode=TRUE)

for (fileName2 in binary.files) {
  # read data
  species_samplePredictions.binary <- read.csv(fileName2, 
                                        header=TRUE,
                                        sep= ",")
  }

# Function multmerg from r bloggers
multmerg = function(mypath) {
  filenames = list.files(path=mypath, full.names = TRUE)
  datalist = lapply(filenames, function(x)
    {read.csv(file=x, header=T)})
    Reduce(function(x,y) {merge(x,y)}, datalist)
}


reduce(merge, list(binary.files))

for (fileName2 in binary.files){
  
  # Read original data
  binary <- read.csv(fileName2,
                     header = TRUE,
                     sep = ",")
  
  # Create new data based on contents of original file
  all.binary <- data.frame(
    File = fileName2,
    Cloglog.prediction.binary = sum(binary$Cloglog.prediction.binary))
  
  # Write new data to separate file
  write.table(all.binary,
              "binary/all_binary.csv",
              append = TRUE,
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  
}

CopyFile <- function(file){
  FileListSource<- paste("C:/Users/Laura/Desktop/MSc Dissertation/output", file, sep="")
  FileListEnd <- paste("C:/Users/Laura/Desktop/MSc Dissertation/binary", file, sep="")
  file.copy(from=FileListSource, to=FileListEnd)
}

sapply(binary.files, CopyFile)
binary.files
for (file in binary.files){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}
bind <- cbind(binary.files(1:123)$Cloglog.prediction.binary)
total.binary <- merge(binary.files, by = "Cloglog.prediction.binary")
