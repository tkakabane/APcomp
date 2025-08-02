######   ORGANIZE NEOTOMA DATA - FOR TERRESTRIAL POLLEN  #######
# Based on the workflow presented in https://open.neotomadb.org/EPD_binder/simple_workflow.html


library(neotoma2)


####    Main directory
# DEFINE DIRECTORY

mainDir <- "C:/........."

####    Name either of an existing folder or the one that
####    will be created within the main directory to store raw Neotoma data
RegionName <- "Neotropics"

options(warn = 1)

#####      STEP 1    ######
#### GET NEOTOMA DATA #####

####    Define coordinates of the area of interest 
cz <- list(geoJSON = '{"type": "Polygon",
        "coordinates": [[
            [-33, -58],
            [-105, -58],
            [-105, 28],
            [-33, 28],
           [-33, -58]]]}')
#https://geojson.io/

cz$sf <- geojsonsf::geojson_sf(cz$geoJSON)[[1]]
#
#get sites with "pollen" data
cz_sites <- neotoma2::get_sites(loc = cz$geoJSON,
                                limit = 99999,
                                datasettype = "pollen")

#display data distribution
neotoma2::plotLeaflet(cz_sites) %>%
  leaflet::addPolygons(map = .,
                       data = cz$sf,
                       color = "green")


write.csv(cz_sites[-1], paste0(mainDir,"Site_metadata_",RegionName,".csv"))


cz_datasets <- cz_sites %>% get_downloads(all_data = TRUE)
allSamp <- samples(cz_datasets)

# Create a filtered version of allSamp
allSamp <- allSamp %>%
  dplyr::filter(
    element == "pollen",                      # Keep only "pollen" elements
    ecologicalgroup %in% c(                 # Include specific ecological groups
      "UPHE", "TRSH", "PALM", 
      "MANG", "AQVP"
    ),
    !is.na(element)                          # Remove rows with NA in 'element'
  )

#Check Remaining Ecological Groups
#unique(allSamp$ecologicalgroup)



#Create folder if it doesn't exist
data_list <- split(allSamp, f = allSamp$datasetid)
dir.create(file.path(mainDir, RegionName), showWarnings = T)
dir.create(file.path(mainDir, RegionName, paste0("Neotoma_Raw_Data ", RegionName)), showWarnings = T)


##### SAve files in a folder with the region name

for (i in 1:length(data_list)) {
  
  clean_name <- gsub("[^\\p{L}0-9_ -]", "", data_list[[i]]$sitename[1], perl = TRUE)
  
  # Check if 'age' has only one unique value
  if (length(unique(data_list[[i]]$age)) > 2) {
    write.csv(data_list[[i]], 
              paste0(mainDir,RegionName,"/",paste0("Neotoma_Raw_Data ",RegionName),"/",
                     i,"_", clean_name,"_", data_list[[i]]$siteid[1], ".csv"))
  }
}