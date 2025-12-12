######   ORGANIZE NEOTOMA DATA - FOR TERRESTRIAL POLLEN  #######
# Based on the workflow presented in https://open.neotomadb.org/EPD_binder/simple_workflow.html
library(neotoma2)
library(paleofire)

####    Main directory
# DEFINE DIRECTORY
mainDir <- "C:/"

####    Name either of an existing folder or the one that
####    will be created within the main directory to store raw Neotoma data
RegionName <- "Region_name"

#Minimum counts to include a sample (all samples with less counts will be excluded)
min.count <- 100

#####      STEP 1    ######
#### GET NEOTOMA DATA #####

####    Define coordinates of the area of interest
## Coordinates can be obtained at https://geojson.io/
cz <- list(geoJSON = '{"type": "Polygon",
        "coordinates": [[
            [-33, -58],
            [-50, -58],
            [-50, -20],
            [-33, -20],
            [-33, -58]
           ]]}')

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

cz_datasets <- cz_sites %>% get_downloads(all_data = TRUE)
allSamp <- samples(cz_datasets)
#Check Ecological Groups
unique(allSamp$ecologicalgroup)


# Create a filtered version of allSamp
allSamp <- allSamp %>%
  dplyr::filter(
    element == "pollen",                      # Keep only "pollen" elements
    ecologicalgroup %in% c(                 # Include specific ecological groups
      "UPHE", "TRSH", "PALM",  #Neotoma classification UPHE - upland herbs; TRSH - trees and shrubs; MANG - MAngrove
      "MANG", "AQVP"           #AQVP - Aquatic vascular plants; PALM - palms
    ),
    !is.na(element)                          # Remove rows with NA in 'element'
  )

#Check Remaining Ecological Groups
unique(allSamp$ecologicalgroup)

#Create folder if it doesn't exist
data_list <- split(allSamp, f = allSamp$datasetid)
dir.create(file.path(mainDir, RegionName), showWarnings = T)
dir.create(file.path(mainDir, RegionName, paste0("Neotoma_Raw_Data ", RegionName)), showWarnings = T)

write.csv(as.data.frame(cz_sites[-1]), paste0(mainDir,RegionName,"/Site_metadata_",RegionName,".csv"))

##### Save files in a folder with the region name
for (i in 1:length(data_list)) {
  
  clean_name <- gsub("[^\\p{L}0-9_ -]", "", data_list[[i]]$sitename[1], perl = TRUE)
  
  # Check if 'age' has only one unique value
  if (length(unique(data_list[[i]]$age)) > 2) {
    write.csv(data_list[[i]], 
              paste0(mainDir,RegionName,"/",paste0("Neotoma_Raw_Data ",RegionName),"/",
                     i,"_", clean_name,"_", data_list[[i]]$siteid[1], ".csv"))
  }
}

############     STEP 2      ###########
######## ORGANIZE NEOTOMA DATA #########

PollenDataL = list.files(file.path(mainDir, RegionName, paste0("Neotoma_Raw_Data ", RegionName)), full.names = T)

for (k in 1:length(PollenDataL)) {
  file_name <- basename(PollenDataL[k])
  file_name <- gsub(".csv$", "", file_name)
  PollenData = read.csv(PollenDataL[k], check.names = FALSE)
  
  PollenData = PollenData[!is.na(PollenData$age), ]
  
  depthData <- PollenData$depth
  
  radiocarb <- PollenData$age
  radiocarb <- as.numeric(radiocarb)
  
  ####Create a matrix for the sample
  categories_taxa <- unique(PollenData$variablename)
  

  if (length(categories_taxa) > 2 && length(unique(radiocarb)) > 2) {
    # REMOVE IF THERE IS ONLY ONE TAXON IN THE RECORD
    # Select only the 'age' and 'depth' columns
    categories_ages <- PollenData[, c("age", "depth")]
    # Identify the unique combinations of 'age' and 'depth'
    categories_ages <- unique(categories_ages)
    # Extract the 'age' values from the unique combinations
    categories_ages <- categories_ages$age
    
    pollen_m <- matrix(nrow = length(categories_taxa),
                       ncol = length(categories_ages))
    colnames(pollen_m) <- categories_ages
    row.names(pollen_m) <- categories_taxa
    
    # Apply cleaning function to the first part of variablename
    #first_namestaxVar <- sapply(strsplit(PollenData$variablename, " "), function(x) x[1])
    first_namestaxVar <- PollenData$variablename
    
    # Loop through each cell in the matrix
    for (i in 1:nrow(pollen_m)) {
      for (j in 1:ncol(pollen_m)) {
        # Check if the row name and column name match a row in the dataframe
        match_row <- first_namestaxVar == row.names(pollen_m)[i] &
          PollenData$age == colnames(pollen_m)[j]
        
        if (any(match_row, na.rm = TRUE)) {
          # If there is a match, fill in the value from the dataframe
          pollen_m[i, j] <- sum(PollenData$value[match_row])
        } else {
          # If there is no match, leave the cell blank
          pollen_m[i, j] <- 0
        }
      }
    }
    
    pollen_m <- cbind(NA, pollen_m)
    colnames(pollen_m) <- c("group", categories_ages)
  
    # Fill in the first column with the ecological groups
    for (i in 1:nrow(pollen_m)) {
      # Find the matching row in the dataframe
      match_row <- first_namestaxVar == row.names(pollen_m)[i]
      
      if (any(match_row, na.rm = TRUE)) {
        # If there is a match, fill in the corresponding ecological group
        pollen_m[i, 1] <- PollenData$ecologicalgroup[match_row][1]
      } else {
        # If there is no match, leave the cell blank
        pollen_m[i, 1] <- NA
      }
    }
    
    AP <- as.data.frame(pollen_m)
    AP <- AP[!(AP$group %in% c("MANG", "AQVP", "SEED")), ]
    
    #AP <- AP[rownames(AP) != "Mauritia", ] # If there are taxa to exclude
    #AP <- AP[rownames(AP) != "Cyperaceae", ] # If there are taxa to exclude

    # Transpose and adjust column names
    AP <- t(AP)
    colnames(AP) <- AP[1, ]
    AP <- as.data.frame(AP[-1, ])
    AP <- as.data.frame(sapply(AP, as.numeric)) # Convert all values to numeric
    AP <- rowsum(t(AP), group = colnames(AP), na.rm = T)
    AP <- t(AP)
      
    sample.sum <- rowSums(AP, na.rm = T)  # Calculate row sums to identify rows with sufficient data
    
    #AP <- AP[which(sample.sum > 100),]  # Subset rows where the sum exceeds 100
    mergedfPerc <- (AP * 100) / sample.sum #; rowSums(mergedfPerc, na.rm = T)  # Calculate relative percentages (normalization by row sums)
    mergedfPerc <- as.data.frame(mergedfPerc)
    
    
    if ("MANG" %in% pollen_m) {
      MANG <- as.data.frame(pollen_m)
      MANG <- as.data.frame(MANG)
      MANG <- t(MANG)
      colnames(MANG) <- MANG[1, ]
      MANG <- as.data.frame(MANG[-1, ])
      MANG <- as.data.frame(sapply(MANG, as.numeric))  # Convert all values to numeric
      MANG <- rowsum(t(MANG), group = colnames(MANG), na.rm = T)
      MANG <- t(MANG)
      
      #sample.sum <- rowSums(MANG, na.rm = T) # Calculate row sums to identify rows with sufficient data
      mergedfPercMANG <- (MANG * 100) / sample.sum
      rowSums(mergedfPercMANG, na.rm = T)  # Subset rows where the sum exceeds 100
      mergedfPercMANG <- as.data.frame(mergedfPercMANG)
      
      categories_depths <- unique(PollenData$depth)
      MANGRad <- rbind(categories_depths, categories_ages, mergedfPercMANG$MANG)
      MANGRad <- t(MANGRad)
      colnames(MANGRad) <- c("Depth", "Age", "MANG")
      MANGRad <- MANGRad[!(sample.sum < min.count), ]
      
    }

    #length(mergedfPerc$UPHE)
    categories_depths <- unique(PollenData$depth)
    UPHERad <- rbind(categories_depths, categories_ages, mergedfPerc$UPHE)
    UPHERad <- t(UPHERad)
    colnames(UPHERad) <- c("Depth", "Age", "UPHE")
    UPHERad <- UPHERad[!(sample.sum < min.count), ]
    
    if ("PALM" %in% pollen_m) {
      TRSHRad <- rbind(categories_depths,
                       categories_ages,
                       mergedfPerc$TRSH + mergedfPerc$PALM)
      TRSHRad <- t(TRSHRad)
      colnames(TRSHRad) <- c("Depth", "Age", "TRSH")
      TRSHRad <- TRSHRad[!(sample.sum < min.count), ]
      
      
    } else {
      TRSHRad <- rbind(categories_depths, categories_ages, mergedfPerc$TRSH)
      TRSHRad <- t(TRSHRad)
      colnames(TRSHRad) <- c("Depth", "Age", "TRSH")
      TRSHRad <- TRSHRad[!(sample.sum < min.count), ]
    }
    
    
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " TRSH")), showWarnings = F)
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " UPHE")), showWarnings = F)
    
    if (length(sample.sum[sample.sum > min.count]) > 4) {
      #REMOVE SITES WITH ONLY ONE REMAINING SAMPLE
  
  if("MANG" %in% PollenData$ecologicalgroup == T){
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " MANG_TRSH")), showWarnings = F)
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " MANG_UPHE")), showWarnings = F)
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " MANG")), showWarnings = F)
    
    write.csv(UPHERad, paste0(mainDir,RegionName, "/",paste0(RegionName, " MANG_TRSH"),"/MANG_TRSH ",file_name,".csv"), row.names = F)
    write.csv(UPHERad, paste0(mainDir,RegionName, "/",paste0(RegionName, " MANG_UPHE"),"/MANG_UPHE ",file_name,".csv"), row.names = F)
    write.csv(MANGRad, paste0(mainDir,RegionName, "/",paste0(RegionName, " MANG"),"/MANG ",file_name,".csv"), row.names = F)
    
    if (max(MANGRad[, 3]) < 15) { # This number indicates the minimum contribution of mangrove taxa to consider the sample
      
      write.csv(TRSHRad, paste0(mainDir,RegionName, "/",paste0(RegionName, " TRSH"),"/TRSH_MANG ",file_name,".csv"), row.names = F)
      write.csv(UPHERad, paste0(mainDir,RegionName, "/",paste0(RegionName, " UPHE"),"/UPHE_MANG ",file_name,".csv"), row.names = F)
    }} else {
      
    write.csv(UPHERad, paste0(mainDir,RegionName, "/",paste0(RegionName, " UPHE"),"/UPHE ",file_name,".csv"), row.names = F)
    write.csv(TRSHRad, paste0(mainDir,RegionName, "/",paste0(RegionName, " TRSH"),"/TRSH ",file_name,".csv"), row.names = F)} 
    } else {
      
      print(paste0(file_name, " - only one sample finding the criteria"))
    }
  } else {
    print(paste0(file_name, " - Not analyzed: Contains only one taxon or a radiocarbon age"))
  }
}

#############      STEP 3        #############
####### PALEOFIRE PACKAGE STARTS HERE ########
#library(paleofire)
#mainDir <-  "C:/"   #Ideally same directory as used in Get_NeotomaData.R

proxy = "TRSH" # CHAR OR TRSH
nameregion = RegionName # NAME OF THE REGION

#Set parameters for generating the composite curve
TarAge = 20 #Pre-binning window value
binhw = 40 #Bin half window value
hw = 1000 #Half window value
basePeriod = c(200, 21000) #Define base period for data treatment



TaddData = list.files(paste0(mainDir,nameregion,"/",nameregion," ",proxy,"/"), full.names = T)

TaddData <- TaddData[basename(TaddData) != "desktop.ini"]
file_name <- basename(TaddData)
file_name <- file_name[file_name != "desktop.ini"]
file_name <- gsub(".csv$", "", file_name)

Tmydata <- pfAddData(
  files = TaddData, 
  type = "NULL")


trans_data <- pfTransform(add = Tmydata)
trans_data1 <- trans_data

trans_data1 <- pfTransform(trans_data1, 
                           QuantType = "NONE",
                           BasePeriod = c(200,21000),
                           method=c("MinMax")
)

trans_data1 <- pfTransform(trans_data1,
                           QuantType = "NONE",
                           alpha = 0.01,
                           method = c("Box-Cox")
)

trans_data1 <- pfTransform(trans_data1, 
                           QuantType = "NONE",
                           method=c("MinMax")
)

trans_data1 <- pfTransform(trans_data1, 
                           QuantType = "NONE",
                           BasePeriod = c(200,21000), 
                           method=c("Z-Score")
)

Data_comp <- pfCompositeLF(
  trans_data1,
  binhw = binhw, #pre-binning
  #pre-binning
  hw = hw,
  nboot = 1000,
  conf = c(0.025, 0.975),
  pseudodata = F,
  tarAge = seq(0, 21000, TarAge) 
)


#
par(mar= c(2,2,1,1))
plot(Data_comp
     ,add = "sitenum"
     #, ylim = c(-1.5,1.5)
)

N_records <- rowSums(!is.na(Data_comp$BinnedData))
comp_data <- cbind(Data_comp$Result, N_records)

write.csv(comp_data,paste0(mainDir,nameregion,"/",nameregion,"_",proxy,"_hw",hw,"_binhw",binhw,"_tarAge",TarAge,".csv"))

