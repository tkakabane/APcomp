######   ORGANIZE NEOTOMA DATA - FOR TERRESTRIAL POLLEN  #######
# Based on the workflow presented in https://open.neotomadb.org/EPD_binder/simple_workflow.html


library(neotoma2)


####    Main directory
# DEFINE DIRECTORY

mainDir <- "C:/"

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


########################################
######## ORGANIZE NEOTOMA DATA #########
########################################


PollenDataL = list.files(file.path(mainDir, RegionName, paste0("Neotoma_Raw_Data ", RegionName)), full.names = T)

# Simplify and homogenize taxa names ########
#clean_taxa_names <- function(names, unique_only = TRUE) {
  #Remove ALL descriptors
  #names <- gsub("[?]", "", names)  # Remove "?"
  #names <- gsub("\\s*\\([^)]+\\)", "", names)  # Remove parentheses and their content
  #names <- gsub("\\b(aff\\.|cf\\.|sp\\.|undiff\\.?)\\s*", "", names, ignore.case = TRUE)  # Remove qualifiers
  #names <- gsub("-type\\b(I{0,3}|\\s.*)?$", "", names, ignore.case = TRUE)  # Remove ALL -type suffixes (including Roman numerals)
  #names <- trimws(names)  # Trim whitespace
  
  # Define replacement patterns and corresponding replacements
  #replacements <- list(
  #  "-type$" = "",
  #  "Alchornea/Conceveibum|Aparisthmium|Conceveibum|Conceveiba|Alchorneopsis" = "Alchornea",
  #  "Chenopodium/Amaranthus|Chenopodium|Arenaria" = "Amaranthus",
  #  "Umbelliferae" = "Apiaceae",
  #  "Asclepiadeae" = "Apocynaceae",
  #  "\\b(Ambrosia|Artemisia|Gnaphalium|Baccharis|Bidens|Eupatorium|Carduoideae|Cichorioideae|Helianthus|Mikania|Mutisioideae|Senecio|Vernonia|Xanthium|Aspilia|Baccharis-type|Asteroideae|Asteraceae-type|Ambrosia-type|Aster)\\b" = "Asteraceae", #Mutisia - arboreal 
  #  "Borreria-type|Mitracarpus|Spermacoce" = "Borreria",
  #  "Crotonoideae" = "Croton",
  #  "Carex|Cyperus|Eleocharis|Fimbristylis|Fuirena|Rhynchospora|Scleria|Cyperaceae-type" = "Cyperaceae",
  #  "Euterpe/Geonoma-type" = "Euterpe",
  #  "Eucryphia/Caldcluvia|Caldcluvia|Caldcluvia/Eucryphia" = "Eucryphia",
  #  "Fabales" = "Faboideae",
  #  "Gomphrena-type|Gomphrena/Pfaffia" = "Gomphrena",
  #  "Mauritia/Mauritiella|Mauritiella|Mauritia carana|Mauritia flexuosa" = "Mauritia",
  #  "Combretum|Combretaceae/Melastomataceae|Melastomataceae/Combretaceae|Combretaceae" = "Melastomataceae",
  #  "mimosoid" = "Mimosoideae",
  #  "Urticaceae/Moraceae-type|Urticaceae/Moraceae|Urticalean|Moraceae|Brosimum|Castilla|Ficus|Maclura|Morus|Sorocea|Trophis" = "Moraceae",
  #  "Eucalyptus|Eugenia|Myrcia|Myrceugenia|Myrciaria|Psidium|Campomanesia|Melaleuca|Eugenia/Melaleuca|Luma|Myrteola|Ugni" = "Myrtaceae",
  #  "Rapanea|Myrsinoideae" = "Myrsine",
  #  "\\b(Alopecurus|Avena|Festuca|Paspalum|Phleum|Poa|Secale)\\b" = "Poaceae",
  #  "Didymopanax-type|Didymopanax" = "Schefflera",
  #  "Tabebuia|Handroanthus" = "Tabebuia/Handroanthus",
  #  "?Lauraceae" =  "Lauraceae"
  #)
  
  
  # Apply taxonomic group replacements (after general cleanup)
  #for (pattern in names(replacements)) {
  #  names <- gsub(pattern, replacements[[pattern]], names)
  #}
  
  # Final standardization
  #names <- sapply(strsplit(names, "/"), function(x) x[1])  # Keep first part before "/"
  #names <- sapply(strsplit(names, " "), function(x) x[1])  # Keep first word
  
  # Return result
  #if (unique_only) {
  #  return(unique(names))
  #} else {
  #  return(names)
  #}
#} # Not used 
#combined_data <- do.call(rbind, lapply(PollenDataL, read.csv))
#first_namestax <- clean_taxa_names(combined_data$variablename, unique_only = F)
#combined_namstax <- cbind(combined_data$variablename,first_namestax,combined_data$ecologicalgroup)
#write.csv(combined_namstax,paste0(mainDir,"All_tables_mergedNeotoma.csv"))
#################


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
    
    #AP <- AP[rownames(AP) != "Mauritia", ]
    #AP <- AP[rownames(AP) != "Cyperaceae", ]
    
    
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
      MANGRad <- MANGRad[!(sample.sum < 100), ]
      
    }

    #length(mergedfPerc$UPHE)
    categories_depths <- unique(PollenData$depth)
    UPHERad <- rbind(categories_depths, categories_ages, mergedfPerc$UPHE)
    UPHERad <- t(UPHERad)
    colnames(UPHERad) <- c("Depth", "Age", "UPHE")
    UPHERad <- UPHERad[!(sample.sum < 100), ]
    
    if ("PALM" %in% pollen_m) {
      TRSHRad <- rbind(categories_depths,
                       categories_ages,
                       mergedfPerc$TRSH + mergedfPerc$PALM)
      TRSHRad <- t(TRSHRad)
      colnames(TRSHRad) <- c("Depth", "Age", "TRSH")
      TRSHRad <- TRSHRad[!(sample.sum < 100), ]
      
      
    } else {
      TRSHRad <- rbind(categories_depths, categories_ages, mergedfPerc$TRSH)
      TRSHRad <- t(TRSHRad)
      colnames(TRSHRad) <- c("Depth", "Age", "TRSH")
      TRSHRad <- TRSHRad[!(sample.sum < 100), ]
    }
    
    
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " TRSH")), showWarnings = F)
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " UPHE")), showWarnings = F)
    
    if (length(sample.sum[sample.sum > 100]) > 4) {
      #REMOVE SITES WITH ONLY ONE REMAINING SAMPLE
  
  if("MANG" %in% PollenData$ecologicalgroup == T){
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " MANG_TRSH")), showWarnings = F)
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " MANG_UPHE")), showWarnings = F)
    dir.create(file.path(mainDir, RegionName, paste0(RegionName, " MANG")), showWarnings = F)
    
    write.csv(UPHERad, paste0(mainDir,RegionName, "/",paste0(RegionName, " MANG_TRSH"),"/MANG_TRSH ",file_name,".csv"), row.names = F)
    write.csv(UPHERad, paste0(mainDir,RegionName, "/",paste0(RegionName, " MANG_UPHE"),"/MANG_UPHE ",file_name,".csv"), row.names = F)
    write.csv(MANGRad, paste0(mainDir,RegionName, "/",paste0(RegionName, " MANG"),"/MANG ",file_name,".csv"), row.names = F)
    
    if (max(MANGRad[, 3]) < 15) {
      
      write.csv(TRSHRad, paste0(mainDir,RegionName, "/",paste0(RegionName, " TRSH"),"/MANG_TRSH ",file_name,".csv"), row.names = F)
      write.csv(UPHERad, paste0(mainDir,RegionName, "/",paste0(RegionName, " UPHE"),"/MANG_UPHE ",file_name,".csv"), row.names = F)
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


#######

############################
############################
library(paleofire)

mainDir <-  "C:/"   #Ideally same directory as used in Get_NeotomaData.R

proxy = "TRSH" # CHAR OR TRSH
nameregion = "Neotropics"


#Set parameters for generating the composite curve
TarAge = 20 #Pre-binning window value
binhw = 40 #Bin half window value
hw = 1000 #Half window value
basePeriod = c(200, 21000)



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

write.csv(comp_data,paste0(mainDir,nameregion,"/",nameregion,"_",proxy,"_1000_40_20_PRSBP.csv"))











library(dplyr)

nameregion <- "Neotropics"

# Folder path containing CSV files
folder_path <- paste0("C:/.../Data used/", nameregion, "/", nameregion, " CHAR/")

# List CSV files
files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Function to read and tag each CSV file by its name (without .csv)
read_and_tag <- function(filepath) {
  filename <- basename(filepath)
  file_id <- sub("\\.csv$", "", filename)  # Remove the ".csv" extension
  
  df <- read.csv(filepath)
  df <- df %>% mutate(file_id = file_id)  # Add filename (no extension) as a column
  
  return(df)
}

# Read and combine all CSV files
all_data <- bind_rows(lapply(files, read_and_tag)) %>%
  dplyr::select(file_id, everything())

write.csv(all_data, paste0(mainDir, nameregion, "All_DATA_CHAR.csv"))
print(all_data)

