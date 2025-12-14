######   ORGANIZE NEOTOMA DATA - FOR TERRESTRIAL POLLEN  #######
# Based on the workflow presented in https://open.neotomadb.org/EPD_binder/simple_workflow.html
library(neotoma2)
library(paleofire)
library(locfit)

####    Main directory
# DEFINE DIRECTORY
mainDir <- "C:/"

####    Name either of an existing folder or the one that
####    will be created within the main directory to store raw Neotoma data
RegionName <- "Neotropics"
min.count <- 100 #Minimum counts to include a sample (all samples with less counts will be excluded)

#####      STEP 1    ######
#### GET NEOTOMA DATA #####

####    Define coordinates of the area of interest
## Coordinates can be obtained at https://geojson.io/
cz <- list(geoJSON = '{"type": "Polygon",
        "coordinates": [[
            [-33, -58],
            [-105, -58],
            [-105, 28],
            [-33, 28],
           [-33, -58]]]}')

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
    !is.na(age), 
    !is.na(element)                          # Remove rows with NA in 'element'
  )

#Check Remaining Ecological Groups
unique(allSamp$ecologicalgroup)

# Simplify and clean taxa names ########
clean_taxa_names <- function(names, unique_only = TRUE) {
  # Standard cleaning operations
  names <- gsub("[?]", "", names)  # Remove "?"
  names <- gsub("\\s*\\([^)]+\\)", "", names)  # Remove parentheses and their content
  names <- gsub("\\b(aff\\.|cf\\.|sp\\.|undiff\\.?)\\s*", "", names, ignore.case = TRUE)  # Remove qualifiers
  names <- gsub("-type\\b(I{0,3}|\\s.*)?$", "", names, ignore.case = TRUE)  # Remove ALL -type suffixes (including Roman numerals)
  names <- trimws(names)  # Trim whitespace
  
  ##### IF HOMOGENIZATION IS DESIRED, THEN REMOVE THE #### AND ADAPT THE EXAMPLES BELOW
  # Define replacement patterns and corresponding replacements
  replacements <- list(
    "-type$" = ""#,
    #"Alchornea/Conceveibum|Aparisthmium|Conceveibum|Conceveiba|Alchorneopsis" = "Alchornea",
    #"Umbelliferae" = "Apiaceae",
    #"\\b(Ambrosia|Artemisia|Gnaphalium|Baccharis|Bidens|Eupatorium|Carduoideae|Cichorioideae|Helianthus|Mikania|Mutisioideae|Senecio|Vernonia|Xanthium|Aspilia|Baccharis-type|Asteroideae|Asteraceae-type|Ambrosia-type|Aster)\\b" = "Asteraceae", #Mutisia - arboreal 
    #"Carex|Cyperus|Eleocharis|Fimbristylis|Fuirena|Rhynchospora|Scleria|Cyperaceae-type" = "Cyperaceae",
    #"Combretum|Combretaceae/Melastomataceae|Melastomataceae/Combretaceae|Combretaceae" = "Melastomataceae",
    #"Urticaceae/Moraceae-type|Urticaceae/Moraceae|Urticalean|Moraceae|Brosimum|Castilla|Ficus|Maclura|Morus|Sorocea|Trophis" = "Moraceae",
  )
  
  # Apply taxonomic group replacements (after general cleanup)
  for (pattern in names(replacements)) {
    names <- gsub(pattern, replacements[[pattern]], names)
  }
  # Keep first part before "/"
  names <- sapply(strsplit(names, "/"), function(x) x[1])
  
  # Keep first word only
  names <- sapply(strsplit(names, " "), function(x) x[1])
  
  # Return result
  if (unique_only) {
    return(unique(names))
  } else {
    return(names)
  }
}

#Create folder if it doesn't exist
data_list <- split(allSamp, f = allSamp$datasetid)
dir.create(file.path(mainDir, RegionName), showWarnings = T)
dir.create(file.path(mainDir, RegionName, paste0("Neotoma_Raw_Data ", RegionName)), showWarnings = T)
cz_sitedf <- as.data.frame(cz_sites[-1])

# Add empty publications column
cz_sitedf$publications <- NA
cz_sitedf$DOI <- NA

# Loop through each siteid
for (i in 1:nrow(cz_sitedf)) {
  site_id <- cz_sitedf$siteid[i]
  
  # Get publications for this site
  pub_result <- get_publications(siteid = site_id)
  
  # Check if successful
  if (length(pub_result)>0) {
    cz_sitedf$publications[i] <- pub_result@publications[[1]]@citation
    cz_sitedf$DOI[i] <- pub_result@publications[[1]]@doi
  } else {
    cz_sitedf$publications[i] <- "Error or No publications"
  }
}

# Save

write.csv(cz_sitedf, paste0(mainDir,RegionName,"/Site_metadata_",RegionName,".csv"))


# Create a data frame with all unique combinations of variablename and ecologicalgroup
unique_combinations <- cbind(taxa = allSamp$variablename, taxa_clean = clean_taxa_names(allSamp$variablename, unique_only = F), ecologicalgroup = allSamp$ecologicalgroup)
unique_combinations <- unique_combinations[!duplicated(unique_combinations[,c('taxa','taxa_clean','ecologicalgroup')]),]
unique_combinations <- unique_combinations[order(unique_combinations[,"taxa"]),]

write.csv(unique_combinations[-1], paste0(mainDir,RegionName,"/Taxon_ecologicalGroup_",RegionName,".csv"))

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

#############      STEP 2   ############
######## ORGANIZE NEOTOMA DATA #########
########################################

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


#############    STEP 3       ###############
####### PALEOFIRE PACKAGE STARTS HERE ########
##############################################

#library(paleofire)
#mainDir <-  "C:/"   #Ideally same directory as used in Get_NeotomaData.R

proxy = "TRSH" # CHAR OR TRSH
nameregion = RegionName # NAME OF THE REGION


#Set parameters for generating the composite curve
TarAge = 50 #Pre-binning window value
binhw = 100 #Bin half window value
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
     , ylim = c(-2,2)
)

N_records <- rowSums(!is.na(Data_comp$BinnedData))
comp_data <- cbind(Data_comp$Result, N_records)

write.csv(comp_data,paste0(mainDir,nameregion,"/",nameregion,"_",proxy,"_hw",hw,"_binhw",binhw,"_tarAge",TarAge,".csv"))


#######   STEP 4   ############
##### SPATIALIZED VIEW #####

spatialtable <- Data_comp$BinnedData
colnames(spatialtable) <- file_name
spatialtable <- as.data.frame(spatialtable)
rownames(spatialtable) <- Data_comp$Result$AGE


# STEP 4.1: Define your time intervals (start and end ages)
my_intervals <- data.frame(
  start = c(0, 4200, 8200, 11700, 12800, 14800,18000),      # Start ages
  end = c(4200, 8200, 11700, 12800, 14800,18000,21000)     # End ages
)
# STEP 4.2: Define your labels
my_labels <- c("LH", "MH", "EH", "YD", "BA", "HS1", "LGM")  # One label per interval



# STEP 4.3: Function to apply intervals to your data
apply_time_intervals <- function(data, ages, intervals, labels = NULL) {
  
  # If no labels provided, create generic ones
  if (is.null(labels)) {
    labels <- paste0("interval_", 1:nrow(intervals))
  }
  
  # Check length matches
  if (length(labels) != nrow(intervals)) {
    stop("Number of labels must match number of intervals")
  }
  
  # Calculate means for each interval
  interval_means <- list()
  
  for (i in 1:nrow(intervals)) {
    start_age <- intervals$start[i]
    end_age <- intervals$end[i]
    
    # Find rows in this interval
    idx <- which(ages >= start_age & ages <= end_age)
    
    if (length(idx) > 0) {
      interval_means[[labels[i]]] <- colMeans(data[idx, , drop = FALSE], na.rm = TRUE)
    } else {
      interval_means[[labels[i]]] <- rep(NA, ncol(data))
    }
  }
  
  # Convert to data frame
  interval_df <- as.data.frame(do.call(rbind, interval_means))
  
  return(interval_df)
}

# STEP 4.4: Apply to your data
interval_means <- apply_time_intervals(
  data = spatialtable,
  ages = Data_comp$Result$AGE,
  intervals = my_intervals,
  labels = my_labels
)


# STEP 4.5: Add to your spatialtable
#spatialtable1 <- rbind(interval_means, spatialtable)
spatialtable1 <- t(interval_means)

# Extract siteid from row names
siteid <- sub(".*_", "", rownames(spatialtable1))

# Prepare sites dataframe with coordinates
sites_latlong <- as.data.frame(cz_sites)[, c("siteid", "lat", "long")]

# Create final dataframe directly
spatialtable_Comb <- data.frame(
  siteid = siteid,
  lat = sites_latlong$lat[match(siteid, sites_latlong$siteid)],
  long = sites_latlong$long[match(siteid, sites_latlong$siteid)],
  as.data.frame(spatialtable1),
  row.names = rownames(spatialtable1)  # Keep original row names
)


write.csv(spatialtable_Comb,paste0(mainDir,nameregion,"/",nameregion,"_",proxy,"_hw",hw,"_binhw",binhw,"_tarAge",TarAge,"_spatialtable_Comb.csv"), na = "")




#### FOR Z-SCORES

# Plot all maps based on your labels
for (label in my_labels) {
  if (label %in% colnames(spatialtable_Comb)) {
    print(
      ggplot(spatialtable_Comb, aes(x = long, y = lat)) +
        borders("world", fill = "gray95", color = "gray40") +
        geom_point(
          aes(fill = .data[[label]], color = .data[[label]]),  # Both fill and color
          size = 4, 
          shape = 21,  # Circle with both fill and color
          stroke = 0.3  # Thickness of the contour line
        ) +
        scale_fill_gradient2(
          low = "brown", 
          mid = "white", 
          high = "green", 
          midpoint = 0,
          na.value = "transparent"
        ) +
        scale_color_gradient2(
          low = "darkred",  # Darker version of brown
          mid = "gray50",   # Darker version of white
          high = "darkgreen",  # Darker version of green
          midpoint = 0,
          na.value = "transparent",
          guide = "none"  # Hide the color legend since fill is enough
        ) +
        theme_minimal() +
        labs(
          title = paste("Distribution of", label),
          fill = label
        ) +
        coord_fixed(
          xlim = range(spatialtable_Comb$long, na.rm = TRUE) + c(-2, 2),
          ylim = range(spatialtable_Comb$lat, na.rm = TRUE) + c(-2, 2)
        )
    )
  } else {
    warning(paste("Label", label, "not found in data"))
  }
}





##### FOR PERCENTAGES



# Plot all maps based on your labels
for (label in my_labels) {
  if (label %in% colnames(spatialtable_Comb)) {
    print(
      ggplot(spatialtable_Comb, aes(x = long, y = lat)) +
        borders("world", fill = "gray95", color = "gray40") +
        geom_point(
          aes(fill = .data[[label]]),  # Use fill for the color gradient
          size = 4, 
          shape = 21, 
          color = "black",  # Color for the contour line
          stroke = 0.1  # Thickness of the contour line
        ) +
        scale_fill_gradient2(
          low = "orange", 
          high = "forestgreen",
          midpoint = 50,
          na.value = "transparent"  # This will make NA values transparent
        ) +
        theme_minimal() +
        labs(
          title = paste("Distribution of", label),
          fill = label  # Change legend title to match the interval
        ) +
        coord_fixed(
          xlim = range(spatialtable_Comb$long, na.rm = TRUE) + c(-2, 2),
          ylim = range(spatialtable_Comb$lat, na.rm = TRUE) + c(-2, 2)
        )
    )
  } else {
    warning(paste("Label", label, "not found in data"))
  }
}
