######## ORGANIZE NEOTOMA DATA #########



####    Main directory
mainDir <- "C:/Users/Thomas/OneDrive/1 - Artigo/Akabane compilations/1 - Data for paper/AP_pollenR/raw data/"
RegionName <- "Neotropics"


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