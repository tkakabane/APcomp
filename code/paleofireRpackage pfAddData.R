#### Function from paleofire R package (Blarquez et al., 2014) (https://github.com/oblarquez/paleofire)

### pfaddData function


pfAddData <- function (files, metadata = NULL, type = "NULL", Int = TRUE, 
          first = NULL, last = NULL, yrInterp = NULL, sep = ",", dec = ".") {
  for (i in 1:length(files)) assign(paste("data", i, sep = ""), 
                                    read.csv(paste(files[i]), sep = sep, dec = dec))
  if (type == "CharAnalysis") {
    for (i in 1:length(files)) {
      temp <- pretreatment(get(paste("data", i, sep = ""))[, 
                                                           1:5], get(paste("data", i, sep = ""))[, 6], Int = Int, 
                           first = first, last = last, yrInterp = yrInterp)
      assign(paste("dataI", i, sep = ""), na.omit(as.data.frame(cbind(10000 + 
                                                                        i, temp$cmI, temp$ybpI, temp$accI))))
    }
    all_data <- do.call(rbind, lapply(1:length(files), function(i) rbind(get(paste("dataI", 
                                                                                   i, sep = "")))))
  }
  if (type == "NULL") {
    all_data <- do.call(rbind, lapply(1:length(files), function(i) rbind(cbind(10000 + 
                                                                                 i, get(paste("data", i, sep = ""))))))
  }
  colnames(all_data) <- c("ID_SITE", "DEPTH", "EST_AGE", "QUANTITY")
  if (is.null(metadata)) {
    meta <- data.frame(ID_SITE = c(1:length(files)) + 10000, 
                       SITE_NAME = rep(NA, length(files)), LATITUDE = rep(NA, 
                                                                          length(files)), LONGITUDE = rep(NA, length(files)))
  }
  if (is.null(metadata) == FALSE) {
    options(warn = -1)
    assign("meta", read.csv(paste(metadata), header = T))
    meta[, 1] <- as.character(meta[, 1])
    class(meta[, 1]) <- "character"
    meta <- cbind(unique(all_data[, 1]), meta)
    colnames(meta)[1] <- "ID_SITE"
  }
  output <- structure(list(data = all_data, metadata = meta))
  class(output) <- "pfAddData"
  return(output)
}
