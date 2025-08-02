library(paleofire)
#

# FUNCTION INCLUDES METHODS PRS - Proportional relative scaling | PRSBP - Proportional relative scaling with a base period


{pfTransform2 <-
    function (ID = NULL, add = NULL, Interpolate = FALSE, Age = NULL, 
              method = "Z-Score", BasePeriod = c(-100, 1e+09), span = 0.3, 
              RunWidth = 500, RunQParam = 0.5, stlYears = 500, type = "BoxCox1964", 
              alpha = 0.01, QuantType = "INFL", MethodType = NULL, verbose = TRUE) {
      paleofiresites <- NULL
      rm(paleofiresites)
      IDChar <- ID
      methods <- c("stl", "Z-Score", "Box-Cox", "LOESS", "MinMax", 
                   "RunMed", "RunMean", "RunMin", "RunMax", "RunQuantile", "PRS", "PRSBP",
                   "SmoothSpline", "Hurdle", "NULL")
      warnmethod <- method[(method %in% methods) == FALSE]
      if (length(warnmethod) != 0) {
        stop(paste(warnmethod, "is not a valid method for pfTransform", 
                   sep = " "))
      }
      types <- c("BoxCox1964", "JohnDraper")
      warntype <- type[(type %in% types) == FALSE]
      if (length(warntype) != 0) {
        stop(paste(warntype, "is not a valid type for pfBoxCox", 
                   sep = " "))
      }
      if (method == "RunMean" || method == "RunMin" || method == 
          "RunMed" || method == "RunMax" || method == "RunQuantile") {
        if ("caTools" %in% rownames(installed.packages()) == 
            FALSE) {
          install.packages("caTools")
        }
        if ("gtools" %in% rownames(installed.packages()) == FALSE) {
          install.packages("gtools")
        }
      }
      if (identical(method, "Hurdle")) {
        install.packages("pscl")
      }
      params <- list(ID = ID, Interpolate = Interpolate, Age = Age, 
                     method = method, BasePeriod = BasePeriod, span = span, 
                     RunWidth = RunWidth, RunQParam = RunQParam, stlYears = stlYears, 
                     type = type, alpha = alpha)
      if (verbose == TRUE) {
        cat("Loading and preparing data...")
        cat("\n")
      }
      if (is.null(ID) == FALSE) {
        if (is.list(ID) & length(ID) == 2) {
          data(paleofiredata, envir = environment())
          data(paleofiresites, envir = environment())
          ID <- ID$id_site
          paleofiredata <- paleofiredata[paleofiredata[, 1] %in% 
                                           ID, ]
          if (is.null(MethodType)) {
            for (i in ID) {
              paleofiredata[paleofiredata$ID_SITE %in% i & 
                              !(paleofiredata$UNIT %in% paleofiresites$pref_units[paleofiresites$id_site == 
                                                                                    i]), 7] <- NA
            }
            paleofiredata <- paleofiredata[!is.na(paleofiredata$TYPE), 
            ]
            if (QuantType == "INFL") {
              for (i in ID) if (!(unique(paleofiredata[paleofiredata[, 
                                                                     1] == i, 7]) %in% "INFL") & is.na(sum(paleofiredata[paleofiredata[, 
                                                                                                                                       1] == i, 2])) == FALSE & sum(paleofiredata[paleofiredata[, 
                                                                                                                                                                                                1] == i, 2]) > 0) {
                infl <- influx(paleofiredata[paleofiredata[, 
                                                           1] == i, ])
                paleofiredata[paleofiredata[, 1] == i, 4] <- c(infl)
              }
            }
          }
          else {
            paleofiredata <- paleofiredata[paleofiredata[, 
                                                         6] %in% MethodType, ]
            for (i in ID) {
              if (length(unique(paleofiredata[paleofiredata[, 
                                                            1] %in% i, 5])) >= 2) {
                paleofiredata[paleofiredata[, 1] %in% i & 
                                !(paleofiredata[, 5] %in% paleofiresites$pref_units[paleofiresites[, 
                                                                                                   1] == i]), 7] <- NA
              }
            }
            paleofiredata <- paleofiredata[!is.na(paleofiredata$TYPE), 
            ]
            cat(IDChar$site_name[!(IDChar$id_site %in% unique(paleofiredata[, 
                                                                            1]))], "\n")
            cat(length(IDChar$site_name) - length(unique(paleofiredata[, 
                                                                       1])), " sites were excluded from the analysis \n")
            ID <- unique(paleofiredata[, 1])
            if (QuantType == "INFL") {
              for (i in ID) if (!(unique(paleofiredata[paleofiredata[, 
                                                                     1] == i, 7]) %in% "INFL") & is.na(sum(paleofiredata[paleofiredata[, 
                                                                                                                                       1] == i, 2])) == FALSE) {
                infl <- influx(paleofiredata[paleofiredata[, 
                                                           1] == i, ])
                paleofiredata[paleofiredata[, 1] == i, 4] <- c(infl)
              }
            }
          }
          if (is.null(add) == FALSE) {
            add$data <- cbind(add$data, UNIT = NA, METHOD = NA, 
                              TYPE = "INFL")
            paleofiredata <- rbind(paleofiredata, add$data)
            ID <- c(ID, unique(add$data[, 1]))
          }
        }
      }
      if (is.null(ID)) {
        paleofiredata <- add$data
        ID <- c(unique(add$data[, 1]))
      }
      if (is.character(ID)) {
        paleofiredata <- read.csv(ID)
        ID <- unique(paleofiredata[, 1])
      }
      if (is.list(ID) & length(ID) > 2) {
        temp <- ID$TransData
        depths <- ID$IntDepths
        age <- ID$Age
        sites <- as.numeric(colnames(temp))
        ids <- matrix(nrow = length(temp[, 1]), ncol = length(temp[1, 
        ]))
        for (i in 1:length(temp[, 1])) {
          ids[i, ] <- sites
        }
        ids <- c(ids)
        data <- c(temp)
        age <- c(age)
        depths <- c(depths)
        if (length(depths) == 0) 
          depths <- rep(NA, length(age))
        paleofiredata <- cbind(ids, depths, age, data)
        ID <- unique(paleofiredata[, 1])
      }
      if (is.matrix(ID)) {
        paleofiredata <- ID
        ID <- unique(paleofiredata[, 1])
      }
      if (Interpolate == TRUE) {
        if (is.null(Age)) {
          res <- matrix(ncol = 1, nrow = length(ID))
          for (k in 1:length(ID)) {
            resT <- diff(paleofiredata[paleofiredata[, 1] == 
                                         ID[k], 3])
            res[k] <- c(median(resT[resT > 0]))
          }
          step <- round(median(res))
          minA <- round(min(paleofiredata[, 3]))
          maxA <- round(max(paleofiredata[, 3]))
          AgeN <- seq(minA, maxA, step)
        }
        if (is.null(Age) == FALSE) {
          AgeN <- Age
          ID <- unique(paleofiredata[, 1])
        }
        rawI <- matrix(nrow = length(AgeN), ncol = length(ID))
        for (k in 1:length(ID)) {
          if (length(paleofiredata[paleofiredata[, 1] == ID[k], 
                                   3]) >= 3) {
            rawI[, k] <- approx(paleofiredata[paleofiredata[, 
                                                            1] == ID[k], 3], paleofiredata[paleofiredata[, 
                                                                                                         1] == ID[k], 4], AgeN, method = "linear")$y
          }
          else {
            print(paste(IDChar$site_name[k], "has < 3 charcoal values and was excluded", 
                        sep = " "))
          }
        }
        depthI <- matrix(nrow = length(AgeN), ncol = length(ID))
        for (k in 1:length(ID)) {
          if (is.na(sum(paleofiredata[paleofiredata[, 1] == 
                                      ID[k], 2])) == F) {
            if (length(paleofiredata[paleofiredata[, 1] == 
                                     ID[k], 3]) >= 3) {
              depthI[, k] <- approx(paleofiredata[paleofiredata[, 
                                                                1] == ID[k], 3], paleofiredata[paleofiredata[, 
                                                                                                             1] == ID[k], 2], AgeN, method = "linear")$y
            }
          }
          else {
            depthI[, k] <- NA
          }
        }
        supp <- c()
        for (i in 1:length(rawI[1, ])) {
          if (sum(!is.na(rawI[, i])) < 3) {
            supp[i] <- 1
          }
          else {
            supp[i] <- 0
          }
        }
        rawI <- rawI[, supp == 0]
        SuppSites <- ID[supp == 1]
        ID <- ID[supp == 0]
        transI <- matrix(nrow = length(AgeN), ncol = length(ID))
        Ages <- matrix(ncol = length(ID), nrow = length(AgeN))
        for (k in 1:length(ID)) {
          Ages[, k] <- c(AgeN)
        }
      }
      if (Interpolate == FALSE) {
        lengths <- matrix(ncol = 1, nrow = length(ID))
        for (k in 1:length(ID)) {
          lengths[k] <- c(length(paleofiredata[paleofiredata[, 
                                                             1] %in% ID[k], 1]))
        }
        m <- max(lengths)
        transI <- matrix(nrow = m, ncol = length(ID))
        rawI <- matrix(nrow = m, ncol = length(ID))
        depthI <- matrix(nrow = m, ncol = length(ID))
        Ages <- matrix(ncol = length(ID), nrow = m)
        for (k in 1:length(ID)) {
          forNA <- m - length(paleofiredata[paleofiredata[, 
                                                          1] %in% ID[k], 3])
          AgeTemp <- c(paleofiredata[paleofiredata[, 1] %in% 
                                       ID[k], 3], rep(NA, forNA))
          Ages[, k] <- c(AgeTemp)
          rawTemp <- c(paleofiredata[paleofiredata[, 1] %in% 
                                       ID[k], 4], rep(NA, forNA))
          rawI[, k] <- c(rawTemp)
          depthTemp <- c(paleofiredata[paleofiredata[, 1] %in% 
                                         ID[k], 2], rep(NA, forNA))
          depthI[, k] <- c(depthTemp)
        }
      }
      if (verbose == TRUE) {
        percent <- seq(10, 100, by = 10)
        values <- round(percent * length(ID)/100)
        cat("Transforming...")
        cat("\n")
        cat("Percentage done: ")
      }
      pb <- txtProgressBar(0, length(method), style = 3)
      for (j in 1:length(method)) {
        methodj <- method[j]
        if (j >= 2) {
          rawI <- transI
        }
        for (k in 1:length(ID)) {
          tmp <- cbind(Ages[, k], rawI[, k])
          tmp <- na.omit(tmp)
          
          if (any(duplicated(tmp[, 1]))) {
            cat("\n")
            cat("\n")
            cat("Duplicate ages in ID:", ID[k], "at ages:", tmp[duplicated(tmp[, 1]), 1], "\n")
          }
          
          if (length(tmp[, 1]) > 3 & ID[k] != 882) {
            if (methodj == "stl") {
              agesI <- seq(tmp[1, 1], tmp[length(tmp[, 1]), 
                                          1], 1)
              forTS <- approx(tmp[, 1], tmp[, 2], agesI)$y
              x <- ts(forTS, start = 1, frequency = stlYears)
              dim(x) <- NULL
              stlResult <- stl(x, "per")$time.series[, 2]
              transI[, k] <- approx(agesI, stlResult, Ages[, 
                                                           k])$y
            }
            if (methodj == "NULL") {
              transI[, k] <- approx(tmp[, 1], tmp[, 2], Ages[, 
                                                             k])$y
            }
            if (methodj == "Z-Score") {
              mu <- mean(tmp[tmp[, 1] >= BasePeriod[1] & 
                               tmp[, 1] <= BasePeriod[2], 2])
              sigma <- sd(tmp[tmp[, 1] >= BasePeriod[1] & 
                                tmp[, 1] <= BasePeriod[2], 2])
              if (is.na(mu) | is.na(sigma) | sigma == 0) {
                transI[, k] <- approx(tmp[, 1], scale(tmp[, 
                                                          2]), Ages[, k])$y
              }
              else {
                transI[, k] <- approx(tmp[, 1], (tmp[, 2] - 
                                                   mu)/sigma, Ages[, k])$y
              }
            }
            if (methodj == "Box-Cox") {
              transI[, k] <- approx(tmp[, 1], pfBoxCox(tmp[, 
                                                           2], alpha = alpha, type = type), Ages[, k])$y
            }
            if (methodj == "LOESS") {
              transI[, k] <- approx(tmp[, 1], predict(loess(tmp[, 
                                                                2] ~ tmp[, 1], span = span)), Ages[, k])$y
            }
            if (methodj == "MinMax") {
              transI[, k] <- approx(tmp[, 1], pfMinMax(tmp[, 
                                                           2]), Ages[, k])$y
            }
            if (methodj == "PRS") {
              transI[, k] <- 100* (approx(tmp[, 1], tmp[, 2] / max(tmp[, 2], na.rm = TRUE)*(sum(tmp[, 2] != 0)/nrow(tmp)), Ages[, k])$y)
            }
            if (methodj == "PRSBP") {
              # Get values within BasePeriod
              base_idx <- tmp[,1] >= BasePeriod[1] & tmp[,1] <= BasePeriod[2]
              base_vals <- tmp[base_idx, 2]
              
              # Calculate fire frequency ONLY within BasePeriod
              non_zero_count <- sum(base_vals != 0, na.rm = TRUE)
              total_count <- sum(!is.na(base_vals))
              fire_freq <- ifelse(total_count > 0, non_zero_count/total_count, 0)
              
              # Calculate min and max from BasePeriod
              min_base <- min(base_vals, na.rm = TRUE)
              max_base <- max(base_vals, na.rm = TRUE)
              
              # Handle cases where rescaling isn't possible
              if (is.infinite(min_base) | is.infinite(max_base) | (max_base) == 0) {
                transI[, k] <- 100* (approx(tmp[, 1], tmp[, 2] / max(tmp[, 2], na.rm = TRUE)*(sum(tmp[, 2] != 0)/nrow(tmp)), Ages[, k])$y)
              } else {
                # Rescale using BasePeriod min/max and multiply by BasePeriod fire frequency
                prs_scaled <- ((tmp[,2]) / (max_base)) * fire_freq
                transI[, k] <- 100* approx(tmp[,1], prs_scaled, Ages[,k])$y
              }
            }
            if (methodj == "RunMed") {
              w <- round(RunWidth/((max(tmp[, 1]) - min(tmp[, 
                                                            1]))/length(tmp[, 1])))
              if (gtools::odd(w)) 
                w <- w
              else w <- w + 1
              transI[, k] <- approx(tmp[, 1], runmed(tmp[, 
                                                         2], w), Ages[, k])$y
            }
            if (methodj == "RunMean") {
              w <- round(RunWidth/((max(tmp[, 1]) - min(tmp[, 
                                                            1]))/length(tmp[, 1])))
              if (gtools::odd(w)) 
                w <- w
              else w <- w + 1
              transI[, k] <- approx(tmp[, 1], caTools::runmean(tmp[, 
                                                                   2], w), Ages[, k])$y
            }
            if (methodj == "RunMin") {
              w <- round(RunWidth/((max(tmp[, 1]) - min(tmp[, 
                                                            1]))/length(tmp[, 1])))
              if (gtools::odd(w)) 
                w <- w
              else w <- w + 1
              transI[, k] <- approx(tmp[, 1], caTools::runmin(tmp[, 
                                                                  2], w), Ages[, k])$y
            }
            if (methodj == "RunMax") {
              w <- round(RunWidth/((max(tmp[, 1]) - min(tmp[, 
                                                            1]))/length(tmp[, 1])))
              if (gtools::odd(w)) 
                w <- w
              else w <- w + 1
              transI[, k] <- approx(tmp[, 1], caTools::runmax(tmp[, 
                                                                  2], w), Ages[, k])$y
            }
            if (methodj == "RunQuantile") {
              w <- round(RunWidth/((max(tmp[, 1]) - min(tmp[, 
                                                            1]))/length(tmp[, 1])))
              if (gtools::odd(w)) 
                w <- w
              else w <- w + 1
              transI[, k] <- approx(tmp[, 1], caTools::runquantile(tmp[, 
                                                                       2], w, RunQParam), Ages[, k])$y
            }
            if (methodj == "SmoothSpline") {
              transI[, k] <- approx(tmp[, 1], smooth.spline(tmp[, 
                                                                1], tmp[, 2], spar = span)$y, Ages[, k])$y
            }
            if (methodj == "Hurdle") {
              tmp[, 2] <- round(pfMinMax(tmp[, 2]) * 100)
              transI[, k] <- approx(tmp[, 1], pscl::hurdle(tmp[, 
                                                               2] ~ tmp[, 1])$fitted.values, Ages[, k])$y
            }
          }
          Sys.sleep(2e-05)
          if (k %in% seq(0, length(ID) * length(method), 1) & 
              verbose == TRUE) {
            setTxtProgressBar(pb, k)
          }
        }
      }
      if (verbose == TRUE) 
        cat("\n")
      colnames(transI) <- ID
      output <- structure(list(Age = structure(Ages, col.names = as.character(ID), 
                                               class = "matrix"), IntDepths = structure(depthI, col.names = as.character(ID), 
                                                                                        class = "matrix"), IntData = structure(rawI, col.names = as.character(ID), 
                                                                                                                               class = "matrix"), TransData = structure(transI, col.names = as.character(ID), 
                                                                                                                                                                        class = "matrix"), Method = method, params = params))
      class(output) <- "pfTransform2"
      return(output)
    }
}
#