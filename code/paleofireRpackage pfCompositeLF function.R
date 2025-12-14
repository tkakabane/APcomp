# pfCompositeLF from paleofire R Package


function (TR, hw = 250, tarAge = NULL, binhw = NULL, nboot = 1000, 
    conf = c(0.05, 0.95), pseudodata = FALSE, verbose = TRUE) 
{
    if (is.matrix(TR) | is.data.frame(TR)) {
        ID <- unique(TR[, 1])
        lgth <- c()
        for (i in 1:length(ID)) {
            lgth[i] <- length(na.omit(TR[TR[, 1] == ID[i], 1]))
        }
        m <- max(lgth)
        Age <- matrix(nrow = m, ncol = length(ID))
        TransData <- matrix(nrow = m, ncol = length(ID))
        for (i in 1:length(ID)) {
            Age[, i] <- c(TR[TR[, 1] == ID[i], 3], rep(NA, m - 
                length(TR[TR[, 1] == ID[i], 3])))
            TransData[, i] <- c(TR[TR[, 1] == ID[i], 4], rep(NA, 
                m - length(TR[TR[, 1] == ID[i], 4])))
        }
        colnames(TransData) <- ID
        colnames(Age) <- ID
        TR <- structure(list(Age = structure(Age, class = "matrix"), 
            TransData = structure(TransData, class = "matrix"), 
            Method = "unspecified"))
    }
    if (is.null(tarAge)) {
        AgeRes <- matrix(nrow = length(TR$Age[, 1]) - 1, ncol = length(TR$Age[1, 
            ]))
        for (i in 1:length(TR$Age[1, ])) {
            AgeRes[, i] <- c(diff(TR$Age[, i]))
        }
        width <- ceiling(median(na.omit(AgeRes))/10) * 10
        binI <- floor(min(na.omit(TR$Age))/10) * 10
        binF <- ceiling(max(na.omit(TR$Age))/10) * 10
        tarAge <- seq(binI, binF, width)
        if (is.null(binhw)) {
            binhw <- (tarAge[2] - tarAge[1])/2
        }
    }
    m <- length(TR$TransData[, 1])
    n <- length(TR$TransData[1, ])
    result <- matrix(ncol = length(TR$Age[1, ]), nrow = length(tarAge))
    if (is.null(binhw)) {
        binhw <- (tarAge[2] - tarAge[1])/2
    }
    if (verbose == TRUE) {
        percent <- seq(1, 100, by = 1)
        values <- round(percent * n/100)
        cat("Prebinning...")
        cat("\n")
    }
    pb <- txtProgressBar(1, n, style = 3)
    for (k in 1:n) {
        if (length(TR$TransData[is.na(TR$TransData[, k]) == FALSE, 
            k]) != 0) {
            for (i in 1:length(tarAge)) {
                t <- na.omit(cbind(as.numeric(TR$Age[, k]), as.numeric(TR$TransData[, 
                  k])))
                result[i, k] <- mean(t[t[, 1] > tarAge[i] - binhw & 
                  t[, 1] < tarAge[i] + binhw, 2])
            }
        }
        Sys.sleep(2e-05)
        if (k %in% values & verbose == TRUE) {
            setTxtProgressBar(pb, k)
        }
    }
    if (verbose == TRUE) {
        cat("\n")
    }
    result[!is.finite(result)] <- NA
    centres <- tarAge
    mboot <- matrix(nrow = length(centres), ncol = nboot)
    set.seed(1)
    if (verbose == TRUE) {
        percent <- seq(1, 100, by = 1)
        values <- round(percent * nboot/100)
        cat("Bootstrapping...")
        cat("\n")
        cat("Percentage done: ")
    }
    pb <- txtProgressBar(1, nboot, style = 3)
    for (i in 1:nboot) {
        ne <- sample(seq(1, length(result[1, ]), 1), length(result[1, 
            ]), replace = TRUE)
        y <- as.vector(result[, ne])
        x <- as.vector(rep(centres, length(ne)))
        dat <- na.omit(data.frame(x = x, y = y))
        if (pseudodata == TRUE) {
            pseudo_n <- round(length(dat$x) * 0.3)
            pseudo_up <- -((dat$x) - rep((min(dat$x)), length(dat$x))) + 
                min(dat$x)
            pseudo_up <- as.data.frame(cbind(pseudo_up, dat$y))
            pseudo_up <- pseudo_up[1:pseudo_n, ]
            pseudo_up <- pseudo_up[length(pseudo_up[, 1]):1, 
                ]
            colnames(pseudo_up) <- c("x", "y")
            pseudo_lo <- -((dat$x) - rep((max(dat$x)), length(dat$x))) + 
                max(dat$x)
            pseudo_lo <- as.data.frame(cbind(pseudo_lo, dat$y))
            pseudo_lo <- pseudo_lo[(length(dat$x) - pseudo_n):length(dat$x), 
                ]
            pseudo_lo <- pseudo_lo[length(pseudo_lo[, 1]):1, 
                ]
            colnames(pseudo_lo) <- c("x", "y")
            dat <- rbind(pseudo_up, dat, pseudo_lo)
        }
        x <- as.vector(dat$x)
        y <- as.vector(dat$y)
        locboot <- locfit(y ~ lp(x, deg = 1, h = hw), maxk = 2000, 
            family = "qrgauss")
        if (is.na(locboot$dp[7]) == FALSE) {
            predboot <- predict(locboot, newdata = centres, se.fit = TRUE)
            mboot[, i] <- predboot$fit
        }
        Sys.sleep(2e-05)
        if (i %in% values & verbose == TRUE) {
            setTxtProgressBar(pb, i)
        }
    }
    cat("\n")
    bootci <- t(apply(mboot, 1, quantile, probs = conf, na.rm = TRUE))
    bootmean <- t(apply(mboot, 1, mean, na.rm = TRUE))
    bootmed <- t(apply(mboot, 1, median, na.rm = TRUE))
    rm(x, y, dat)
    y <- c(result)
    x <- rep(centres, length(result[1, ]))
    dat <- as.data.frame(cbind(x, y))
    dat <- na.omit(dat[order(x), ])
    if (pseudodata == TRUE) {
        pseudo_n <- round(length(dat$x) * 0.3)
        pseudo_up <- -((dat$x) - rep((min(dat$x)), length(dat$x))) + 
            min(dat$x)
        pseudo_up <- as.data.frame(cbind(pseudo_up, dat$y))
        pseudo_up <- pseudo_up[1:pseudo_n, ]
        pseudo_up <- pseudo_up[length(pseudo_up[, 1]):1, ]
        colnames(pseudo_up) <- c("x", "y")
        pseudo_lo <- -((dat$x) - rep((max(dat$x)), length(dat$x))) + 
            max(dat$x)
        pseudo_lo <- as.data.frame(cbind(pseudo_lo, dat$y))
        pseudo_lo <- pseudo_lo[(length(dat$x) - pseudo_n):length(dat$x), 
            ]
        pseudo_lo <- pseudo_lo[length(pseudo_lo[, 1]):1, ]
        colnames(pseudo_lo) <- c("x", "y")
        dat <- rbind(pseudo_up, dat, pseudo_lo)
    }
    x <- as.vector(dat$x)
    y <- as.vector(dat$y)
    locbootA <- locfit(y ~ lp(x, deg = 1, h = hw), maxk = 500, 
        family = "qrgauss")
    predbootA <- predict(locbootA, newdata = centres, se.fit = TRUE)
    locfitAll <- predbootA$fit
    result2 <- as.data.frame(cbind(centres, locfitAll, t(bootmean), 
        bootci))
    colnames(result2) <- c("AGE", "LocFit", "MEAN(of_boot)", 
        as.character(conf))
    colnames(result) <- colnames(TR$TransData)
    colnames(result) <- colnames(TR$TransData)
    output <- structure(list(BinnedData = structure(result, row.names = as.character(centres), 
        col.names = colnames(TR$TransData), class = "matrix"), 
        Result = result2, mboot = mboot, BinCentres = centres, 
        BinWidth = binhw * 2, nboot = nboot, halfwidth = hw, 
        conf = conf, locfitAll = locfitAll, BootMed = structure(bootmed, 
            row.names = as.character(centres), col.names = c("Median"), 
            class = "matrix"), BootMean = structure(bootmean, 
            row.names = as.character(centres), col.names = c("Mean"), 
            class = "matrix"), BootCi = structure(bootci, row.names = as.character(centres), 
            class = "matrix")))
    class(output) <- "pfCompositeLF"
    return(output)
    if (verbose == TRUE) {
        cat("\n")
    }
}
