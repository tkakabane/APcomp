############################
library(paleofire)

mainDir <-  "C:/........"   #Ideally same directory as used in Get_NeotomaData.R

proxy = "TRSH" # CHAR OR TRSH
nameregion = "Neotropics"


#Set parameters for generating the composite curve
TarAge = 20 #Res
binhw = 40 #Bin half width value
hw = 1000 #Smoothing half width value
basePeriod = c(200, 21000)



TaddData = list.files(paste0(mainDir,nameregion,"/",nameregion," ",proxy,"/"), full.names = T)

TaddData <- TaddData[basename(TaddData) != "desktop.ini"]
file_name <- basename(TaddData)
file_name <- file_name[file_name != "desktop.ini"]
file_name <- gsub(".csv$", "", file_name)

trans_data1 <- pfAddData(
  files = TaddData, 
  type = "NULL")

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


write.csv(
  comp_data,
  paste0(
    mainDir,
    "/",
    nameregion,
    "_",
    proxy,
    "_",
    binhw,
    "_",
    TarAge,
    "_",
    hw,
    ".csv"
  ))



####
