source("./code/fill_peaks_functions.R")

ms_features <- read.csv("./data/all_features_with_mzrt.csv")

files <- list.files("K://XCMS_ANALYSIS/All_Inga_Project/2_Nov_All_Inga_XCMS_Directories/", recursive = TRUE, pattern = "*mzXML", full.names = TRUE)
files <- grep("Sample", files, value = TRUE)

feature_tics <- list("sample" = character(length(files)*nrow(ms_features)), "feature" = character(length(files)*nrow(ms_features)), "TIC" = numeric(length(files)*nrow(ms_features)), stringsAsFactors = FALSE)
current_row <- 0
for(i in 1:length(files)) {
  sample <- sub(".mzXML", "", unlist(strsplit(files[i], split="/"))[10])
  print(sample)
  msfile <- openMSfile(files[i])
  scandf <- create_spec_df(msfile)
  for(j in 1:nrow(ms_features)) {
    if(min(abs(ms_features$mz[j] - as.numeric(names(scandf)))) < 0.01) {
      TIC <- findpeak(scandf, mz = ms_features$mz[j], rt_minutes = ms_features$rt[j])
      if(TIC > 0) {
        current_row <- current_row + 1
        feature_tics[current_row, ] <- c(sample, ms_features$feature_number[j], TIC)
      }
    }
  }
}


feature_tics_2 <- feature_tics[feature_tics$sample != "",]
tail(feature_tics)
table(feature_tics_2$sample)
head(feature_tics, 30)
