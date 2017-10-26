setwd("K:/GY_LAB_FILES/github_repositories/chem_similarity/")
source("./code/fill_peaks_functions.R")
library(xcms)

files <- list.files("K://XCMS_ANALYSIS/All_Inga_Project/2_Nov_All_Inga_XCMS_Directories/", recursive = TRUE, pattern = "*mzXML", full.names = TRUE)
files <- grep("Sample", files, value = TRUE)

sample_IS <- data.frame("file" = sapply(1:length(files), function(i) unlist(strsplit(files[i], split="/"))[10]), "IS_TIC" = numeric(length(files)), stringsAsFactors = FALSE)

for(i in 1:nrow(sample_IS)) {
  sample <- sample_IS$file[i]
  print(sample)
  msfile <- openMSfile(files[i])
  scandf <- create_spec_df(msfile, rt_min = 21, rt_max = 24)
    if(min(abs(283.06 - as.numeric(names(scandf)))) < 0.01) {
      TIC <- findpeak(scandf, mz = 283.06, rt_minutes = 22.5, override_peak_protection = TRUE)
      if(TIC > 0) {
        sample_IS$IS_TIC[i] <- TIC
    }
  }
}

head(sample_IS)
range(sample_IS$IS_TIC)
hist(sample_IS$IS_TIC)
hist(sample_IS$IS_TIC[sample_IS$IS_TIC<100000])
table(sample_IS$IS_TIC)

species_IS <- data.frame(compound_tic_table,cSplit(compound_tic_table, 'compound_sample', sep="_", type.convert=FALSE))

raw_data <- readMSData(files[i], msLevel = 1)
extractMsData(raw_data, rt = c(21,24), mz = c(283.05, 283.07))
extractMsData
