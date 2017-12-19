### can I just use getPekas() function from XCMS to fill peaks?
library(xcms)

# load list of all features with associated mz and rt
ms_features <- read.csv("./data/samples_with_features_mzrt.csv")
features <- data.frame("feature_number" = unique(ms_features$feature_number), mz = numeric(length(unique(ms_features$feature_number))), rt = numeric(length(unique(ms_features$feature_number))))
for(i in 1:nrow(features)) {
  features$mz[i] <- mean(ms_features[ms_features$feature_number == features$feature_number[i], "mz"])
  features$rt[i] <- mean(ms_features[ms_features$feature_number == features$feature_number[i], "rt"])
}

features$mzmin <- features$mz*0.999975
features$mzmax <- features$mz*1.000025
features$rtmin <- features$rt*60 - 30
features$rtmax <- features$rt*60 + 30
ms_features_matrix <- as.matrix(features[, c("mzmin","mzmax","rtmin","rtmax")])

# load file that associates sample names with blank names
sample_blank <- read.csv("K:/XCMS_ANALYSIS/All_Inga_Project/sample_blank_association.csv")

sites <- c("BCI","FG", "LA", "Manaus", "Tiputini")

# if you want to run a specific species, fill out these two vectors with species and associated site.
# Then change sites[j] in line under 'for(k in 1:length(species)) to sites[k] and run all but outermost loop
species <- c("LA4")
sites <- c("LA")

for(j in 1:length(sites)) {
  species <- list.files(paste("K:/XCMS_ANALYSIS/All_Inga_Project/2_Nov_All_Inga_XCMS_Directories/",sites[j],"/",sep=""))
  species <- species[!startsWith(species, "X")]
  species <- species[!endsWith(species, "undiluted")]
  for(k in 1:length(species)) {
    samples <- list.files(paste("K:/XCMS_ANALYSIS/All_Inga_Project/2_Nov_All_Inga_XCMS_Directories/",sites[j],"/",species[k],"/Sample/",sep=""), full.names = TRUE)
    for(l in 1:length(samples)) {
    # fill peaks for sample
      print(samples[l])
      sample_xcms <- xcmsRaw(samples[l], profstep = 1, profmethod = , includeMSn = FALSE, mslevel = 1)
sample_peaks <- getPeaks(object = sample_xcms, peakrange = ms_features_matrix, step = 0.05)

sample_peaks_1 <- as.data.frame(sample_peaks)
sample_peaks_2 <- sample_peaks_1[sample_peaks_1$into >= 1000,]
sample_peaks_2$mass_matches <- character(nrow(sample_peaks_2))
sample_peaks_2$feature_number <- numeric(nrow(sample_peaks_2))
sample_peaks_2$feature_mz <- numeric(nrow(sample_peaks_2))
sample_peaks_2$feature_rt <- numeric(nrow(sample_peaks_2))
sample_peaks_2$actual_mz <- numeric(nrow(sample_peaks_2))
sample_peaks_2$actual_rt <- numeric(nrow(sample_peaks_2))
sample_peaks_2$TIC <- numeric(nrow(sample_peaks_2))
# make sure mass actually matches
for(i in 1:nrow(sample_peaks_2)) {
  feature_number <- features[features$mzmin == sample_peaks_2$mzmin[i] & features$rtmin == sample_peaks_2$rtmin[i], "feature_number"]
  sample_peaks_2$feature_number[i] <- feature_number
  sample_peaks_2$feature_mz[i] <- features[features$feature_number == feature_number, "mz"]
  sample_peaks_2$feature_rt[i] <- features[features$feature_number == feature_number, "rt"]
  #peakscan <- as.numeric(sample_xcms@scanindex[which.min(abs(sample_xcms@scantime - sample_peaks_2$rt[i]))])
  peakscan <- which.min(abs(sample_xcms@scantime - sample_peaks_2$rt[i]))
  test2 <- getScan(sample_xcms, peakscan, mzrange = c(sample_peaks_2$mzmin[i],sample_peaks_2$mzmax[i]))
  if(length(test2) == 0) {
    sample_peaks_2$mass_matches[i] <- "mismatch"
    next }
  # use 25 ppm error?
  if(abs(sample_peaks_2$feature_mz[i] - test2[1,1]) * 1000000 / sample_peaks_2$feature_mz[i]  < 25) {
    intensities <- sapply((peakscan-20):(peakscan+20), function(x) {
      tempscan = getScan(sample_xcms, x, mzrange = c(sample_peaks_2$mzmin[i],sample_peaks_2$mzmax[i]))
      ifelse(nrow(tempscan) == 0, 0, tempscan[1,2]) } )
    midpoint <- length(intensities)/2 + 0.5
    if(intensities[midpoint] / intensities[1] > 1.2 & 
       intensities[midpoint] / intensities[length(intensities)] > 1.2 &
       sum((midpoint-3):(midpoint+3) %in% which(intensities == 0)) <= 1) {
    diffs <- which(intensities == 0) - midpoint
    if(sum(diffs<0) == 0) first_scan = -20
    else first_scan <- diffs[diffs<0][sum(diffs<0)]
    if(sum(diffs>0) == 0) last_scan = 20
    else last_scan <- diffs[diffs>0][1]
      sample_peaks_2$mass_matches[i] <- "match"
    sample_peaks_2$actual_mz[i] <- test2[1,1]
    sample_peaks_2$actual_rt[i] <- sample_peaks_2$rt[i]/60
    sample_peaks_2$TIC[i] <- sum(intensities[(midpoint+first_scan):(midpoint+last_scan)])}
    else sample_peaks_2$mass_matches[i] <- "noise" }
  else sample_peaks_2$mass_matches[i] <- "mismatch"
}
sample_peaks_3 <- sample_peaks_2[sample_peaks_2$mass_matches == "match", c("feature_number","feature_mz","feature_rt", "into", "actual_mz", "actual_rt", "TIC")]

sample_name <- unlist(strsplit(samples[l], split = "/"))[8]
blank_name <- as.character(sample_blank[sample_blank$newname == sample_name, "Blank"])
# also check associated blank for same peaks
blank_xcms <- xcmsRaw(paste("K:/Lab_Map/Database/2015_Relational_DB/DATA_STORAGE/UPLC_MS_DATA/Data_Not_Active_Projects_etc/10_Converted_Data/2014_all_Inga.PRO/mzXML/Blank/", blank_name, ".mzXML", sep = ""), profstep = 1, profmethod = , includeMSn = FALSE, mslevel = 1)
blank_findpeaks <- getPeaks(object = blank_xcms, peakrange = ms_features_matrix, step = 0.05)
blank_findpeaks_1 <- as.data.frame(blank_findpeaks)
blank_findpeaks_2 <- blank_findpeaks_1[blank_findpeaks_1$into > 100,]
blank_findpeaks_2$mass_matches <- character(nrow(blank_findpeaks_2))
blank_findpeaks_2$feature_number <- numeric(nrow(blank_findpeaks_2))
blank_findpeaks_2$feature_mz <- numeric(nrow(blank_findpeaks_2))
blank_findpeaks_2$feature_rt <- numeric(nrow(blank_findpeaks_2))
blank_findpeaks_2$TIC <- numeric(nrow(blank_findpeaks_2))
# make sure mass actually matches
for(i in 1:nrow(blank_findpeaks_2)) {
  feature_number <- features[features$mzmin == blank_findpeaks_2$mzmin[i] & features$rtmin == blank_findpeaks_2$rtmin[i], "feature_number"]
  blank_findpeaks_2$feature_number[i] <- feature_number
  blank_findpeaks_2$feature_mz[i] <- features[features$feature_number == feature_number, "mz"]
  blank_findpeaks_2$feature_rt[i] <- features[features$feature_number == feature_number, "rt"]
  peakscan <- as.numeric(blank_xcms@scanindex[which.min(abs(blank_xcms@scantime - blank_findpeaks_2$rt[i]))])
  peakscan <- which.min(abs(blank_xcms@scantime - blank_findpeaks_2$rt[i]))
  test2 <- getScan(blank_xcms, peakscan, mzrange = c(blank_findpeaks_2$mzmin[i],blank_findpeaks_2$mzmax[i]))
  if(length(test2) == 0) {
    blank_findpeaks_2$mass_matches[i] <- "mismatch"
    next }
  if(abs(blank_findpeaks_2$feature_mz[i] - test2[1,1]) * 1000000 / blank_findpeaks_2$feature_mz[i]  < 25) {
    intensities <- sapply((peakscan-20):(peakscan+20), function(x) {
      tempscan = getScan(blank_xcms, x, mzrange = c(blank_findpeaks_2$mzmin[i],blank_findpeaks_2$mzmax[i]))
      ifelse(nrow(tempscan) == 0, 0, tempscan[1,2]) } )
    midpoint <- length(intensities)/2 + 0.5
      diffs <- which(intensities == 0) - midpoint
      if(sum(diffs<0) == 0) first_scan = -20
      else first_scan <- diffs[diffs<0][sum(diffs<0)]
      if(sum(diffs>0) == 0) last_scan = 20
      else last_scan <- diffs[diffs>0][1]
      blank_findpeaks_2$mass_matches[i] <- "match"
      blank_findpeaks_2$TIC[i] <- sum(intensities[(midpoint+first_scan):(midpoint+last_scan)])}
  else blank_findpeaks_2$mass_matches[i] <- "mismatch"
}
blank_findpeaks_3 <- blank_findpeaks_2[blank_findpeaks_2$mass_matches == "match",c("feature_number","feature_mz","feature_rt", "TIC")]

# delete all peaks that are not at least 5x as abundant in blank as in sample.
sample_blank_peaks <- merge(sample_peaks_3, blank_findpeaks_3[,c("feature_number", "TIC")], by = "feature_number", all.x = TRUE, all.y = FALSE)
sample_blank_peaks$TIC.y[is.na(sample_blank_peaks$TIC.y)] <- 0
sample_peaks <- sample_blank_peaks[sample_blank_peaks$TIC.x / sample_blank_peaks$TIC.y > 5, ]

sample_peaks_1 <- data.frame("feature_number" = sample_peaks$feature_number, "TIC" = sample_peaks$TIC.x, "actual_mz" = sample_peaks$actual_mz, "actual_rt" = sample_peaks$actual_rt, "sample_name" = unlist(strsplit(sample_name, split = "[.]"))[1])

if(j == 1 & k == 1 & l == 1) {
  write.table(sample_peaks_1, "./results/all_inga_filled_features_ppm_2017_12_05.csv", sep = ",", append = FALSE, row.names = FALSE, col.names = TRUE)
}
else {
write.table(sample_peaks_1, "./results/all_inga_filled_features_ppm_2017_12_05.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE) }
    }}}
