library(mzR)

# load mzxml of interest
msfile <- openMSfile("K:/XCMS_ANALYSIS/All_Inga_Project/2_Nov_All_Inga_XCMS_Directories/BCI/IngA/Sample/IngA_1619.mzXML")
scandf <- create_spec_df(msfile)

# function first removes scans we definitely aren't interested in (early,late,etc), then puts all scans into a single dataframe where rows are scans, rownames are RT (in minutes), columns are masses, and values are TIC
create_spec_df <- function(msfile, rt_min = 0, rt_max = 40) {
ms1scans <- header(msfile)[header(msfile)$retentionTime <= rt_max*60 & header(msfile)$retentionTime >= rt_min*60 & header(msfile)$msLevel == 1,]
# create dataframe with rows as scans and columns as masses
scandf <- data.frame()
for(k in 1:length(ms1scans$seqNum)) {
current_scan <- mzR::peaks(msfile, ms1scans$seqNum[k])
current_scan <- current_scan[current_scan[,2]>=200,,drop=FALSE]
if(nrow(current_scan) < 1) next
scandf[as.character(ms1scans$retentionTime[k]),] <- rep(0, ncol(scandf))
for(i in 1:nrow(current_scan)) {
if(ncol(scandf) > 0) {
if(min(abs(current_scan[i,1] - as.numeric(names(scandf)))) <= 0.015) {
  scandf[as.character(ms1scans$retentionTime[k]),which.min(abs(current_scan[i,1] - as.numeric(names(scandf))))] <- current_scan[i,2]
}
  else {
    scandf[,as.character(current_scan[i,1])] <- rep(0,nrow(scandf))
    scandf[as.character(ms1scans$retentionTime[k]),as.character(current_scan[i,1])] <- current_scan[i,2]
  }
  } 
else {
  scandf[,as.character(current_scan[i,1])] <- rep(0,nrow(scandf))
  scandf[as.character(ms1scans$retentionTime[k]),as.character(current_scan[i,1])] <- current_scan[i,2]
}}
}
return(scandf)
}

# Function that returns the total TIC for a given mz/rt
# requres scandf = output of create_spec_df(), mz, rt_minutes, and minscans (minimum number of consecutive scans in which a mass must appear to count as a peak. default is 5.)
findpeak <- function(scandf, mz, rt_minutes, minscans = 5, override_peak_protection = FALSE) {
curr_peak <- scandf[, which(abs(mz - as.numeric(names(scandf))) <= 0.01), drop = FALSE]
if(ncol(curr_peak) > 1) {
  curr_peak <- curr_peak[, which.min(abs(mz - as.numeric(names(curr_peak)))), drop=FALSE]
}
if(ncol(curr_peak) == 0) {
  return(0) 
  break }
if(sum(curr_peak > 0)/nrow(curr_peak) > 0.1 & override_peak_protection == FALSE) {
  return(0)
  break}
for(i in 2:(nrow(curr_peak)-1)) {
  if(curr_peak[i-1,] > 0 & curr_peak[i+1,] > 0) curr_peak[i,] <- mean(c(curr_peak[i-1,],curr_peak[i+1,]))
}
zero_scans <- c(min(which(curr_peak != 0))-1, setdiff(min(which(curr_peak != 0)):max(which(curr_peak != 0)),which(curr_peak != 0)), max(which(curr_peak != 0))+1)
mass_peaks <- data.frame("firstScan" = character(0), "lastScan" = character(0), "peakScan" = character(0), stringsAsFactors = FALSE)
for(x in 2:length(zero_scans)) {
  if(zero_scans[x]-zero_scans[x-1] - 1 > minscans) {
    current_row <- nrow(mass_peaks) + 1
    mass_peaks[current_row,"firstScan"] <- row.names(curr_peak)[zero_scans[x-1] + 1]
    mass_peaks[current_row,"lastScan"] <- row.names(curr_peak)[zero_scans[x] - 1]
    mass_peaks[current_row,"peakScan"] <- row.names(curr_peak)[zero_scans[x-1] + which.max(curr_peak[(zero_scans[x-1]+1):(zero_scans[x]-1),])]
  }}
if(sum(abs(as.numeric(mass_peaks$peakScan) - rt_minutes*60)<30)>0){
peak_to_use <- which.min(abs(as.numeric(mass_peaks$peakScan) - rt_minutes*60))
return(sum(curr_peak[which(row.names(curr_peak)==mass_peaks$firstScan[peak_to_use]):which(row.names(curr_peak)==mass_peaks$lastScan[peak_to_use]),]))
}
else return(0)
}

