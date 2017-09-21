library(mzR)
library(MSnbase)

# Function finds weighted average of multiple msms scans. Requires MSnExp object containing MSMS files for species containing peak of interest as well as mz and rt of peak.
avg.msms.spec <- function(data, rt, mz) {
  rownames <- row.names(fData(ms2spec))[abs(fData(ms2spec)$retentionTime - rt*60) <= 20 & abs(fData(ms2spec)$precursorMZ - mz) <= 0.01]
  spec1 <- lapply(1:length(rownames), function(j) as(data[[rownames[j]]], "data.frame"))
  spec2 <- do.call(rbind, spec1)
  spec3 <- spec2[order(spec2$i, decreasing = T),]
  row.names(spec3) <- 1:nrow(spec3)
  
  rows_done = numeric()
  spec4 <- array(numeric(), dim=c(0,2))
  while(length(setdiff(1:nrow(spec3), rows_done)) > 0) {
    i = min(setdiff(1:nrow(spec3), rows_done))
    to_merge = setdiff(which(abs(spec3[,1] - spec3[i,1]) <= 0.25), rows_done)
    
    mz = weighted.mean(spec3[to_merge, 1], spec3[to_merge, 2])
    TIC = sum(spec3[to_merge, 2])
    
    spec4 <- rbind(spec4, c(mz, TIC))
    
    rows_done <- c(rows_done, to_merge)
  }
  spec5 <- spec4[order(spec4[,1]),]
  spec5 <- spec5[spec5[,2] >= 500,,drop=FALSE]
  return(spec5)
}
