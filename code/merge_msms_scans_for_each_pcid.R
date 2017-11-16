library(RMySQL)
library(splitstackshape)
source("./code/average_multiple_msms_scans.R")

compound_table <- read.csv("./data/filled_compound_table_2017_11_15.csv")
compound_table$sample <- as.character(compound_table$sample)
compound_table$species <-  sapply(1:nrow(compound_table), function(x) unlist(strsplit(compound_table$sample[x], split="_"))[1])
compound_table_2 <- aggregate(compound_table$TIC, by=list(compound_table$species, compound_table$compound), FUN=mean)
names(compound_table_2) <- c("species","compound","TIC")

compound_feature <- read.csv("./data/compound_feature_table.csv")
compound_max_feature <- aggregate(compound_feature$TIC, by=list(Category=compound_feature$compound_number),FUN=max)
names(compound_max_feature) <- c("compound","TIC_max")
compound_max_feature_2 <- merge(compound_max_feature, compound_feature, by.x=c("compound","TIC_max"), by.y=c("compound_number","TIC"))

compound_table_3 <- merge(compound_table_2, compound_max_feature_2, by="compound", all.x = FALSE)
head(compound_table_3)


# if available, get spectrum for each compound
species <- unique(compound_table_3$species)
msms_spec <- vector("list", max(compound_table_3$compound))

for(i in 1:length(species)) {
  print(species[i])
  if(exists("mydb")) dbDisconnect(mydb)
  mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
  files <- dbGetQuery(mydb, paste("SELECT File_Name, Project_Name FROM UPLC_Results WHERE species_code = '",species[i],"' AND Project_Name LIKE '%DDA%'", sep=""))
  if(nrow(files) < 1) next
  file.paths <- sapply(1:nrow(files), function(x) paste("K:/Lab_Map/Database/2015_Relational_DB/DATA_STORAGE/UPLC_MS_DATA/Data_Not_Active_Projects_etc/10_Converted_Data/",files$Project_Name[x],"/mzXML/Sample/",files$File_Name[x],".mzXML", sep=""))
  peakdata <- extract.peakdata(file.paths, MSlevel = 2)
  spec_i <- compound_table_3[compound_table_3$species == species[i], ]
  for(k in 1:nrow(spec_i)) {
    rt_k <- spec_i$rt[k]
    mz_k <- spec_i$mz[k]
    current_spec <- avg.msms.spec(file.paths, peakdata, rt = rt_k, mz = mz_k)
    if(is.null(current_spec)) next
    if(nrow(current_spec) < 5) next
    if(max(msms_spec[[spec_i$compound[k]]][,2]) > max(current_spec[,2])) next
    msms_spec[[spec_i$compound[k]]] <- current_spec
  }
}

# which compounds don't have associated spectra?
no_spec <- which(sapply(1:length(msms_spec), function(x) is.null(msms_spec[[x]])))
no_spec <- no_spec[no_spec > 2000]
for(i in no_spec) {
  print(i)
  compound_table_i <- compound_table_3[compound_table_3$compound == i,]
  if(nrow(compound_table_i) == 0) next
  features_i <- compound_feature[compound_feature$compound_number == i,]
  features_i$relTIC <- features_i$TIC / max(features_i$TIC)
  features_i <- features_i[features_i$relTIC > 0.33 & features_i$relTIC != 1,]
  if(nrow(features_i) == 0) next 
  features_i <- features_i[order(features_i$relTIC, decreasing=TRUE),]
  compound_matched = FALSE
  for(j in 1:length(compound_table_i$species)) {
    curr_species <- compound_table_i$species[j]
    if(exists("mydb")) dbDisconnect(mydb)
    mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
    files <- dbGetQuery(mydb, paste("SELECT File_Name, Project_Name FROM UPLC_Results WHERE species_code = '",curr_species,"' AND Project_Name LIKE '%DDA%'", sep=""))
    if(nrow(files) < 1) next
    file.paths <- sapply(1:nrow(files), function(x) paste("K:/Lab_Map/Database/2015_Relational_DB/DATA_STORAGE/UPLC_MS_DATA/Data_Not_Active_Projects_etc/10_Converted_Data/",files$Project_Name[x],"/mzXML/Sample/",files$File_Name[x],".mzXML", sep=""))
    peakdata <- extract.peakdata(file.paths, MSlevel = 2)
    for(k in 1:nrow(features_i)) {
      rt_k <- features_i$rt[k]
      mz_k <- features_i$mz[k]
      current_spec <- avg.msms.spec(file.paths, peakdata, rt = rt_k, mz = mz_k)
      if(is.null(current_spec)) next
      if(nrow(current_spec) < 5) next
      msms_spec[[features_i$compound_number[k]]] <- current_spec
      compound_matched = TRUE
      break
    }
    if(compound_matched == TRUE) break
  }
}

# save msms spec files
saveRDS(msms_spec, "./results/msms_spectra_list_2017_11_15.rds")

sum(sapply(1:length(msms_spec), function(x) !is.null(msms_spec[[x]])))
compound_feature_2[compound_feature_2$compound == 4,]


