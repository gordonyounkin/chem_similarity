library(RMySQL)
library(splitstackshape)
source("./code/average_multiple_msms_scans.R")

compound_pcid <- read.csv("K:/DDA/all_inga/compound_pcid_2017_08_01.csv")
pcid <- read.csv("K:/DDA/all_inga/all_pc_id_2017_08_01.csv")
pc_id_spec2 <- merge(pcid, compound_pcid, by = "PC_ID", all.x = TRUE)
head(pc_id_spec2,30)

# select highest percent TIC for each PC_ID
pc_ids <- unique(pc_id_spec2$PC_ID)
pc_id_spec <- pc_id_spec2[-(1:nrow(pc_id_spec2)), ]
for(i in 1:length(pc_ids)) {
  pc_id_i <- pc_id_spec2[pc_id_spec2$PC_ID == pc_ids[i], ]
  pc_id_spec <- rbind(pc_id_spec, pc_id_i[which.max(pc_id_i$"Percent_TIC"), ])
}
pc_id_spec$MZ_RT <- as.character(pc_id_spec$MZ_RT)
pc_id_spec$MS_MS_Spec_ID <- as.character(pc_id_spec$MS_MS_Spec_ID)

# if available, get spectrum for each pc_id
species <- unique(pc_id_spec$species)
msms_spec <- vector("list", nrow(pc_id_spec))
for(i in 1:length(species)) {
  if(exists("mydb")) dbDisconnect(mydb)
  mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
  files <- dbGetQuery(mydb, paste("SELECT File_Name, Project_Name FROM UPLC_Results WHERE species_code = '",species[i],"' AND Project_Name LIKE '%DDA%'", sep=""))
  if(nrow(files) < 1) next
  file.paths <- sapply(1:nrow(files), function(x) paste("K:/Lab_Map/Database/2015_Relational_DB/DATA_STORAGE/UPLC_MS_DATA/Data_Not_Active_Projects_etc/10_Converted_Data/",files$Project_Name[x],"/mzXML/Sample/",files$File_Name[x],".mzXML", sep=""))
  peakdata <- extract.peakdata(file.paths, MSlevel = 2)
  spec_i <- pc_id_spec[pc_id_spec$species == species[i], ]
  for(k in 1:nrow(spec_i)) {
    print(as.character(spec_i$PC_ID[k]))
    rt_k <- as.numeric(unlist(strsplit(spec_i$MZ_RT[k], split = "_"))[2])
    mz_k <- as.numeric(unlist(strsplit(spec_i$MZ_RT[k], split = "_"))[1])
    current_spec <- avg.msms.spec(file.paths, peakdata, rt = rt_k, mz = mz_k)
    if(is.null(current_spec)) next
    if(nrow(current_spec) < 5) next
    msms_spec[[as.numeric(row.names(spec_i)[k])]] <- current_spec
    pc_id_spec[row.names(spec_i)[k], "MS_MS_Spec_ID"] <- as.numeric(row.names(spec_i)[k])
    pc_id_spec[row.names(spec_i)[k], "MSMS_TIC"] <- sum(current_spec[,2])
  }
}

# choose spectrum with highest TIC per compound
pc_id_spec_2 <- pc_id_spec[pc_id_spec$MS_MS_Spec_ID != "Null", ]

# save msms spec files
saveRDS(pc_id_spec_2, "./results/pcid_spectra_ids_2017_09_22.rds")
saveRDS(msms_spec, "./results/msms_spectra_list_2017_09_22.rds")





