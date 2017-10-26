# load pc_id_spec data and actual spectra created by merge_msms_scans_for_each_pcid.R
pc_id_spec_2 <- readRDS("./results/pcid_spectra_ids_2017_09_22.rds")
msms_spec <- readRDS("./results/msms_spectra_list_2017_09_22.rds")

compound_spec <- pc_id_spec_2[-(1:nrow(pc_id_spec_2)),]
for(i in 1:length(unique(pc_id_spec_2$compound_number))) {
  compound_i <- unique(pc_id_spec_2$compound_number)[i]
  compound_spec_i <- pc_id_spec_2[pc_id_spec_2$compound_number==compound_i,]
  compound_spec[i,] <- compound_spec_i[which.max(compound_spec_i$MSMS_TIC),]
}
compound_spec <- as.data.frame(cSplit(compound_spec,'MZ_RT', sep="_", type.convert= FALSE))
compound_spec$MS_MS_Spec_ID <- as.numeric(compound_spec$MS_MS_Spec_ID)
names(compound_spec)[names(compound_spec)=="MZ_RT_1"] <- "mz"
names(compound_spec)[names(compound_spec)=="MZ_RT_2"] <- "rt"
compound_spec$mz <- as.numeric(compound_spec$mz)
compound_spec$rt <- as.numeric(compound_spec$rt)
head(compound_spec)
unique(compound_spec$compound_number)

# choose compounds you are interested in
compound_spec_all <- compound_spec
compound_spec_sap <- compound_spec[compound_spec$mz >=580 & compound_spec$rt >= 18,]
compound_spec_phen <- compound_spec[!compound_spec_all$compound_number %in% compound_spec_sap$compound_number, ]
libhits075 <- read.csv("K:/DDA/in_silico_fragmentation/library_hits_over_075.csv")
libhits095 <- libhits075[libhits075$MQScore >= 0.95,]
compound_spec_high_cosine <- compound_spec[compound_spec$compound_number %in% libhits095$compound, ]
tyrosine_gallate <- compound_spec_all[compound_spec_all$compound_number == 4978, ,drop=FALSE]

# turn list of compounds into mgf file
# be sure to name mgf and put it in folder you want
compound_spec <- tyrosine_gallate
write_txt_file <- file("K:/DDA/in_silico_fragmentation/tyrosine_gallate_real_msms.mgf")
writeLines(
  unlist(lapply(1:nrow(compound_spec), function(i)
    c("BEGIN IONS",paste("PEPMASS=",compound_spec$mz[i],sep=""),"CHARGE=1-",paste("SCANS=",compound_spec$compound_number[i],sep=""),sapply(1:nrow(msms_spec[[compound_spec$MS_MS_Spec_ID[i]]]), function(x) paste(msms_spec[[compound_spec$MS_MS_Spec_ID[i]]][x,1],msms_spec[[compound_spec$MS_MS_Spec_ID[i]]][x,2],sep="\t")),"END IONS")))
  , write_txt_file)
close(write_txt_file)