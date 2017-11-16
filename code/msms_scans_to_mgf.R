# load pc_id_spec data and actual spectra created by merge_msms_scans_for_each_pcid.R
msms_spec <- readRDS("./results/msms_spectra_list_2017_11_15.rds")
comps_with_spec <- which(sapply(1:length(msms_spec), function(x) !is.null(msms_spec[[x]])))

# load file with mz/rt data for each compound
compound_feature <- read.csv("./data/compound_feature_table.csv")
compound_max_feature <- aggregate(compound_feature$TIC, by=list(Category=compound_feature$compound_number),FUN=max)
names(compound_max_feature) <- c("compound","TIC_max")
compound_max_feature_2 <- merge(compound_max_feature, compound_feature, by.x=c("compound","TIC_max"), by.y=c("compound_number","TIC"))
compound_max_feature_3 <- compound_max_feature_2[compound_max_feature_2$compound %in% comps_with_spec, ]


# choose compounds you are interested in
all_compounds <- compound_max_feature_3
sap_compounds <- compound_max_feature_3[compound_max_feature_3$mz >=580 & compound_max_feature_3$rt >= 18,]
phen_compounds <- compound_max_feature_3[!compound_max_feature_3$compound %in% sap_compounds$compound, ]

# turn list of compounds into mgf file
# be sure to name mgf and put it in folder you want
compounds_to_use <- phen_compounds
write_txt_file <- file("./results/phen_spectra_2017_11_16.mgf")
writeLines(
  unlist(lapply(1:nrow(compounds_to_use), function(i)
    c("BEGIN IONS",paste("PEPMASS=",compounds_to_use$mz[i],sep=""),"CHARGE=1-",paste("SCANS=",compounds_to_use$compound[i],sep=""),sapply(1:nrow(msms_spec[[compounds_to_use$compound[i]]]), function(x) paste(msms_spec[[compounds_to_use$compound[i]]][x,1],msms_spec[[compounds_to_use$compound[i]]][x,2],sep="\t")),"END IONS")))
  , write_txt_file)
close(write_txt_file)
