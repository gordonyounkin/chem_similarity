library(foreach)
library(doParallel)
source("./code/create_pairwiseComps_sampsByComps.R")
source("./code/chem_similarity_function.R")

pairwise.comps.phen <- make_pairwisecomps("./data/phenolics_network_2017_11_16/")
pairwise.comps.sap <- make_pairwisecomps("./data/saponin_network_2017_11_16/")

# load file with compounds that are contaminants and remove:
contaminants <- read.csv("./data/contaminants_2017_11_16.csv")
non_contaminants <- names(pairwise.comps.phen)[!names(pairwise.comps.phen) %in% contaminants$compound_number]
pairwise.comps.phen <- pairwise.comps.phen[non_contaminants, non_contaminants]
non_contaminants_sap <- names(pairwise.comps.sap)[!names(pairwise.comps.sap) %in% contaminants$compound_number]
pairwise.comps.sap <- pairwise.comps.sap[non_contaminants_sap, non_contaminants_sap]

sampsByCompounds <- make_sampsByCompounds("./data/filled_compound_table_major_features_2017_11_30.csv", by_species = FALSE)
sampsByCompoundsLog <- log(sampsByCompounds)
sampsByCompoundsLog[sampsByCompoundsLog <= 0] <- 0
write.csv(t(sampsByCompounds), "./data/filled_samps_by_comps_2017_11_30.csv")

# without log TIC
sampsByCompoundsSap <- sampsByCompounds[,names(sampsByCompounds) %in% names(pairwise.comps.sap)]
sampsByCompoundsPhen <- sampsByCompounds[,names(sampsByCompounds) %in% names(pairwise.comps.phen)]

# with log TIC
sampsByCompoundsSapLOG <- sampsByCompoundsLog[,names(sampsByCompoundsLog) %in% names(pairwise.comps.sap)]
sampsByCompoundsPhenLOG <- sampsByCompoundsLog[,names(sampsByCompoundsLog) %in% names(pairwise.comps.phen)]

sampsCompsStandSap <- standardizeByRow(sampsByCompoundsSapLOG)
sampsCompsStandPhen <- standardizeByRow(sampsByCompoundsPhenLOG)

# PARALLELIZE CHEM SIMILARITY CODE
cores = detectCores()
cl <- makeCluster(cores[1]) #not to overload your computer
registerDoParallel(cl)

# Saponins
similarity_matrix_sap <- foreach(i = 1:nrow(sampsCompsStandSap), .combine = rbind) %:% foreach(j = 1:nrow(sampsCompsStandSap)) %dopar% {
  if(i <= j) {
  chemical_similarity_single(row.names(sampsCompsStandSap)[i], row.names(sampsCompsStandSap)[j], sampsCompsStandSap, pairwise.comps.sap)
  }
  else NA
}
simmatrixbackup <- similarity_matrix_sap
similarity_matrix_sap <- as.data.frame(similarity_matrix_sap)
names(similarity_matrix_sap) <- row.names(sampsCompsStandSap)
row.names(similarity_matrix_sap) <- row.names(sampsCompsStandSap)
for(i in 1:ncol(similarity_matrix_sap)) similarity_matrix_sap[,i] <- unlist(similarity_matrix_sap[,i])
for(i in 1:nrow(similarity_matrix_sap)) {
  for(j in i:nrow(similarity_matrix_sap)) {
    similarity_matrix_sap[j,i] <- similarity_matrix_sap[i,j]
  }
}

write.csv(similarity_matrix_sap, "./data/BCI_test/BCI_sap_sim_filled_comps_LOG_2017_11_19.csv")


# Phenolics
similarity_matrix_phen <- foreach(i = 1:nrow(sampsCompsStandPhen), .combine = rbind) %:% foreach(j = 1:nrow(sampsCompsStandPhen)) %dopar% {
  if(i <= j) {
    chemical_similarity_single(row.names(sampsCompsStandPhen)[i], row.names(sampsCompsStandPhen)[j], sampsCompsStandPhen, pairwise.comps.phen)
  }
  else NA
}
similarity_matrix_phen <- as.data.frame(similarity_matrix_phen)
names(similarity_matrix_phen) <- row.names(sampsCompsStandPhen)
row.names(similarity_matrix_phen) <- row.names(sampsCompsStandPhen)
for(i in 1:ncol(similarity_matrix_phen)) similarity_matrix_phen[,i] <- unlist(similarity_matrix_phen[,i])
for(i in 1:nrow(similarity_matrix_phen)) {
  for(j in i:nrow(similarity_matrix_phen)) {
    similarity_matrix_phen[j,i] <- similarity_matrix_phen[i,j]
  }
}

write.csv(similarity_matrix_phen, "./data/BCI_test/BCI_phen_sim_filled_comps_LOG_2017_11_19.csv")

