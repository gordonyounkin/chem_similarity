library(foreach)
source("./Code/create_pairwiseComps_sampsByComps.R")
source("./Code/chem_similarity_function.R")

pairwise.comps <- make_pairwisecomps("./data/all_compounds_merged_spectra/all_compounds_merged_spectra")


sampsByCompounds <- make_sampsByCompounds("./data/rand.samples.csv", samps_to_remove = c("COJR", "Zygl"), by_species = TRUE)
sampsByCompounds1 <- sampsByCompounds[,names(sampsByCompounds) %in% names(pairwise.comps)]
sampsCompsStand <- standardizeByRow(sampsByCompounds1)

# PARALLELIZE CHEM SIMILARITY CODE
cores = detectCores()
cl <- makeCluster(cores[1]) #not to overload your computer

registerDoParallel(cl)

similarity_matrix <- foreach(i = 1:nrow(sampsCompsStand), .combine = rbind) %:% foreach(j = 1:nrow(sampsCompsStand)) %dopar% {
  if(i <= j) {
  chemical_similarity_single(row.names(sampsCompsStand)[i], row.names(sampsCompsStand)[j], sampsCompsStand, pairwise.comps)
  }
  else NA
}

similarity_matrix <- as.data.frame(similarity_matrix)
names(similarity_matrix) <- row.names(sampsCompsStand)
row.names(similarity_matrix) <- row.names(sampsCompsStand)

for(i in 1:nrow(similarity_matrix)) {
  for(j in i:nrow(similarity_matrix)) {
    similarity_matrix[j,i] <- similarity_matrix[i,j]
  }
}

for(i in 1:ncol(similarity_matrix)) similarity_matrix[,i] <- unlist(similarity_matrix[,i])

chem_sim_matrix_162spec <- similarity_matrix[names(Inga_phy_dist_match_chem),names(Inga_phy_dist_match_chem)]
write.csv(similarity_matrix, "./data/similarity_randspecies.csv")

