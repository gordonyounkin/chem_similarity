sampsCompsStand <- standardizeByRow(sampsByCompounds)
sample_1 = row.names(sampsCompsStand)[720]
sample_2 = row.names(sampsCompsStand)[346]

sampsByCompounds <- sampsCompsStand
pairwise.comps <- pairwise.comps.phen
# build function that takes sampsByCompounds + pairwise.comps and returns the x compounds that contribute most to similarity
function(sample_1, sample_2, sampsByCompounds, pairwise.comps) {
  samp_1_comps = sampsByCompounds[row.names(sampsByCompounds) == sample_1, ,drop=FALSE]
  samp_1_comps = samp_1_comps[,samp_1_comps > 0 & names(samp_1_comps) %in% names(pairwise.comps)]
  samp_2_comps = sampsByCompounds[row.names(sampsByCompounds) == sample_2, ]
  samp_2_comps = samp_2_comps[,samp_2_comps > 0 & names(samp_2_comps) %in% names(pairwise.comps)]
  pairwise.comps <- pairwise.comps[names(samp_1_comps), names(samp_2_comps)]
  shared_tic = matrix(data = NA, nrow = ncol(samp_1_comps), ncol = ncol(samp_2_comps))
  for(i in 1:ncol(samp_1_comps)) {
    for(j in 1:ncol(samp_2_comps)) {
      shared_tic[i,j] = min(samp_1_comps[1,i], samp_2_comps[1,j])
    }
  }  
  contribution = pairwise.comps * shared_tic
  compounds_of_interest <- order(contribution, decreasing = TRUE)[1:5]
  important_samp1_comps <- unique(sapply(1:length(compounds_of_interest), function(x) row.names(contribution)[compounds_of_interest[x] %% nrow(contribution)]))
  important_samp2_comps <- unique(sapply(1:length(compounds_of_interest), function(x) names(contribution)[ceiling(compounds_of_interest[x] / nrow(contribution))]))
  
}


