library(gdata)

# Calculates chemical similarity between samples using our newest method (as of 10/10/2017). Requires sampsCompsStand (dataframe with samples as rows and compounds are TICS. Each sample is normalized to 1) and pairwise.comps (compound x compound similarity matrix)
chemical_similarity <- function(sampsCompsStand, pairwise.comps) {
  nspp = nrow(sampsCompsStand)
  pairwise.spp = as.data.frame(matrix(0,nrow = nspp, ncol = nspp))
  names(pairwise.spp) = row.names(sampsCompsStand)
  row.names(pairwise.spp) = row.names(sampsCompsStand)
  for (i in 1:nspp)
  {
    spp1 = as.character(row.names(sampsCompsStand)[i])
    i_compounds = which(sampsCompsStand[i,] != 0)
    for (j in i:nspp)
    {
      spp2 = as.character(row.names(sampsCompsStand)[j])
      cat("Comparing ", spp1, " to ", spp2, "\n", sep = "")
      j_compounds = which(sampsCompsStand[j,] != 0)
      if(length(i_compounds) == 0 | length(j_compounds) == 0) {
        pairwise.spp[i,j] = pairwise.spp[j,i] = 0
        next
      }
      ij_compounds = sort(unique(c(i_compounds, j_compounds)))
      i_tics = sampsCompsStand[i, ij_compounds]
      j_tics = sampsCompsStand[j, ij_compounds]
      shared_tics = sapply(1:length(ij_compounds), function(x) min(i_tics[x], j_tics[x]))
      i_tics = i_tics - shared_tics
      j_tics = j_tics - shared_tics
      similarity = sum(shared_tics)
      current_compounds = data.matrix(pairwise.comps[names(i_tics),names(j_tics)])
      while(sum(i_tics) > 0.0001 & sum(j_tics) > 0.0001) {
        i_tics = i_tics[, i_tics > 0.0001, drop = FALSE]
        j_tics = j_tics[, j_tics > 0.0001, drop = FALSE]
        if(length(i_tics) < 1 | length(j_tics) < 1) break
        current_compounds = data.matrix(pairwise.comps[names(i_tics),names(j_tics)])
        current_cos = max(current_compounds)
        if(current_cos == min(current_compounds)) {
          similarity = similarity + min(sum(i_tics), sum(j_tics)) * current_cos
          break
        }
        current_comp_locations = which(current_compounds == current_cos)
        i.numbers = current_comp_locations %% nrow(current_compounds)
        i.numbers[i.numbers==0] = nrow(current_compounds)
        j.numbers = ceiling(current_comp_locations/nrow(current_compounds))
        
        i.repeats = sapply(1:length(i.numbers), function(x) sum(i.numbers == i.numbers[x]))
        j.repeats = sapply(1:length(j.numbers), function(x) sum(j.numbers == j.numbers[x]))
        ij.min = sapply(1:length(i.numbers), function(x) min(i_tics[i.numbers[x]]/i.repeats[x], j_tics[j.numbers[x]]/j.repeats[x]))
        
        similarity = similarity + sum(ij.min * current_cos)
        
        for(k in 1:length(i.numbers)) {
          i_tics[i.numbers[k]] = i_tics[i.numbers[k]] - ij.min[k]
          j_tics[j.numbers[k]] = j_tics[j.numbers[k]] - ij.min[k]
        }
        current_compounds[current_comp_locations] <- 0
      }
      pairwise.spp[i,j] = pairwise.spp[j,i] = similarity
    }
  }
  return(pairwise.spp)
}
