library(RMySQL)
source("./code/create_pairwiseComps_sampsByComps.R")
source("./code/chem_similarity_function.R")

pairwise.comps.all <- make_pairwisecomps("K:/DDA/all_inga/all_compounds_merged_spectra/", through_network = FALSE)
pairwise.comps.through <- make_pairwisecomps("K:/DDA/all_inga/all_compounds_merged_spectra/", through_network = TRUE)

sampsByCompounds <- make_sampsByCompounds("K:/DDA/all_inga/compound_tic_2017_08_02.csv", samps_to_remove = c("COJR", "Zygl"), by_species = TRUE)

sampsByCompounds <- sampsByCompounds[, names(sampsByCompounds) %in% names(pairwise.comps.through)]

sampsCompsStand <- standardizeByRow(sampsByCompounds)

# get percent extracted (for method 20) from database
dbDisconnect(mydb)
mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')

extr.pct <- dbGetQuery(mydb, "SELECT AVG(Percent_Extracted) as percent_extracted, species_code FROM (
SELECT extraction_weight.Extraction_Number, extraction_weight.Percent_Extracted, Extraction.species_code FROM `extraction_weight`
  LEFT JOIN Extraction ON extraction_weight.Extraction_Number = Extraction.Extraction_Number) a 
  GROUP BY a.species_code")

# get percent tyrosine from database
tyr.pct <- dbGetQuery(mydb, "SELECT AVG(percent_tyrosine) as percent_tyrosine, species_code FROM (
SELECT Tyrosine.Extraction_Number, SUBSTRING_INDEX(Tyrosine.Percent_Tyrosine, '%', 1)*0.01 as percent_tyrosine, Extraction.species_code FROM `Tyrosine`
                      LEFT JOIN Extraction ON Tyrosine.Extraction_Number = Extraction.Extraction_Number
                      WHERE Percent_Tyrosine != 'ND') a 
                      GROUP BY a.species_code")

# calculate percent each species has in phenolics/saponins/tryosine
comp.class.pcts <- data.frame("species_code" = row.names(sampsByCompounds), "phenSumTIC" = sapply(1:nrow(sampsByCompoundsPhen), function(x) sum(sampsByCompoundsPhen[x,])), "sapSumTIC" = sapply(1:nrow(sampsByCompoundsSap), function(x) sum(sampsByCompoundsSap[x,])))
comp.class.pcts$phenTICpct <- comp.class.pcts$phenSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)
comp.class.pcts$sapTICpct <- comp.class.pcts$sapSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)
comp.class.pcts <- merge(comp.class.pcts, extr.pct, by = "species_code", all.x = TRUE, all.y = FALSE)
comp.class.pcts <- merge(comp.class.pcts, tyr.pct, by = "species_code", all.x = TRUE, all.y = FALSE)
comp.class.pcts$percent_tyrosine[is.na(comp.class.pcts$percent_tyrosine)] <- 0
comp.class.pcts$tyr.final.pct <- comp.class.pcts$percent_tyrosine / (comp.class.pcts$percent_extracted * 0.01 + comp.class.pcts$percent_tyrosine)
comp.class.pcts$phensap.final.pct <- comp.class.pcts$percent_extracted * 0.01 / (comp.class.pcts$percent_extracted * 0.01 + comp.class.pcts$percent_tyrosine)
comp.class.pcts$phen.final.pct <- comp.class.pcts$phenTICpct * comp.class.pcts$phensap.final.pct
comp.class.pcts$sap.final.pct <- comp.class.pcts$sapTICpct * comp.class.pcts$phensap.final.pct

# add tyrosine data into sedio's dataset (tyrosine is compound 358)
sampsCompsStand$"358" <- 0
for(i in 1:nrow(comp.class.pcts)) {
  current_row <- which(row.names(sampsCompsStand) == comp.class.pcts$species_code[i])
  sampsCompsStand[current_row, ] <- sampsCompsStand[current_row, ] * comp.class.pcts$phensap.final.pct[i] / sum(sampsCompsStand[current_row, ])
  sampsCompsStand[current_row, "358"] <- comp.class.pcts$tyr.final.pct[i]
} 

# Sedio's calculations
nspp = nrow(sampsByCompounds)
ncomps = ncol(sampsByCompounds)
pairwise.spp = as.data.frame(matrix(0,nrow = nspp, ncol = nspp))
names(pairwise.spp) = row.names(sampsByCompounds)
row.names(pairwise.spp) = row.names(sampsByCompounds)

pairwise.comps = pairwise.comps.through
diags = pairwise.spp
for (k in 1:nspp){
  sppX = as.character(row.names(sampsCompsStand)[k])
  cat("Comparing ", sppX, " to itself", "\n", sep = "")
  species_comps = as.numeric(sampsCompsStand[k,]) > 0
  diags[k,k] = sum(((outer(as.numeric(sampsCompsStand[k,][species_comps]), as.numeric(sampsCompsStand[k,][species_comps])))*pairwise.comps[species_comps, species_comps]),na.rm = T)
}


for (i in 1:nspp)
  #for (i in 1:3)
{
  spp1 = as.character(row.names(sampsCompsStand)[i])
  i_comps = sampsCompsStand[i,] > 0
  for (j in i:nspp)
    #for (j in 1:3)
  {
    spp2 = as.character(row.names(sampsCompsStand)[j])
    j_comps = sampsCompsStand[j,] > 0
    cat("Comparing ", spp1, " to ", spp2, "\n", sep = "")
    pairwise.spp[i,j] = pairwise.spp[j,i] = sum(((outer(as.numeric(sampsCompsStand[i,][i_comps]), as.numeric(sampsCompsStand[j,][j_comps])))*pairwise.comps[i_comps,j_comps]), na.rm = T)/max(diags[i,i], diags[j,j]) #JW
  }
}

mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
for(i in 1:nrow(pairwise.spp)) {
  species_name <- dbGetQuery(mydb, paste("SELECT Species_name from Species WHERE species_code = '", names(pairwise.spp)[i], "'", sep = ""))
  names(pairwise.spp)[i] <- paste(names(pairwise.spp)[i], species_name, sep = "_")
  row.names(pairwise.spp)[i] <- paste(row.names(pairwise.spp)[i], species_name, sep = "_")
}

write.csv(pairwise.spp, "K:/DDA/all_inga/pairwise.spp_SEDIOMETHOD_tru_network_2017_10_24.csv")

####
library(vegan)
library(pvclust)

species_dist <- vegdist(pairwise.spp,"euclidean")  				# Distance matrix
species_single <- hclust(species_dist, method="single")					# Single linkage clustering
species_ward <- hclust(species_dist, method="ward.D")					# Ward clustering
species_complete <- hclust(species_dist, method="complete")				# Complete linkage clustering
species_centroid <- hclust(species_dist, method="centroid")				# Centroid clustering
species_median <- hclust(species_dist, method="median")				        # Median clustering
#Comparison between the distance matrix and binary matrices representing partitions 
coph1 <- cophenetic(species_single)							# Compute Patristic distances		
coph2 <- cophenetic(species_ward)
coph3 <- cophenetic(species_complete)
coph4 <- cophenetic(species_centroid)
coph5 <- cophenetic(species_median)
a <-cor(coph1, species_dist)								#Cophenetic correlations
b <-cor(coph2, species_dist)
c <-cor(coph3, species_dist)
d <-cor(coph4, species_dist)
e <-cor(coph5, species_dist)
method <- c("single","ward","complete", "centroid", "median")
mhc<- method[which.max(c(a,b,c,d,e))]
result_species <- pvclust(pairwise.spp, method.hclust=mhc, method.dist="correlation", use.cor="pairwise.complete.obs", nboot=1000,parallel=T,)


dev.new()
plot(result_species, cex = 0.1)
dev.copy2pdf(file = "K:/DDA/all_inga/2017_10_25_new_chem_dend.pdf", width = 11, height = 2)
dev.off()

dev.new()
jpeg(filename = "K:/DDA/all_inga/2017_10_25_new_chem_dend.jpg", width = 20000, height = 4000, res = 400)
plot(result_species)
dev.off()


# compare our method with Sedio's
pairwise.spp.sedio <- read.csv("K:/DDA/all_inga/pairwise.spp_SEDIOMETHOD_2017_10_12.csv")
pairwise.spp.new <- read.csv("K:/DDA/all_inga/2017_10_10_chem_similarity_with_tyrosine.csv")
pairwise.spp.sedio.thru <- read.csv("K:/DDA/all_inga/pairwise.spp_SEDIOMETHOD_tru_network_2017_10_24.csv")
row.names(pairwise.spp.sedio) <- pairwise.spp.sedio$X
pairwise.spp.sedio <- pairwise.spp.sedio[,names(pairwise.spp.sedio) != "X"]
row.names(pairwise.spp.new) <- pairwise.spp.new$X
pairwise.spp.new <- pairwise.spp.new[,names(pairwise.spp.new) != "X"]
row.names(pairwise.spp.sedio.thru) <- pairwise.spp.sedio.thru$X
pairwise.spp.sedio.thru <- pairwise.spp.sedio.thru[,names(pairwise.spp.sedio.thru) != "X"]
pairwise.spp.sedio.thru[1:5,1:5]

method.comparison <- data.frame("species" = character(0), "new" = numeric(0), "sedio" = numeric(0),"sedio_not_thru" = numeric(0), "ncomps_a" = numeric(0), "ncomps_b" = numeric(0), "diag_a" = numeric(0), "diag_b" = numeric(0), stringsAsFactors = FALSE)
current_row <- 0
ncomps <- sapply(1:nrow(sampsByCompounds), function(x) sum(as.numeric(sampsByCompounds[x,]) > 0))
for(i in 1:(nrow(pairwise.spp.sedio)-1)) {
  for(j in (i+1):nrow(pairwise.spp.sedio)) {
   current_row = current_row + 1 
   method.comparison[current_row, "species"] <- paste(names(pairwise.spp.sedio)[i], names(pairwise.spp.sedio)[j], sep = "_")
   method.comparison[current_row, "new"] <- pairwise.spp.new[i,j]
   method.comparison[current_row, "sedio"] <- pairwise.spp.sedio.thru[i,j]
   method.comparison[current_row, "sedio_not_thru"] <- pairwise.spp.sedio[i,j]
   method.comparison[current_row, "ncomps_a"] <- ncomps[i]
   method.comparison[current_row, "ncomps_b"] <- ncomps[j]
   method.comparison[current_row, "diag_a"] <- diags[i,i]
   method.comparison[current_row, "diag_b"] <- diags[j,j] 
  }
}

plot(method.comparison$sedio ~ method.comparison$sedio_not_thru)

plot(method.comparison$new ~ method.comparison$sedio, ylab = "new method", xlab = "sedio method")
lm.new.sedio <- lm(method.comparison$new ~ method.comparison$sedio)
abline(lm.new.sedio)
summary(lm.new.sedio)

hist(method.comparison$new)
hist(method.comparison$sedio)

method.comparison$difference <- method.comparison$new - method.comparison$sedio
hist(method.comparison$difference)
method.comparison[method.comparison$difference >= 0.33,]
method.comparison[method.comparison$difference <= -0.29, ]
method.comparison[abs(method.comparison$difference) <= 0.001, ]
method.comparison[abs(method.comparison$new - 0.2) <=0.01 & abs(method.comparison$sedio - 0.4) <=0.01,]

hist(as.numeric(pairwise.spp.new[3,]))
hist(as.numeric(pairwise.spp.sedio[3,]))

# closely related
plot(as.numeric(pairwise.spp.new$N42_thibaudiana) ~ as.numeric(pairwise.spp.new$LA4_thibaudiana))
lm1 <- lm(as.numeric(pairwise.spp.new$N42_thibaudiana) ~ as.numeric(pairwise.spp.new$LA4_thibaudiana))
summary(lm1)

# very far
plot(as.numeric(pairwise.spp.new$N42_thibaudiana) ~ as.numeric(pairwise.spp.new$N22_acreana))

# middle distance
plot(as.numeric(pairwise.spp.new$N42_thibaudiana) ~ as.numeric(pairwise.spp.new$LA60_laurina))

# middle distance
plot(as.numeric(pairwise.spp.new$N42_thibaudiana) ~ as.numeric(pairwise.spp.new$N35_paraensis))

# does number of compounds influence similarity score?
plot(c(method.comparison$new, method.comparison$new) ~ c(method.comparison$ncomps_a, method.comparison$ncomps_b))
plot(c(method.comparison$sedio, method.comparison$sedio) ~ c(method.comparison$ncomps_a, method.comparison$ncomps_b))

ncomps_score <- data.frame("new" = c(method.comparison$new, method.comparison$new), "diag" = c(method.comparison$diag_a, method.comparison$diag_b), "sedio" = c(method.comparison$sedio, method.comparison$sedio), "ncomps" = c(method.comparison$ncomps_a, method.comparison$ncomps_b))
ncomps_score_max <- data.frame("ncomps" = unique(ncomps_score$ncomps), "max_new" = sapply(1:length(unique(ncomps_score$ncomps)), function(x) max(ncomps_score[ncomps_score$ncomps == unique(ncomps_score$ncomps)[x],"new"])), "max_sedio" = sapply(1:length(unique(ncomps_score$ncomps)), function(x) max(ncomps_score[ncomps_score$ncomps == unique(ncomps_score$ncomps)[x],"sedio"])), "mean_new" = sapply(1:length(unique(ncomps_score$ncomps)), function(x) mean(ncomps_score[ncomps_score$ncomps == unique(ncomps_score$ncomps)[x],"new"])), "mean_sedio" = sapply(1:length(unique(ncomps_score$ncomps)), function(x) mean(ncomps_score[ncomps_score$ncomps == unique(ncomps_score$ncomps)[x],"sedio"])))
diag_max <- data.frame("diag" = sapply(1:nrow(diags), function(x) diags[x,x]), "max_new" = sapply(1:length(unique(ncomps_score$diag)), function(x) max(ncomps_score[ncomps_score$diag == unique(ncomps_score$diag)[x],"new"])), "max_sedio" = sapply(1:length(unique(ncomps_score$diag)), function(x) max(ncomps_score[ncomps_score$diag == unique(ncomps_score$diag)[x],"sedio"])), "mean_new" = sapply(1:length(unique(ncomps_score$diag)), function(x) mean(ncomps_score[ncomps_score$diag == unique(ncomps_score$diag)[x],"new"])), "mean_sedio" = sapply(1:length(unique(ncomps_score$diag)), function(x) mean(ncomps_score[ncomps_score$diag == unique(ncomps_score$diag)[x],"sedio"])))

# effect of number of compounds
plot(ncomps_score_max$max_new ~ ncomps_score_max$ncomps)
plot(ncomps_score_max$max_sedio ~ ncomps_score_max$ncomps)

plot(ncomps_score_max$mean_new ~ ncomps_score_max$ncomps, ylim = c(0,0.35), main = "new method", ylab = "mean similarity score", xlab = "number of compounds")
lm.mean.new <- lm(ncomps_score_max$mean_new ~ ncomps_score_max$ncomps)
abline(lm.mean.new)
summary(lm.mean.new)
plot(ncomps_score_max$mean_sedio ~ ncomps_score_max$ncomps, main = "new method", ylab = "mean similarity score", xlab = "number of compounds")
lm.mean.sedio <- lm(ncomps_score_max$mean_sedio ~ ncomps_score_max$ncomps)
abline(lm.mean.sedio)
summary(lm.mean.sedio)

plot(ncomps_score_max$mean_new - ncomps_score_max$mean_sedio ~ ncomps_score_max$ncomps)

# what about sum of number of compounds of two species being compared
method.comparison$ncomps_ab <- method.comparison$ncomps_a + method.comparison$ncomps_b
plot(method.comparison$new ~ method.comparison$ncomps_ab)
plot(method.comparison$sedio ~ method.comparison$ncomps_ab)

# what about difference of number of compounds of two species being compared
method.comparison$ncomps_aminb <- abs(method.comparison$ncomps_a - method.comparison$ncomps_b)
plot(method.comparison$new ~ method.comparison$ncomps_aminb)
plot(method.comparison$sedio ~ method.comparison$ncomps_aminb)

# what about diagonal from sedio's method (Rao's quadratic entropy)
# max score
plot(diag_max$max_new ~ diag_max$diag)
plot(diag_max$max_sedio ~ diag_max$diag)
# mean score
par(mar = c(4,4,2,0.5))
summary(lm.mean.new)
plot(diag_max$mean_new ~ diag_max$diag, main = "new method", xlab = "Rao's quadratic entropy", ylab = "average similarity score")
lm.mean.new <- lm(diag_max$mean_new ~ diag_max$diag)
abline(lm.mean.new)
plot(diag_max$mean_sedio ~ diag_max$diag, main = "sedio method", xlab = "Rao's quadratic entropy", ylab = "average similarity score")
lm.mean.sedio <- lm(diag_max$mean_sedio ~ diag_max$diag)
abline(lm.mean.sedio)
summary(lm.mean.sedio)
plot(diag_max$mean_new - diag_max$mean_sedio ~ diag_max$diag)
