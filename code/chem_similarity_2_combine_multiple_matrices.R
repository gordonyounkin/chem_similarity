library(RMySQL)
library(plyr)
source("./code/create_pairwiseComps_sampsByComps.R")
source("./code/chem_similarity_function.R")

# load pairwise.comps/sampsByCompounds files
#    pairwise.comps.sap <- make_pairwisecomps("K:/DDA/all_inga/sap_network_merged_spec/")
#    pairwise.comps.phen <- make_pairwisecomps("K:/DDA/all_inga/phen_network_merged_spec/")

#    sampsByCompounds <- make_sampsByCompounds("./data/compound_tic_2017_08_02.csv", samps_to_remove = c("COJR", "Zygl"), by_species = FALSE)

#    ampsByCompoundsSap <- sampsByCompounds[,names(sampsByCompounds) %in% names(pairwise.comps.sap)]
#    sampsByCompoundsPhen <- sampsByCompounds[,names(sampsByCompounds) %in% names(pairwise.comps.phen)]

# load any similarity files you want to combine
chem_similarity_phen <- read.csv("./results/phen_sim_filled_comps2_LOG_2017_11_20.csv")
row.names(chem_similarity_phen) <- chem_similarity_phen$X
chem_similarity_phen <- chem_similarity_phen[,names(chem_similarity_phen) != "X"]
chem_similarity_sap <- read.csv("./results/sap_sim_filled_comps2_LOG_2017_11_20.csv")
row.names(chem_similarity_sap) <- chem_similarity_sap$X
chem_similarity_sap <- chem_similarity_sap[,names(chem_similarity_sap) != "X"]

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
                      WHERE Percent_Tyrosine != 'ND'
AND LCASE(NOTES) NOT LIKE '%drop%') a 
                      GROUP BY a.species_code")

# calculate percent each species has in phenolics/saponins/tryosine
comp.class.pcts <- data.frame("sample" = row.names(sampsByCompounds), 
                              "phenSumTIC" = sapply(1:nrow(sampsByCompoundsPhen), function(x) sum(sampsByCompoundsPhen[x,])), 
                              "sapSumTIC" = sapply(1:nrow(sampsByCompoundsSap), function(x) sum(sampsByCompoundsSap[x,])), 
                              "species_code" = sapply(1:nrow(sampsByCompoundsPhen), 
                                                      function(x) unlist(strsplit(row.names(sampsByCompoundsPhen)[x], split = "_"))[1]),
                              stringsAsFactors=FALSE)
comp.class.pcts$logPhen <- log(comp.class.pcts$phenSumTIC)
comp.class.pcts$logSap <- log(comp.class.pcts$sapSumTIC)
comp.class.pcts$phenTICpct <- comp.class.pcts$phenSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)
comp.class.pcts$sapTICpct <- comp.class.pcts$sapSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)
comp.class.pcts <- join(comp.class.pcts, extr.pct, by="species_code", type="left", match="all")
comp.class.pcts <- join(comp.class.pcts, tyr.pct, by="species_code", type="left", match="all")
comp.class.pcts$percent_tyrosine[is.na(comp.class.pcts$percent_tyrosine)] <- 0
comp.class.pcts$tyr.final.pct <- comp.class.pcts$percent_tyrosine / (comp.class.pcts$percent_extracted * 0.01 + comp.class.pcts$percent_tyrosine)
comp.class.pcts$phensap.final.pct <- comp.class.pcts$percent_extracted * 0.01 / (comp.class.pcts$percent_extracted * 0.01 + comp.class.pcts$percent_tyrosine)
comp.class.pcts$phen.final.pct <- comp.class.pcts$phenTICpct * comp.class.pcts$phensap.final.pct
comp.class.pcts$sap.final.pct <- comp.class.pcts$sapTICpct * comp.class.pcts$phensap.final.pct


write.csv(comp.class.pcts, "./results/phen_sap_tyr_sample_percents_2017_11_20.csv")

pairwise.phen.percent <- outer(comp.class.pcts$phen.final.pct, comp.class.pcts$phen.final.pct, FUN = function(X,Y) (X+Y)/2)
pairwise.phen.1mindiff <- outer(comp.class.pcts$phen.final.pct, comp.class.pcts$phen.final.pct, FUN = function(X,Y) 1-abs(X-Y))
pairwise.phen.min <- outer(comp.class.pcts$phen.final.pct, comp.class.pcts$phen.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.sap.percent <- outer(comp.class.pcts$sap.final.pct, comp.class.pcts$sap.final.pct, FUN = function(X,Y) (X+Y)/2)
pairwise.sap.1mindiff <- outer(comp.class.pcts$sap.final.pct, comp.class.pcts$sap.final.pct, FUN = function(X,Y) 1-abs(X-Y))
pairwise.sap.min <- outer(comp.class.pcts$sap.final.pct, comp.class.pcts$sap.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.tyr.percent <- outer(comp.class.pcts$tyr.final.pct, comp.class.pcts$tyr.final.pct, FUN = function(X,Y) (X+Y)/2)
pairwise.tyr.1mindiff <- outer(comp.class.pcts$tyr.final.pct, comp.class.pcts$tyr.final.pct, FUN = function(X,Y) 1-abs(X-Y))
pairwise.tyr.min <- outer(comp.class.pcts$tyr.final.pct, comp.class.pcts$tyr.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))

pairwise.spp <- chem_similarity_phen*pairwise.phen.percent*pairwise.phen.1mindiff + chem_similarity_sap*pairwise.sap.percent*pairwise.sap.1mindiff + pairwise.tyr.percent*pairwise.tyr.1mindiff


# try not doing 1-difference when combining compound classes...it seems to be splitting some species
#    pairwise.spp.mincompclass <- chem_similarity_phen*pairwise.phen.min + chem_similarity_sap*pairwise.sap.min + pairwise.tyr.min
#    write.csv(pairwise.spp.mincompclass, "./data/BCI_test/BCI_sim_matrix_MINCOMP.csv")

# try just using average of investment in each compound class
#    pairwise.spp.avgcompclass <- chem_similarity_phen*pairwise.phen.percent + chem_similarity_sap*pairwise.sap.percent + pairwise.tyr.percent
#    write.csv(pairwise.spp.avgcompclass, "K:/GY_LAB_FILES/github_repositories/chem_similarity/results/2017_10_31_pairwise.samps.avgcompclass.csv")


dbDisconnect(mydb)
mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
for(i in 1:nrow(pairwise.spp)) {
  species_name <- dbGetQuery(mydb, paste("SELECT Species_name from Species WHERE species_code = '", unlist(strsplit(row.names(pairwise.spp)[i], split = "_"))[1], "'", sep = ""))
  names(pairwise.spp)[i] <- paste(names(pairwise.spp)[i], species_name, sep = "_")
  row.names(pairwise.spp)[i] <- paste(row.names(pairwise.spp)[i], species_name, sep = "_")
}

write.csv(pairwise.spp, "./results/combined_similarity_matrix_2017_11_20.csv")


# create chem similarity tree from similarity matrix (similarity matrix should be named 'pairwise.spp')
library(vegan)
library(pvclust)

species_dist <- vegdist(pairwise.spp,"euclidean")      		# Distance matrix
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
method <- c("single","ward.D","complete", "centroid", "median")
mhc<- method[which.max(c(a,b,c,d,e))]
result_samples <- pvclust(pairwise.spp, method.hclust=mhc, method.dist="correlation", use.cor="pairwise.complete.obs", nboot=1000,parallel=T)

# save tree--make sure to give it a name
dev.new()
plot(result_samples, cex=1.66, cex.pv=1, lwd=1, float = 0.003)
dev.copy2pdf(file = "./results/chem_dendrogram_all_samps_2017_11_20.pdf", width = 200, height = 20)
dev.off()

