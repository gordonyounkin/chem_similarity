library(RMySQL)
source("./code/create_pairwiseComps_sampsByComps.R")
source("./code/chem_similarity_function.R")

pairwise.comps.sap <- make_pairwisecomps("K:/DDA/all_inga/sap_network_merged_spec/")
pairwise.comps.phen <- make_pairwisecomps("K:/DDA/all_inga/phen_network_merged_spec/")

sampsByCompounds <- make_sampsByCompounds("K:/DDA/all_inga/compound_tic_2017_08_02.csv", samps_to_remove = c("COJR", "Zygl"), by_species = FALSE)

sampsByCompoundsSap <- sampsByCompounds[, names(sampsByCompounds) %in% names(pairwise.comps.sap)]
sampsByCompoundsPhen <- sampsByCompounds[, names(sampsByCompounds) %in% names(pairwise.comps.phen)]

sampsCompsStandSap <- standardizeByRow(sampsByCompoundsSap)
sampsCompsStandPhen <- standardizeByRow(sampsByCompoundsPhen)

chem_similarity_sap <- chemical_similarity(sampsCompsStandSap, pairwise.comps.sap)
chem_similarity_phen <- chemical_similarity(sampsCompsStandPhen, pairwise.comps.phen)

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

pairwise.phen.percent <- outer(comp.class.pcts$phen.final.pct, comp.class.pcts$phen.final.pct, FUN = function(X,Y) (X+Y)/2)
pairwise.phen.1mindiff <- outer(comp.class.pcts$phen.final.pct, comp.class.pcts$phen.final.pct, FUN = function(X,Y) 1-abs(X-Y))
pairwise.sap.percent <- outer(comp.class.pcts$sap.final.pct, comp.class.pcts$sap.final.pct, FUN = function(X,Y) (X+Y)/2)
pairwise.sap.1mindiff <- outer(comp.class.pcts$sap.final.pct, comp.class.pcts$sap.final.pct, FUN = function(X,Y) 1-abs(X-Y))
pairwise.tyr.percent <- outer(comp.class.pcts$tyr.final.pct, comp.class.pcts$tyr.final.pct, FUN = function(X,Y) (X+Y)/2)
pairwise.tyr.1mindiff <- outer(comp.class.pcts$tyr.final.pct, comp.class.pcts$tyr.final.pct, FUN = function(X,Y) 1-abs(X-Y))

pairwise.spp <- chem_similarity_phen*pairwise.phen.percent*pairwise.phen.1mindiff + chem_similarity_sap*pairwise.sap.percent*pairwise.sap.1mindiff + pairwise.tyr.percent*pairwise.tyr.1mindiff

for(i in 1:nrow(pairwise.spp)) {
  dbDisconnect(mydb)
  mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
  species_name <- dbGetQuery(mydb, paste("SELECT Species_name from Species WHERE species_code = '", names(pairwise.spp)[i], "'", sep = ""))
  names(pairwise.spp)[i] <- paste(names(pairwise.spp)[i], species_name, sep = "_")
  row.names(pairwise.spp)[i] <- paste(row.names(pairwise.spp)[i], species_name, sep = "_")
}

write.csv(pairwise.spp, "K:/DDA/all_inga/2017_10_10_chem_similarity_with_tyrosine.csv", row.names = TRUE)



