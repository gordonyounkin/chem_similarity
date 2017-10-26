library(reshape2)
library(splitstackshape)

dbDisconnect(mydb)
mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')

# set working directory to folder with clean network from GNPS ('view raw spectra')
setwd("K:/DDA/all_inga/all_inga_merged_spectra_2017_09_22/")
# load .pairsinfo file from networkedges folder
network_clean <- read.delim(file = "./networkedges/1209a7bf5c764e7bac63a8716269ba97.pairsinfo",header=FALSE, sep="\t")
# load .tsv file from current directory
compound_cluster_match <- read.delim(file = "METABOLOMICS-SNETS-31662883-view_raw_spectra-main.tsv",header=TRUE, sep="\t")[,c("cluster.index","ScanNumber")]
head(network_clean)
names(network_clean) <- c("CLUSTERID1","CLUSTERID2","DeltaMZ","V4","Cosine","OtherScore","V7")
head(compound_cluster_match)

# convert network from clusterid to scan number (same as compound number)
names(compound_cluster_match) <- c("CLUSTERID1","NODE1")
compound_network <- merge(network_clean, compound_cluster_match, by="CLUSTERID1",all.x=TRUE)
names(compound_cluster_match) <- c("CLUSTERID2","NODE2")
compound_network <- merge(compound_network, compound_cluster_match, by="CLUSTERID2",all.x=TRUE)
head(compound_network)
compound_network <- compound_network[,c("NODE1","NODE2","DeltaMZ","Cosine")]

head(compound_network)
length(unique(c(as.vector(compound_network$NODE1),as.vector(compound_network$NODE2))))
# only includes compounds that are connected to another compound in the network

# which compounds are in this network
nodes <- sort(as.numeric(unique(c(compound_network$NODE1, compound_network$NODE2))))

# create sampsByCompounds table
# load compound tic table (created in python)
compound_tic_table <- read.csv("K:/DDA/all_inga/compound_tic_2017_08_02.csv") 
sampsByCompounds <- dcast(compound_tic_table, compound_number  ~ compound_sample, value.var = "TIC")
head(compound_tic_table)

# normalize by TIC of associated RTI
samples <- row.names(sampsByCompounds)
samp_RTI <- data.frame("sample" = samples, "filename" = sapply(1:length(samples), function(x) paste(unlist(strsplit(samples[x], split = "_"))[2:length(unlist(strsplit(samples[x], split = "_")))], collapse = "_")))
samp_RTI$RTI <- character(nrow(samp_RTI))
for(i in 1:nrow(samp_RTI)) {
  samp_RTI$RTI[i] <- dbGetQuery(mydb, paste("SELECT Asoc_RTI from UPLC_Results WHERE File_Name LIKE '%",samp_RTI$filename[i],"'", sep = ""))[1,1]}

samp_RTI$rti_tic <- numeric(nrow(samp_RTI))
for(i in 1:nrow(samp_RTI)) {
  samp_RTI$rti_tic[i] <- sum(dbGetQuery(mydb, paste("SELECT TIC_into FROM RTI_QC WHERE RTI = '",samp_RTI$RTI[i],"'",sep = ""))$TIC_into)
}
samp_RTI <- samp_RTI[samp_RTI$rti_tic!=0, ]
samp_RTI$rti_multiplier <- mean(samp_RTI$rti_tic) / samp_RTI$rti_tic
# for each sample, multiply by rti multiplier calculated above.
compound_tic_table_2 <- merge(compound_tic_table, samp_RTI[c("sample","rti_multiplier")], by.x = "compound_sample", by.y = "sample", all.x = TRUE, all.y = FALSE)
compound_tic_table_2$adj_rti <- compound_tic_table_2$TIC * compound_tic_table_2$rti_multiplier

### if you want pairwise species instead of sample
compound_tic_table_species <- data.frame(compound_tic_table_2,cSplit(compound_tic_table_2, 'compound_sample', sep="_", type.convert=FALSE))
compound_tic_table_species_1 <- aggregate(compound_tic_table_species$adj_rti, by=list(Category=compound_tic_table_species$compound_number,compound_tic_table_species$compound_sample_1),FUN=mean)
names(compound_tic_table_species_1) <- c( "compound_number","species","TIC")
sampsBySpecies <- dcast(compound_tic_table_species_1, compound_number  ~ species, value.var = "TIC")
## only if you want a pairwise species matrix instead of sample matrix
sampsByCompounds <-sampsBySpecies 
##

sampsByCompounds <- as.data.frame(sampsByCompounds)
sampsByCompounds <- sampsByCompounds[sampsByCompounds$compound_number %in% nodes, ]
sampsByCompounds <- t(sampsByCompounds)
sampsByCompounds[is.na(sampsByCompounds)] <- 0
sampsByCompounds <- as.data.frame(sampsByCompounds)
row.names(sampsByCompounds)
names(sampsByCompounds) <- sampsByCompounds["compound_number",]
sampsByCompounds <- sampsByCompounds[-which(row.names(sampsByCompounds)=="compound_number"),]

tyrosine <- sampsByCompounds[, "358",drop = FALSE]
sampsByCompounds$'358'

# load mz & rt data for compounds
pc_id_spec <- read.csv("K:/DDA/all_inga/all_inga_compound_msms_spec_2017_08_02.csv")
compound_spec <- pc_id_spec[-(1:nrow(pc_id_spec)),]
for(i in 1:length(unique(pc_id_spec$compound_number))) {
  compound_i <- unique(pc_id_spec$compound_number)[i]
  compound_spec_i <- pc_id_spec[pc_id_spec$compound_number==compound_i,]
  compound_spec[i,] <- compound_spec_i[which.max(compound_spec_i$totIonCurrent),]
}
compound_spec <- as.data.frame(cSplit(compound_spec,'MZ_RT', sep="_", type.convert= FALSE))
names(compound_spec)[names(compound_spec)=="MZ_RT_1"] <- "mz"
names(compound_spec)[names(compound_spec)=="MZ_RT_2"] <- "rt"
compound_spec$mz <- as.numeric(compound_spec$mz)
compound_spec$rt <- as.numeric(compound_spec$rt)
names(compound_spec)[names(compound_spec) == "compound_number"] <- "compound"
#separate into saponins and phenolics
sap.comps <- compound_spec[compound_spec$rt > 18 & compound_spec$mz > 580, "compound"]
phen.comps <- compound_spec$compound[!compound_spec$compound %in% sap.comps]

TICsums <- data.frame("species_code" = row.names(sampsByCompounds), "TICsum" = rowSums(sampsByCompounds), "phensum" = rowSums(sampsByCompounds[,as.numeric(names(sampsByCompounds)) %in% phen.comps]), "sapsum" = rowSums(sampsByCompounds[,as.numeric(names(sampsByCompounds)) %in% sap.comps]))

## get stuff from database
dbDisconnect(mydb)
mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
ingadb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga', host='mysql.chpc.utah.edu')

tyr.pct <- dbGetQuery(mydb, "SELECT Tyrosine.Extraction_Number, Tyrosine.Percent_Tyrosine, Extraction.species_code FROM `Tyrosine`
  LEFT JOIN Extraction ON Tyrosine.Extraction_Number = Extraction.Extraction_Number
  WHERE Percent_Tyrosine != 'ND'")
tyr.pct$Percent_Tyrosine <- as.numeric(sub("%", "", tyr.pct$Percent_Tyrosine))
head(tyr.pct)
tyr.pct.avg <- aggregate(tyr.pct$Percent_Tyrosine, by=list(Category=tyr.pct$species_code),FUN=mean)
names(tyr.pct.avg) <- c("species_code", "percent_tyrosine")

extr.pct <- dbGetQuery(mydb, "SELECT extraction_weight.Extraction_Number, extraction_weight.Percent_Extracted, Extraction.species_code FROM `extraction_weight`
LEFT JOIN Extraction ON extraction_weight.Extraction_Number = Extraction.Extraction_Number")
extr.pct.avg <- aggregate(extr.pct$Percent_Extracted, by=list(Category=extr.pct$species_code),FUN=mean, na.rm = TRUE)
names(extr.pct.avg) <- c("species_code", "percent_extracted")


extr.weights <- dbGetQuery(ingadb, "SELECT extraction_weight.`Chem#`, extraction_weight.Processed_dry_weight_g as DW_g, extraction_weight.Percent_of_Volume_that_was_fractionated as percent_fractionated, extraction_weight.`ODS:_60%_MeOH_Fraction_g` as phenolics_weight,
                           `extraction_weight`.`ODS:_100%_MeOH_Fractiong` as saponin_weight, extraction_weight.Percent_recovery FROM `extraction_weight`")
chem.spcode <- dbGetQuery(mydb, "SELECT `Chem#`, Species_code FROM chemistry")

extr.weights <- merge(extr.weights, chem.spcode, by = "Chem#", all = FALSE)
head(extr.weights)
extr.weights$sap_phen_weight <- (extr.weights$phenolics_weight + extr.weights$saponin_weight) / extr.weights$DW_g
extr.weights$sap_weight <- extr.weights$saponin_weight / extr.weights$DW_g
extr.weights$phen_weight <- extr.weights$phenolics_weight / extr.weights$DW_g

extr.weights <- extr.weights[extr.weights$percent_fractionated %in% c(NA,1),]

extr.weights.avg <- aggregate(extr.weights$sap_phen_weight, by=list(Category=extr.weights$Species_code),FUN=mean, na.rm = TRUE)
names(extr.weights.avg) <- c("species_code", "phensap_pct")
phen.weights.avg <- aggregate(extr.weights$phen_weight, by=list(Category=extr.weights$Species_code),FUN=mean, na.rm = TRUE)
names(phen.weights.avg) <- c("species_code", "phen_pct")
sap.weights.avg <- aggregate(extr.weights$sap_weight, by=list(Category=extr.weights$Species_code),FUN=mean, na.rm = TRUE)
names(sap.weights.avg) <- c("species_code", "sap_pct")

ext.pcts <- merge(extr.weights.avg, phen.weights.avg, by = "species_code", all = FALSE)
ext.pcts <- merge(ext.pcts, sap.weights.avg, by = "species_code", all = FALSE)
ext.pcts <- ext.pcts[!is.na(ext.pcts$phensap_pct), ]

ext.pcts <- merge(ext.pcts, TICsums, by = "species_code", all = FALSE)
ext.sums <- merge(ext.pcts, extr.pct.avg, by = "species_code", all = FALSE)
head(ext.sums)


# using method 20 extraction weight
plot(ext.sums$TICsum ~ ext.sums$percent_extracted, xlim = c(0,80))
lm.ext.sums <- lm(ext.sums$TICsum ~ ext.sums$percent_extracted)
lm.ext.sums.origin <- lm(ext.sums$TICsum ~ ext.sums$percent_extracted - 1)
abline(lm.ext.sums)
abline(lm.ext.sums.origin)
summary(lm.ext.sums)
summary(lm.ext.sums.origin)

# method 20 extraction weight with curved fit
nls(TICsum ~ a*(percent_extracted-c)**b, data = ext.sums, start = list(a = 10000, b = 2, c=20))
curve(10000 * (x-20) ** 2, add = TRUE)

# using phenolics + saponins extraction weight
plot(ext.sums$TICsum ~ ext.sums$phensap_pct, xlim = c(0,0.5))
lm.sapphen <- lm(ext.sums$TICsum ~ ext.sums$phensap_pct - 1)
abline(lm.sapphen)
summary(lm.sapphen)

pct.to.tic <- unname(coef(lm.sapphen)[1])
# change tyrosine column (column 358) to reflect appropriate TIC
for(i in 1:nrow(sampsBycompounds)) {
if(row.names(sampsByCompounds)[i] %in% tyr.pct.avg$species_code) {
  sampsByCompounds$'358'[i] <- tyr.pct.avg[tyr.pct.avg$species_code == row.names(sampsByCompounds)[i], "percent_tyrosine"] * 0.01 * pct.to.tic
}
else {
  sampsByCompounds$'358'[i] <- 0
}}

sampsCompsStand <- sampsByCompounds
for(i in 1:nrow(sampsCompsStand)) {
  sampsCompsStand[i, ] <- sampsCompsStand[i, ] / sum(sampsCompsStand[i, ])
}
sampsCompsStand$'358'

# using just phenolics extraction weight to predict phenolics
plot(ext.sums$phensum ~ ext.sums$phen_pct, xlim = c(0,0.4))
lm.phen <- lm(ext.sums$phensum ~ ext.sums$phen_pct)
abline(lm.phen)
summary(lm.phen)

# using just saponin extraction weight to predict saponins
plot(ext.sums$sapsum ~ ext.sums$sap_pct, xlim = c(0,0.15))
lm.sap <- lm(ext.sums$sapsum ~ ext.sums$sap_pct)
abline(lm.sap)
summary(lm.sap)




