library(reshape2)
library(splitstackshape)

#dbDisconnect(mydb)
#mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')

# load .pairsinfo file from networkedges folder from clean network from GNPS ('View raw spectra')
network_clean <- read.delim(file = "./data/all_compounds_all_libraries/networkedges_selfloop/cc7bd1ebef404cc5a6f9a4e5cb024698.pairsinfo",header=FALSE, sep="\t")
# load .tsv file from main folder in same clean network
compound_cluster_match <- read.delim(file = "./data/all_compounds_all_libraries/METABOLOMICS-SNETS-940a973a-view_raw_spectra-main.tsv",header=TRUE, sep="\t")[,c("cluster.index","ScanNumber")]
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
# remove compounds that are known to be contaminants
contaminants <- read.csv("./data/contaminants_2017_11_16.csv")
nodes <- nodes[!nodes %in% contaminants$compound_number]

# create sampsByCompounds table
# load compound tic table (created in python)
compound_tic_table <- read.csv("./data/compound_table_2018_01_11.csv") 
sampsByCompounds <- dcast(compound_tic_table, compound_number  ~ compound_sample, value.var = "TIC")
head(compound_tic_table)

### if you want pairwise species instead of sample
compound_tic_table_species <- data.frame(compound_tic_table,cSplit(compound_tic_table, 'compound_sample', sep="_", type.convert=FALSE))
compound_tic_table_species_1 <- aggregate(compound_tic_table_species$TIC, by=list(Category=compound_tic_table_species$compound_number,compound_tic_table_species$compound_sample_1),FUN=mean)
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
setdiff(nodes, names(sampsByCompounds))

# remove compounds that were removed from analysis by fill compounds protocol
compound_network <- compound_network[compound_network$NODE1 %in% names(sampsByCompounds) & compound_network$NODE2 %in% names(sampsByCompounds), ]

# for normalized by sample attribute file
sampsCompsStand <- sampsByCompounds
for(i in 1:nrow(sampsByCompounds)){
  sampsCompsStand[i,] = sampsByCompounds[i,]/sum(sampsByCompounds[i,])
}
attribute_file <- data.frame(ceiling(t(sampsCompsStand)*100))
attribute_file$compound <- as.numeric(row.names(attribute_file))

# for attribute file with raw TIC values
attribute_file <- data.frame(t(sampsByCompounds))
attribute_file$compound <- as.numeric(row.names(attribute_file))

# load library hits from GNPS. Download 'View Library Hits' and load .tsv file in main folder
libhits <- read.delim("./data/all_compounds_all_libraries/db_matches.tsv", header = TRUE, sep = "\t", skip = 0, stringsAsFactors = FALSE)
names(libhits)
libhits <- libhits[,c("X.Scan.", "Compound_Name", "INCHI", "Smiles", "PI")]
# convert cluster.index to compound number
names(compound_cluster_match) <- c("X.Scan.","compound")
libhits <- merge(libhits, compound_cluster_match, by="X.Scan.",all.x=TRUE)
libhits <- libhits[, names(libhits) != "X.Scan."]
libhits$INCHI <- as.factor(libhits$INCHI)
libhits$Smiles <- as.factor(libhits$Smiles)
libhits$match_type <- character(nrow(libhits))

# for UNPD matches, need to load key file to convert UNPD key to compound name
unpdkey <- read.csv("K:/DDA/in_silico_fragmentation/UNPD/Prep Files/Database_Key_UNPD_plus_INGA.csv", stringsAsFactors=FALSE)
names(unpdkey)[1] <- "Name"
for(i in 1:nrow(libhits)) {
  if(startsWith(libhits$Compound_Name[i], "UNPD")) {
    libhits$Compound_Name[i] <- unpdkey[unpdkey$ID==libhits$Compound_Name[i], "Name"]
    libhits$match_type[i] <- "in silico UNPD"
  }
  else {
    if(libhits$PI[i] == "TAKursar") libhits$match_type[i] <- "in silico"
    else libhits$match_type[i] <- "GNPS databases"
  }
}

# load mz & rt data for compounds
compound_feature <- read.csv("./data/compound_feature_table.csv")
compound_max_feature <- aggregate(compound_feature$TIC, by=list(Category=compound_feature$compound_number),FUN=max)
names(compound_max_feature) <- c("compound","TIC_max")
compound_max_feature_2 <- merge(compound_max_feature, compound_feature, by.x=c("compound","TIC_max"), by.y=c("compound_number","TIC"))

# merge library hits and mz/rt to attribute file
attribute_file <- merge(attribute_file, libhits, by = "compound", all.x = TRUE, all.y = FALSE)
attribute_file <- merge(attribute_file, compound_max_feature_2[,c("compound", "mz", "rt")], by = "compound", all.x = TRUE, all.y = FALSE)
attribute_file$Compound_Name <- as.character(attribute_file$Compound_Name)
attribute_file$Smiles <- as.character(attribute_file$Smiles)

# for compounds that don't have any matches, match to known compound database based on mass only
known_comps <- read.csv("K:/Lab_Map/METHODS FOR CHEMISTRY_pdfs_protocols/6_our database for MS/6b_Known_compounds/Known_Compounds_020.csv", stringsAsFactors=FALSE)
known_comps <- known_comps[!known_comps$CS_MM %in% c("","#VALUE!"), ]
known_comps$CS_MM <- as.numeric(as.character(known_comps$CS_MM))

for(i in 1:nrow(attribute_file)) {
 if(is.na(attribute_file$Compound_Name[i])) {
   current_mass <- attribute_file$mz[i] + 1.007
   if(sum(abs(current_mass - known_comps$CS_MM) < 0.01) > 0) {
     attribute_file$Compound_Name[i] <- paste(known_comps[which(abs(current_mass - known_comps$CS_MM) == min(abs(current_mass - known_comps$CS_MM), na.rm = TRUE)),"putative.name"], collapse = ",")
     attribute_file$Smiles[i] <- paste(known_comps[which(abs(current_mass - known_comps$CS_MM) == min(abs(current_mass - known_comps$CS_MM), na.rm = TRUE)),"SMILES"], collapse = ",")
     attribute_file$match_type[i] <- "mass_only" 
   }}}

# write files to use in cytoscape
write.table(attribute_file,"./data/cytoscape_all_inga_C18_neg/attribute_file_3.out", sep=",",row.names=FALSE, col.names = TRUE)
write.table(compound_network, "./data/cytoscape_all_inga_C18_neg/compound_network.tsv", sep = "\t", row.names=FALSE)


# list of all compounds id'ed by GNPS:
ided_comps <- attribute_file[!is.na(attribute_file$match_type), c("compound", "Compound_Name","Smiles","match_type","mz","rt")]
write.csv(ided_comps, "./results/compounds_ided_by_gnps.csv", row.names = FALSE)
names(attribute_file)
