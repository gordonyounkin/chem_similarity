library(RMySQL)

chem_similarity <- read.csv("./results/combined_similarity_matrix_only_mincompclass_2017_11_20.csv")
chem_similarity <- chem_similarity[,names(chem_similarity) != "X"]
row.names(chem_similarity) <- names(chem_similarity)

species_codes <- unique(sapply(1:nrow(chem_similarity), function(x) unlist(strsplit(row.names(chem_similarity)[x], split="_"))[1]))

pairwise.spp <- as.data.frame(matrix(0,nrow = length(species_codes), ncol = length(species_codes)))
names(pairwise.spp) = row.names(pairwise.spp) = species_codes
for(i in 1:length(species_codes)) {
  for(j in i:length(species_codes)) {
    curr_data <- chem_similarity[startsWith(row.names(chem_similarity),paste(species_codes[i],"_",sep="")),
                                 startsWith(row.names(chem_similarity),paste(species_codes[j],"_",sep=""))]
    pairwise.spp[i,j] = pairwise.spp[j,i] = mean(as.matrix(curr_data))
}}

dbDisconnect(mydb)
mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
for(i in 1:nrow(pairwise.spp)) {
  species_name <- dbGetQuery(mydb, paste("SELECT Species_name from Species WHERE species_code = '", row.names(pairwise.spp)[i], "'", sep = ""))
  names(pairwise.spp)[i] <- paste(names(pairwise.spp)[i], species_name, sep = "_")
  row.names(pairwise.spp)[i] <- paste(row.names(pairwise.spp)[i], species_name, sep = "_")
}


write.csv(pairwise.spp, "./results/chem_sim_by_species_mincompclass_2017_11_29.csv", row.names=TRUE)

### extract subset of chemical similarity matrix
pairwise.spp <- read.csv("./results/chem_sim_by_species_1minusdiff_2017_11_29.csv")
pairwise.spp <- pairwise.spp[,names(pairwise.spp) != "X"]
row.names(pairwise.spp) <- names(pairwise.spp)

pairwise.la <- pairwise.spp[startsWith(row.names(pairwise.spp), "LA"), startsWith(names(pairwise.spp), "LA")]

write.csv(pairwise.la, "./results/LA_chem_sim_by_species_1minusdiff_2017_11_30.csv")
