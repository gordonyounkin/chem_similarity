library(RMySQL)

chem_similarity <- read.csv("./results/combined_similarity_matrix_3RT_2017_12_07.csv")
chem_similarity <- chem_similarity[,names(chem_similarity) != "X"]
row.names(chem_similarity) <- names(chem_similarity)

# remove bad samples
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "LA20b_1745"),!startsWith(names(chem_similarity), "LA20b_1745")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "LA9_696"),!startsWith(names(chem_similarity), "LA9_696")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "N1_1250"),!startsWith(names(chem_similarity), "N1_1250")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "T82_1291"),!startsWith(names(chem_similarity), "T82_1291")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "LA7_691_2"),!startsWith(names(chem_similarity), "LA7_691_2")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "M18_612_2"),!startsWith(names(chem_similarity), "M18_612_2")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "T86_1257"),!startsWith(names(chem_similarity), "T86_1257")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "N4_1473"),!startsWith(names(chem_similarity), "N4_1473")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "N4_1472"),!startsWith(names(chem_similarity), "N4_1472")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "N31_1471"),!startsWith(names(chem_similarity), "N31_1471")]
chem_similarity <- chem_similarity[!startsWith(names(chem_similarity), "LA23_1479"),!startsWith(names(chem_similarity), "LA23_1479")]


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


write.csv(pairwise.spp, "./results/species_similarity_3RT_2017_12_11.csv", row.names=TRUE)

### extract subset of chemical similarity matrix
pairwise.spp <- read.csv("./results/chem_sim_by_species_1minusdiff_2017_11_29.csv")
pairwise.spp <- pairwise.spp[,names(pairwise.spp) != "X"]
row.names(pairwise.spp) <- names(pairwise.spp)

pairwise.la <- pairwise.spp[startsWith(row.names(pairwise.spp), "LA"), startsWith(names(pairwise.spp), "LA")]

write.csv(pairwise.la, "./results/LA_species_similarity_3RT_2017_12_11.csv")
