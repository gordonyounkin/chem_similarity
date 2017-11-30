library(RMySQL)

chem_similarity <- read.csv("./results/combined_similarity_matrix_only_compclassavg_2017_11_20.csv")
chem_similarity <- chem_similarity[,names(chem_similarity) != "X"]
row.names(chem_similarity) <- names(chem_similarity)


all_samps.df <-  data.frame(t(combn(names(chem_similarity),2)),chemdist=t(chem_similarity)[lower.tri(chem_similarity)])
all_samps.df$X1 <- as.character(all_samps.df$X1)
all_samps.df$X2 <- as.character(all_samps.df$X2)
all_samps.df$spec_code_1 <- sapply(1:nrow(all_samps.df), function(x) unlist(strsplit(all_samps.df$X1[x], split="_"))[1])
all_samps.df$spec_code_2 <- sapply(1:nrow(all_samps.df), function(x) unlist(strsplit(all_samps.df$X2[x], split="_"))[1])
all_samps.df$spec_codes <- sapply(1:nrow(all_samps.df), function(x) paste(all_samps.df$spec_code_1[x], all_samps.df$spec_code_2[x], sep="_"))
head(all_samps.df)

spec_codes <- unique(c(as.character(all_samps.df$spec_code_1), as.character(all_samps.df$spec_code_2)))
spec_code_species <- data.frame("species_code" = spec_codes, "species_name" = character(length(spec_codes)), stringsAsFactors = FALSE)
dbDisconnect(mydb)
mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')
for(i in 1:nrow(spec_code_species)) {
  species_name <- dbGetQuery(mydb, paste("SELECT Species_name from Species WHERE species_code = '", spec_code_species$species_code[i], "'", sep = ""))
  spec_code_species$species_name[i] <- species_name
}
table(unlist(spec_code_species$species_name))
spec_code_species[spec_code_species$species_code == "M63", "species_name"] <- "umbellifera"

all_samps.df$species_1 <- unlist(sapply(1:nrow(all_samps.df), function(x) spec_code_species[spec_code_species$species_code == all_samps.df$spec_code_1[x], "species_name"]))
all_samps.df$species_2 <- unlist(sapply(1:nrow(all_samps.df), function(x) spec_code_species[spec_code_species$species_code == all_samps.df$spec_code_2[x], "species_name"]))

all_samps.df$sample_pair_type <- character(nrow(all_samps.df))
for(i in 1:nrow(all_samps.df)) {
  if(all_samps.df$spec_code_1[i] == all_samps.df$spec_code_2[i]) all_samps.df$sample_pair_type[i] <- "same_spec_same_site"
  else {
    if(is.na(all_samps.df$species_1[i]) | is.na(all_samps.df$species_2[i])) {
      all_samps.df$sample_pair_type[i] <- "all_others"
    }
    else {
      if(all_samps.df$species_1[i] == all_samps.df$species_2[i]) all_samps.df$sample_pair_type[i] <- "same_spec_diff_site"
      else all_samps.df$sample_pair_type[i] <- "all_others"
    }
  }
}

dev.new()
boxplot(all_samps.df[all_samps.df$sample_pair_type == "all_others", "chemdist"],
  all_samps.df[all_samps.df$sample_pair_type == "same_spec_diff_site", "chemdist"],
  all_samps.df[all_samps.df$sample_pair_type == "same_spec_same_site", "chemdist"],
  names = c("all_others", "same species diff site", "same species/site"),
  main = "similarity calcs using ln(TIC)",
  ylab = "similarity score")
dev.copy2pdf(file = "./results/sim_scores_by_sample_pair_type_LOGTIC_boxplot.pdf")
dev.off()


all_samps.df_1mindiff 
all_samps.df_mincompclass 
all_samps.df_compclassavg 

all_samps.df_1mindiff[all_samps.df_1mindiff$sample_pair_type == "same_spec_same_site" & all_samps.df_1mindiff$chemdis < 0.5, ]

write.csv(all_samps.df[,names(all_samps.df) != "sample_pair_type"], "./results/all_samples_chem_scores_long_format.csv")
