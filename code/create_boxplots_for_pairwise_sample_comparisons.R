library(RMySQL)

chem_similarity <- read.csv("./results/combined_similarity_matrix_3RT_2017_12_07.csv")
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
spec_code_species[spec_code_species$species_code == "T05", "species_name"] <- "sapindoides"
spec_code_species[spec_code_species$species_code %in% c("N65","N67", "M63"), "species_name"] <- "umbellifera"

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

# throw out outliers??
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "LA20b_1745"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "LA20b_1745"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "LA9_696"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "LA9_696"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "N1_1250"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "N1_1250"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "T82_1291"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "T82_1291"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "LA7_691_2"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "LA7_691_2"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "M18_612_2"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "M18_612_2"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "T86_1257"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "T86_1257"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "N4_1473"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "N4_1473"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "N4_1472"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "N4_1472"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "N31_1471"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "N31_1471"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X2, "LA23_1479"),]
all_samps.df <- all_samps.df[!startsWith(all_samps.df$X1, "LA23_1479"),]



# what do these similarity scores look like for species in sawfly analysis?
sawfly_species <- c("T67a","IngC", "N42","N29","N37","N38","LA9","T63","T19","N13","N35","N6","T59","LA43","T29","N31","T69","T82","N41","T65","N26","T60","N7","N66","IngG","IngU","T57","N5","LA3","T74","LA32","T23","N21","T76","T35","T33","LA33","N23","N10","LA35","N8","IngF","T20","N25","LA60")
sawfly_samps <- all_samps.df[all_samps.df$spec_code_1 %in% sawfly_species & all_samps.df$spec_code_2 %in% sawfly_species, ]

all_samps.df_3RT <- all_samps.df

all_samps.df <- all_samps.df_SQRTSQRT
all_samps.df <- all_samps.df_4RT
all_samps.df <- all_samps.df_3RT
all_samps.df <- all_samps.df_SQRT
all_samps.df <- all_samps.df_LOG
#dev.new()
boxplot(all_samps.df[all_samps.df$sample_pair_type == "all_others", "chemdist"],
  all_samps.df[all_samps.df$sample_pair_type == "same_spec_diff_site", "chemdist"],
  all_samps.df[all_samps.df$sample_pair_type == "same_spec_same_site", "chemdist"],
  names = c("all_others", "same species diff site", "same species/site"),
  main = "similarity calcs using sqrt(TIC) cw sqrt",
  ylab = "similarity score",
  ylim = c(0,1))
#dev.copy2pdf(file = "./results/sim_scores_by_sample_pair_type_LOGTIC_boxplot.pdf")
#dev.off()

hist(all_samps.df[all_samps.df$sample_pair_type == "all_others", "chemdist"], xlim=c(0,1), main="Different species similarity using sqrt(TIC)")
hist(all_samps.df[all_samps.df$sample_pair_type == "same_spec_diff_site", "chemdist"])
hist(all_samps.df[all_samps.df$sample_pair_type == "same_spec_same_site", "chemdist"], xlim = c(0,1), ylim=c(0,500),main = "Same species/site similarity calcs using 4rt(TIC) cw sqrt")

plot(all_samps.df_SQRT$chemdist ~ all_samps.df_LOG$chemdist)
plot(all_samps.df_SQRT[all_samps.df_SQRT$sample_pair_type == "same_spec_same_site", "chemdist"] ~ all_samps.df_SQRTSQRT[all_samps.df_SQRTSQRT$sample_pair_type == "same_spec_same_site", "chemdist"])

curve(x*1, add=TRUE, lwd=2, col="red")

mean(all_samps.df[all_samps.df$sample_pair_type == "same_spec_same_site", "chemdist"])

all_samps.df[all_samps.df$sample_pair_type == "same_spec_same_site" & all_samps.df$chemdist < 0.5, ]
all_samps.df[all_samps.df$sample_pair_type == "all_others" & all_samps.df$chemdist > 0.8, ]
all_samps.df[all_samps.df$sample_pair_type == "same_spec_diff_site" & all_samps.df$chemdist < 0.2, ]

write.csv(all_samps.df[,names(all_samps.df) != "sample_pair_type"], "./results/all_samples_chem_scores_long_format.csv")

# How much variation is there among chem of samples from same accession?
head(all_samps.df)
all_samps_same_species <- all_samps.df[all_samps.df$sample_pair_type == "same_spec_same_site",]
species <- unique(all_samps_same_species$spec_code_1)
species_variation <- data.frame("species" = species, "avg_score" = numeric(length(species)), "score_stdev" = numeric(length(species)))
for(i in 1:length(species)) {
  species_scores <- all_samps_same_species[all_samps_same_species$spec_code_1 == species[i], ]
  species_variation$avg_score[i] <- mean(species_scores$chemdist)
  species_variation$score_stdev[i] <- sd(species_scores$chemdist)
}

hist(species_variation$avg_score)
hist(species_variation$score_stdev)

species_variation[species_variation$avg_score < 0.6,]
species_variation[species_variation$score_stdev > 0.15,]

all_samps_same_species[all_samps_same_species$spec_code_1 == "IngU", ]
hist(all_samps_same_species[all_samps_same_species$spec_code_1 == "IngU", "chemdist"])
comp.class.pcts[startsWith(comp.class.pcts$sample, "LA3_"),]
chem_similarity_phen[startsWith(names(chem_similarity_phen), "LA10_"), startsWith(names(chem_similarity_phen), "T58")]
chem_similarity_sap[startsWith(names(chem_similarity_sap), "LA10_"), startsWith(names(chem_similarity_sap), "T58")]
chem_similarity[startsWith(names(chem_similarity_phen), "T05_"), startsWith(names(chem_similarity_phen), "T58")]

plot(species_variation$avg_score ~ species_variation$score_stdev)


