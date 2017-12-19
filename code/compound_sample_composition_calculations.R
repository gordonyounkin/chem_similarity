source("./code/create_pairwiseComps_sampsByComps.R")
source("./code/chem_similarity_function.R")

filled_comps <- read.csv("./data/filled_compound_table_2017_12_05.csv", stringsAsFactors=FALSE)
filled_comps <- read.csv("./data/filled_compound_table_major_features_2017_11_30.csv", stringsAsFactors=FALSE)
head(filled_comps)
table(filled_comps$compound_sample)
range(table(filled_comps$compound_sample))
table(filled_comps$compound_sample)[table(filled_comps$compound_sample)==508]
mean(table(filled_comps$compound_sample))
length(unique(filled_comps$compound_number))

# calculate relative abundance of compounds within each sample
filled_comps$relTIC <- numeric(nrow(filled_comps))
for(i in 1:length(unique(filled_comps$compound_sample))) {
  sample <- unique(filled_comps$compound_sample)[i]
  temp <- filled_comps[filled_comps$compound_sample == sample, ]
  maxtic <- max(temp$TIC)
  for(j in 1:nrow(temp)) {
    filled_comps[filled_comps$compound_sample == sample & filled_comps$compound_number == temp$compound_number[j], "relTIC"] <- temp$TIC[j] / maxtic
  }
}
filled_comps_trim <- filled_comps[filled_comps$relTIC > 0.001,]
table(filled_comps_trim$compound_sample)
range(table(filled_comps_trim$compound_sample))
table(filled_comps_trim$compound_sample)[table(filled_comps_trim$compound_sample)==32]
mean(table(filled_comps_trim$compound_sample))
length(unique(filled_comps_trim$compound_number))

# do species tend to share a lot of compounds?
filled_comps_species <- filled_comps[startsWith(filled_comps$compound_sample, "LA4_"),]


hist(filled_comps_species$relTIC)

table(filled_comps_species$compound_sample)
table(filled_comps_trim$compound_sample)

length(table(filled_comps_species$compound_number)[table(filled_comps_species$compound_number) %in% c(1)])
sum(filled_comps_species[filled_comps_species$compound_number %in% as.numeric(names(table(filled_comps_species$compound_number)[table(filled_comps_species$compound_number) %in% c(4)])), "TIC"])
hist(log(filled_comps_species[filled_comps_species$compound_number %in% as.numeric(names(table(filled_comps_species$compound_number)[table(filled_comps_species$compound_number) %in% c(1)])) & filled_comps_species$compound_sample == "LA4_1440", "TIC"]))


filled_comps[filled_comps$compound_sample == "N34_1336",]

hist(table(filled_comps$compound_number), breaks = 50, main = "With filled peaks", xlab = "Number of samples compound is found in")

# which compounds are super common
table(filled_comps$compound_number)[table(filled_comps$compound_number) > 600]


old_comp_table <- read.csv("./data/compound_tic_2017_08_02.csv")
head(old_comp_table)
table(old_comp_table$compound_sample)
range(table(old_comp_table$compound_sample))
mean(table(old_comp_table$compound_sample))
table(old_comp_table$compound_sample)[table(old_comp_table$compound_sample)==224]

hist(table(old_comp_table$compound_number), breaks = 50, main = "XCMS output", xlab = "Number of samples compound is found in")


filled_sampsbycomps <- make_sampsByCompounds("./data/filled_compound_table_2017_12_05.csv",by_species = FALSE)
filled_sampsbycomps_log <- log(filled_sampsbycomps)
filled_sampsbycomps_log[filled_sampsbycomps_log<=0] <- 0
filled_sampscompsstand <- standardizeByRow(filled_sampsbycomps)
filled_sampsbycomps_log <- standardizeByRow(filled_sampsbycomps_log)
filled_sampsbycomps_sqrt <- sqrt(filled_sampsbycomps)
filled_sampcompsstand_sqrt <- standardizeByRow(filled_sampsbycomps_sqrt)
filled_sampsbycomps_4rt <- filled_sampsbycomps**(1/4)
filled_sampscompsstand_4rt <- standardizeByRow(filled_sampsbycomps_4rt)

curr_sample <- unlist(filled_sampscompsstand[2,,drop=TRUE])
curr_sample <- sort(curr_sample, decreasing = TRUE)
plot(as.numeric(as.vector(curr_sample[curr_sample>0])) ~ as.numeric(1:sum(curr_sample>0)))

# without any transformation
row.names(filled_sampsbycomps)[504]
curr_sample <- unlist(filled_sampscompsstand[504,,drop=TRUE])
curr_sample <- sort(curr_sample, decreasing = TRUE)
plot(as.numeric(as.vector(curr_sample[curr_sample>0])) ~ as.numeric(1:sum(curr_sample>0)), main = "raw TIC")
range(curr_sample[curr_sample>0])
hist(curr_sample[curr_sample>0])

# log of TIC
curr_sample <- unlist(filled_sampsbycomps_log[504,,drop=TRUE])
curr_sample <- sort(curr_sample, decreasing = TRUE)
plot(as.numeric(as.vector(curr_sample[curr_sample>0])) ~ as.numeric(1:sum(curr_sample>0)), main = "log(TIC)")
range(curr_sample[curr_sample>0])
hist(curr_sample[curr_sample>0])

# square root of TIC
curr_sample <- unlist(filled_sampcompsstand_sqrt[504,,drop=TRUE])
curr_sample <- sort(curr_sample, decreasing = TRUE)
plot(as.numeric(as.vector(curr_sample[curr_sample>0])) ~ as.numeric(1:sum(curr_sample>0)), main = "TIC^(1/2)")
range(curr_sample[curr_sample>0])
hist(curr_sample[curr_sample>0])

# 4th root of TIC
curr_sample <- unlist(filled_sampscompsstand_4rt[504,,drop=TRUE])
curr_sample <- sort(curr_sample, decreasing = TRUE)
plot(as.numeric(as.vector(curr_sample[curr_sample>0])) ~ as.numeric(1:sum(curr_sample>0)), main = "TIC^(1/4)")
range(curr_sample[curr_sample>0])
hist(curr_sample[curr_sample>0])

old_sampsbycomps <- make_sampsByCompounds("./data/compound_tic_2017_08_02.csv", by_species = FALSE)




filled_sampsbycomps[1,]
