library(mzR)
library(RMySQL)
library(splitstackshape)
library(xcms)
library(reshape2)
source("./code/create_pairwiseComps_sampsByComps.R")

mydb = dbConnect(MySQL(), user='u6009010', password='3UaUhf7a', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')

comp_table <- make_sampsByCompounds("./data/polar_compound_table_2017_12_15.csv")
comp_table_2 <- as.matrix(comp_table)
comp_table_2[comp_table_2 < 1000] <- 0
comp_table_2[comp_table_2 != 0] <- 1
comp_table_3 <- as.data.frame(t(comp_table_2))
comp_table_3$keep <- character(nrow(comp_table_3))

# only include compounds that are present in at least 3/5 of the samples of any one species
for(i in 1:ncol(comp_table_3)) {
  EN <- gsub("DN","", unlist(strsplit(names(comp_table_3[i]), split="_"))[1])
  species_code <- dbGetQuery(mydb, paste("SELECT species_code FROM `Extraction` WHERE Extraction_Number = ", EN, sep=""))
  names(comp_table_3)[i] <- paste("DN", EN, "_", species_code, sep="")
}

species_compounds <- data.frame("species" = sapply(1:ncol(comp_table_3), function(x) unlist(strsplit(names(comp_table_3)[x], split="_"))[2]), "present" = numeric(ncol(comp_table_3)))

for(i in 1:nrow(comp_table_3)) {
  species_compounds$present <- comp_table_3[i, names(comp_table_3) != "keep", drop=TRUE]
  abundances  <- table(species_compounds[species_compounds$present == 1, "species"])
  proportions <- abundances / table(species_compounds$species)
  if(max(proportions) >= 0.6 & max(abundances > 2)) comp_table_3$keep[i] <- TRUE
  else comp_table_3$keep[i] <- FALSE
  }

comp_table_4 <- comp_table_3[comp_table_3$keep == TRUE, names(comp_table_3) != "keep"]

rowSums(comp_table_4)

samples <- data.frame("name"=names(comp_table_4), "ncomps"=colSums(comp_table_4>0), "removable"=TRUE, stringsAsFactors = FALSE)

while( sum(!samples$removable) <= (nrow(samples)-1) & nrow(samples)>1 ) {
  print(nrow(samples))
  for(i in 1:nrow(samples)) {
    samples$removable[i] <- sum(rowSums(data.frame(comp_table_4[,-i]))>0) == sum(rowSums(comp_table_4)>0)
  }
  if( sum(!samples$removable) <= (nrow(samples)-1) ) {
    comp_table_4 <- comp_table_4[,names(comp_table_4)!=samples[samples$removable==TRUE,"name"][which.min(samples[samples$removable==TRUE,"ncomps"])]]
    samples <- samples[rownames(samples)!=samples[samples$removable==TRUE,"name"][which.min(samples[samples$removable==TRUE,"ncomps"])],]
  }
}


head(comp_table_4)
rowSums(comp_table_4)





comp_table_5 <- data.frame(comp_table_4[,names(comp_table_4)%in%samples$name])
names(comp_table_5) <- names(comp_table_4[names(comp_table_4)%in%samples$name])
row.names(comp_table_5) <- row.names(comp_table_4)
comp_table_5$maxSample <- sapply(1:nrow(comp_table_5), function(x) names(comp_table_5)[which.max(comp_table_5[x,])])
comp_table_6 <- data.frame("Compound_Number"=row.names(comp_table_5),"sample"=comp_table_5$maxSample, "species"=species[n])
all_samp_targets <- rbind(all_samp_targets, comp_table_6)


msms_sample_targets <- merge(dda_misses, all_samp_targets, by=c("Compound_Number","species"), all.x=FALSE, all.y=TRUE)
msms_sample_targets$mz <- as.numeric(msms_sample_targets$mz)
head(msms_sample_targets)
unique(msms_sample_targets$species)
table(msms_sample_targets$species)
unique(msms_sample_targets$sample)
sum(table(msms_sample_targets$sample)>32)

# separate into saponins and non-saponins
sap_targets <- msms_sample_targets[msms_sample_targets$mz >= 550 & msms_sample_targets$rt >= 960, ]
nonsap_targets <- msms_sample_targets[msms_sample_targets$mz <= 550 | msms_sample_targets$rt <= 960, ]

# create target lists for MSMS - SAPONINS
msms_samples <- unique(sap_targets$sample)
for(i in 1:length(msms_samples)) {
  targets <- sap_targets[sap_targets$sample==msms_samples[i],c("mz","rt")]
for(k in 1:ceiling(nrow(targets)/32)) {
  if(32*k<=nrow(targets)) {
    assign(paste("targets_",k,sep=""),targets[((k*32)-31):(32*k),])
  }
  else {
    assign(paste("targets_",k,sep=""),targets[((k*32)-31):nrow(targets),])  
  }
  #### write csv of target peaks into correct folder
  write.csv(get(paste("targets_",k,sep="")),paste("K://DDA//FG//msms_targets//saponins//",msms_samples[i],"_MSMS_sap_targets_",k,".csv",sep=""),row.names=FALSE)
}
}

# create target lists for MSMS - NONSAPONINS
msms_samples <- unique(nonsap_targets$sample)
for(i in 1:length(msms_samples)) {
  targets <- nonsap_targets[nonsap_targets$sample==msms_samples[i],c("mz","rt")]
  for(k in 1:ceiling(nrow(targets)/32)) {
    if(32*k<=nrow(targets)) {
      assign(paste("targets_",k,sep=""),targets[((k*32)-31):(32*k),])
    }
    else {
      assign(paste("targets_",k,sep=""),targets[((k*32)-31):nrow(targets),])  
    }
    #### write csv of target peaks into correct folder
    write.csv(get(paste("targets_",k,sep="")),paste("K://DDA//FG//msms_targets//non_saponins//",msms_samples[i],"_MSMS_targets_",k,".csv",sep=""),row.names=FALSE)
  }
}


# where are these extractions located
extr_locs <- dbGetQuery(mydb, paste("SELECT species_code, Extraction_Number, Box_Number FROM `Extraction`
                        WHERE Extraction_Number IN (",paste(sapply(1:length(unique(nonsap_targets$sample)), function(x) strsplit(unique(as.character(nonsap_targets$sample))[x],split="_")[[1]][2]),collapse=", "),")",sep=""))
head(extr_locs)
table(extr_locs$species_code)
table(dda_misses$species)
write.csv(extr_locs, "K:/DDA/FG/msms_targets/nonsap_FG_extraction_locs_MSMS_2017_06_27.csv", row.names=FALSE)

# for saponins
extr_locs <- dbGetQuery(mydb, paste("SELECT species_code, Extraction_Number, Box_Number FROM `Extraction`
                        WHERE Extraction_Number IN (",paste(sapply(1:length(unique(sap_targets$sample)), function(x) strsplit(unique(as.character(sap_targets$sample))[x],split="_")[[1]][2]),collapse=", "),")",sep=""))
head(extr_locs)
table(extr_locs$species_code)
table(dda_misses$species)
write.csv(extr_locs, "K:/DDA/FG/msms_targets/sap_extraction_locs_MSMS_2017_06_29.csv", row.names=FALSE)



