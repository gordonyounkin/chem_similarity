library(mzR)
library(RMySQL)
library(splitstackshape)
library(xcms)
library(reshape2)
source("./code/create_pairwiseComps_sampsByComps.R")

comp_table <- read.csv("./data/polar_compound_table_2017_12_15.csv")
head(comp_table)
comp_table_2 <- reshape2::dcast(comp_table, compound_number ~ TIC)
comp_table <- make_sampsByCompounds("./data/polar_compound_table_2017_12_15.csv")


missTIC_2 <- dcast(missTIC, ct_compound_number ~ ct_compound_sample)
row.names(missTIC_2) <- missTIC_2$ct_compound_number
missTIC_2.1 <- missTIC_2[,2:ncol(missTIC_2)]
missTIC_3 <- missTIC_2.1[sapply(1:nrow(missTIC_2.1), function(x) max(missTIC_2.1[x,]))>=5000,]
if(nrow(missTIC_3) < 1) next
missTIC_4 <- data.frame(missTIC_3>=5000)
missTIC_4.1 <- missTIC_4

samples <- data.frame("name"=unique(missTIC$ct_compound_sample), "ncomps"=colSums(missTIC_4), "removable"=TRUE)

while( sum(!samples$removable) <= (nrow(samples)-1) & nrow(samples)>1 ) {
  for(i in 1:nrow(samples)) {
    samples$removable[i] <- sum(rowSums(data.frame(missTIC_4.1[,-i]))>0) == sum(rowSums(missTIC_4.1)>0)
  }
  if( sum(!samples$removable) <= (nrow(samples)-1) ) {
    missTIC_4.1 <- missTIC_4.1[,names(missTIC_4.1)!=samples[samples$removable==TRUE,"name"][which.min(samples[samples$removable==TRUE,"ncomps"])]]
    samples <- samples[rownames(samples)!=samples[samples$removable==TRUE,"name"][which.min(samples[samples$removable==TRUE,"ncomps"])],]
  }
}
missTIC_5 <- data.frame(missTIC_3[,names(missTIC_3)%in%samples$name])
names(missTIC_5) <- names(missTIC_3[names(missTIC_3)%in%samples$name])
row.names(missTIC_5) <- row.names(missTIC_3)
missTIC_5$maxSample <- sapply(1:nrow(missTIC_5), function(x) names(missTIC_5)[which.max(missTIC_5[x,])])
missTIC_6 <- data.frame("Compound_Number"=row.names(missTIC_5),"sample"=missTIC_5$maxSample, "species"=species[n])
all_samp_targets <- rbind(all_samp_targets, missTIC_6)
}

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



