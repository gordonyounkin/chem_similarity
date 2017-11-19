# load chem similarity file
pairwise.spp <- read.csv("./results/2017_10_31_pairwise.samps.avgcompclass.csv")
row.names(pairwise.spp) <- pairwise.spp$X
pairwise.spp <- pairwise.spp[,names(pairwise.spp)!="X"]

# create chem similarity tree from similarity matrix (similarity matrix should be named 'pairwise.spp')
library(vegan)
library(pvclust)

species_dist <- vegdist(pairwise.spp,"euclidean")      		# Distance matrix
species_single <- hclust(species_dist, method="single")					# Single linkage clustering
species_ward <- hclust(species_dist, method="ward.D")					# Ward clustering
species_complete <- hclust(species_dist, method="complete")				# Complete linkage clustering
species_centroid <- hclust(species_dist, method="centroid")				# Centroid clustering
species_median <- hclust(species_dist, method="median")				        # Median clustering
#Comparison between the distance matrix and binary matrices representing partitions 
coph1 <- cophenetic(species_single)							# Compute Patristic distances		
coph2 <- cophenetic(species_ward)
coph3 <- cophenetic(species_complete)
coph4 <- cophenetic(species_centroid)
coph5 <- cophenetic(species_median)
a <-cor(coph1, species_dist)								#Cophenetic correlations
b <-cor(coph2, species_dist)
c <-cor(coph3, species_dist)
d <-cor(coph4, species_dist)
e <-cor(coph5, species_dist)
method <- c("single","ward.D","complete", "centroid", "median")
mhc<- method[which.max(c(a,b,c,d,e))]
result_samples <- pvclust(pairwise.spp, method.hclust=mhc, method.dist="correlation", use.cor="pairwise.complete.obs", nboot=1000,parallel=T)

# save as .pdf--make sure to give it a name
dev.new()
plot(result_samples, cex=1.66, cex.pv=1, lwd=1, float = 0.003)
dev.copy2pdf(file = "./data/BCI_test/BCI_chem_dendrogram_filled_LOG.pdf", width = 50, height = 20)
dev.off()

# save as .tre
