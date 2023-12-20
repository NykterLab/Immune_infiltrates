


d=read.table("All_samples.txt")
viable_stained <- d$stainInfo=="stained" & grepl("Live",d$fileID)

channels <- c("FSC_A","SSC_A",
              "CD66b","CD4","CD19",
              "CD14","CD8","CD45","CD3")

cluster_PhenoGraph <- Rphenograph::Rphenograph(data = d[viable_stained,channels], k=60)
cluster_PhenoGraph_vector <- as.vector(membership(cluster_PhenoGraph[[2]]))

out <- d[viable_stained,]
d$clusterID <- cluster_PhenoGraph_vector

write.table(d$clusterID,"Viable_stained_data.txt")

#########
#Filter out CD45 negative and recluster

d <- read.table("viable_stained_data.txt",header=T)

cluster_median <- aggregate(.~clusters,data = d[,c(1:11,16)],median)
#Threshold based on manula inspection
pick <- cluster_median$clusters[which(cluster_median$CD45 > 1)]

data <- d[d$clusters %in% pick,]

channels <- c("FSC_A","SSC_A",
              "CD66b","CD4","CD19",
              "CD14","CD8","CD45","CD3")

cluster_PhenoGraph <- Rphenograph::Rphenograph(data = data[,channels], k=60)
cluster_PhenoGraph <- as.vector(membership(cluster_PhenoGraph[[2]]))

data_tsne <-  Rtsne(X = data[,channels],initial_dims= length(channels),check_duplicates = FALSE)
data_tsne <- data_tsne$Y

data$old_clusters <- data$clusters
data$clusters <- cluster_PhenoGraph
data$tsne1 <- data_tsne[,1]
data$tsne2 <- data_tsne[,2]

write.table(data,"CD45_data.txt")