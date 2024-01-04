source("import_facs.r")

dir <- "../../fcs/"
files <- list.files(dir,recursive = TRUE,full.names = TRUE,pattern="*fcs")

# files <- c(
#   "../../fcs//1.- 255/CD45+ S R14JH255_004_Dead gate.fcs",
#   "../../fcs//1.- 255/CD45+ S R14JH255_004_Live cells.fcs",
#   "../../fcs//1.- 255/CD45+ U R14JH255_003_Dead gate.fcs",
#   "../../fcs//1.- 255/CD45+ U R14JH255_003_Live cells.fcs",
#   "../../fcs//2.- 258/CD45+ S R14JH258_004_Dead gate.fcs", 
#   "../../fcs//2.- 258/CD45+ S R14JH258_004_Live cells.fcs",
#   "../../fcs//2.- 258/CD45+ U R14JH258_003_Dead gate.fcs", 
#   "../../fcs//2.- 258/CD45+ U R14JH258_003_Live cells.fcs",
#   "../../fcs//3.- 264/CD45+ S R14JH264_004_Dead gate.fcs", 
#   "../../fcs//3.- 264/CD45+ S R14JH264_004_Live cells.fcs",
#   "../../fcs//3.- 264/CD45+ U R14JH264_003_Dead gate.fcs", 
#   "../../fcs//3.- 264/CD45+ U R14JH264_003_Live cells.fcs",
#   "../../fcs//4.- 009/CD45+ S R14JH009_L1_Dead gate.fcs",  
#   "../../fcs//4.- 009/CD45+ S R14JH009_L1_Live cells.fcs", 
#   "../../fcs//4.- 009/CD45+ S R14JH009_L2_Dead gate.fcs",  
#   "../../fcs//4.- 009/CD45+ S R14JH009_L2_Live cells.fcs", 
#   "../../fcs//4.- 009/CD45+ U R14JH009_L1_Dead gate.fcs",  
#   "../../fcs//4.- 009/CD45+ U R14JH009_L1_Live cells.fcs", 
#   "../../fcs//4.- 009/CD45+ U R14JH009_L2_Dead gate.fcs",  
#   "../../fcs//4.- 009/CD45+ U R14JH009_L2_Live cells.fcs", 
#   "../../fcs//5.- 261/CD45+ S R14JH261_L1_Dead gate.fcs",  
#   "../../fcs//5.- 261/CD45+ S R14JH261_L1_Live cells.fcs", 
#   "../../fcs//5.- 261/CD45+ S R14JH261_L2_Dead gate.fcs",  
#   "../../fcs//5.- 261/CD45+ S R14JH261_L2_Live cells.fcs", 
#   "../../fcs//5.- 261/CD45+ U R14JH261_L1_Dead gate.fcs",  
#   "../../fcs//5.- 261/CD45+ U R14JH261_L1_Live cells.fcs", 
#   "../../fcs//5.- 261/CD45+ U R14JH261_L2_Dead gate.fcs",  
#   "../../fcs//5.- 261/CD45+ U R14JH261_L2_Live cells.fcs" 
# )

strsplits <- unlist(lapply(strsplit(files,"/"),function(x) x[length(x)]))
strsplits <- strsplit(strsplits," ")
stained <- unlist(lapply(strsplits, function(x) x[[2]]))
fileIDs <- unlist(lapply(strsplits, function(x) x[[3]]))
fileIDs <- paste0(fileIDs,"_",stained)
stained <- stained == "S"
stainStatus <- rep("stained",length(stained))
stainStatus[!stained] <- "unstained"

sampleIDs <- c("R14JH255","R14JH258","R14JH264","R14JH261_L1","R14JH261_L2","R14JH251","R14JH236","R14JH-227_L1","R14JH-227_L2","R14JH009_L1","R14JH009_L2")
#which(grepl(sampleIDs[1],fileIDs))
samples <- rep(0,length(fileIDs))
for(i in sampleIDs){
  pick <- grepl(i,fileIDs)
  samples[pick] <- i
}

flowGroup <- rep(1,length(files))
#merges the viable, dead-gated and unstained data of one particular sample to ensure same transformation within a sample
dataList <- mergeAutoLgcl(files,pick = "all",
                           fileIDs = fileIDs,
                           stainStatus = stainStatus,
                           flowGroup = flowGroup,
                           sampleIDs=samples)


data <- data.table::rbindlist(dataList,fill=TRUE)
data <- as.data.frame(data)

colnames(data) <- c("FSC_A","FSC_H","SSC_A",
            "CD66b","CD4","CD19",
            "CD14","viability","CD8","CD45","CD3","stainInfo","fileID","sampleID","groupID")


write.table(data,"../../data/All_samples.txt")

####
#clustering
####
library(Rphenograph)
data=read.table("../../data/All_samples.txt")
viable_stained <- data$stainInfo=="stained" & grepl("Live",data$fileID)

channels <- c("FSC_A","SSC_A",
              "CD66b","CD4","CD19",
              "CD14","CD8","CD45","CD3")

cluster_PhenoGraph <- Rphenograph::Rphenograph(data = data[viable_stained,channels], k=60)
cluster_PhenoGraph_vector <- as.vector(membership(cluster_PhenoGraph[[2]]))

out <- data[viable_stained,]
out$clusters <- cluster_PhenoGraph_vector

write.table(out,"../../data/viable_stained_data.txt")

#########
#Filter out CD45 negative and recluster
channels <- c("FSC_A","SSC_A",
              "CD66b","CD4","CD19",
              "CD14","CD8","CD45","CD3")

d <- read.table("../../data/viable_stained_data.txt",header=T)
cluster_median <- aggregate(.~clusters,data = d[,c(channels,"clusters")],median)
#Threshold based on manula inspection
pick <- cluster_median$clusters[which(cluster_median$CD45 > 1)]

data <- d[d$clusters %in% pick,]

cluster_PhenoGraph <- Rphenograph::Rphenograph(data = data[,channels], k=60)
cluster_PhenoGraph <- as.vector(membership(cluster_PhenoGraph[[2]]))

data_tsne <-  Rtsne(X = data[,channels],initial_dims= length(channels),check_duplicates = FALSE)
data_tsne <- data_tsne$Y

data$old_clusters <- data$clusters
data$clusters <- cluster_PhenoGraph
data$tsne1 <- data_tsne[,1]
data$tsne2 <- data_tsne[,2]

write.table(data,"../../data/CD45_data.txt")

