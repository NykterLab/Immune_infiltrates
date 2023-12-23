## Bulk-RNA sequencing data analysis

## deconvolution of immune cell type proportions,
## calculating immune cluster activities,
## expressions of immune geneset
## and correlations 



# TPM normalization function
tpm <- function(counts,len) {
  x <- counts/(len) # gene length in kilobases
  return(t(t(x)*1e6/colSums(x)))
}

# load expression matrix of in-house and TCGA samples
counts = read.table()

#load reference cell types from Joyce et al.
ref_file <- "BrainTIME_RawCounts.csv"
ref <- read.csv(ref_file, header=TRUE)
clin_ref = read.csv("BrainTIME_ClinicalAnnotation.csv", header=TRUE)
clin_ref = clin_ref[!is.na(clin_ref$Histology), ]
patients = clin_ref[clin_ref$Histology == "Glioblastoma", ]
patients = as.vector(patients$Patient)

patient_cols = c(1)
for (i in 1:ncol(ref)){
  for (j in 1:length(patients)){
    if (grepl(as.character(patients[j]), as.character(colnames(ref)[i]), fixed = TRUE)){
      patient_cols[length(patient_cols)+1] = i
    }
  }
}
ref_gbm = ref[,c(patient_cols)]



# load lengths of the genes in expression matrix
lengths = read.table()
lengths = as.vector(lengths)

lengths_ref = read.table()
lengths_ref = as.vector(lengths_ref)

# TPM NORMALIZATION
tpm_normalized_counts <- tpm(counts, lengths)

#TPM NORM OF REFERENCES
rownames(ref_gbm) = ref_gbm[,1]
ref_gbm = ref_gbm[,-1]
tpm_ref <- tpm(ref_gbm, lengths_ref)

#merging the data
tpm_all = merge(tpm_normalized_counts, tpm_ref, by=0)
rownames(tpm_all) = tpm_all[,1]
tpm_all = tpm_all[,-1]

# QUANTILE NORMALIZATION: counts
library(preprocessCore)
quantile_counts_joined <- normalize.quantiles(as.matrix(tpm_all), copy=FALSE)
quantile_counts_gbm <- normalize.quantiles(as.matrix(tpm_normalized_counts), copy=FALSE)
quantile_counts_ref <- normalize.quantiles(as.matrix(tpm_ref), copy=FALSE)

# Log2 transformation: ref + counts
quantile_counts_log2_joined <- log2(quantile_counts_joined+1)
quantile_counts_log2_gbm <- log2(quantile_counts_gbm+1)
quantile_counts_log2_ref <- log2(quantile_counts_ref+1)
quantile_counts_log2 <- quantile_counts_log2_gbm

#Samples (12 in-house samples and 156 TCGA samples)
all_gbm_samples <- quantile_counts_log2[,1:168]
samples <- quantile_counts_log2[,1:168]
dim(samples)

all_gbm_samples_joined <- quantile_counts_log2_joined[,1:168]
samples_joined <- quantile_counts_log2_joined[,1:168]
dim(samples_joined)
ref_normalized_joined <- quantile_counts_log2_joined[,169:236]

# DECONVOLUTION

# Limit matrices to immune cluster genes
geneset <- read.table("deconvolution_genes.txt", header=FALSE)
geneset <- as.character(geneset$V1)

samples2 <- all_gbm_samples_joined[rownames(all_gbm_samples_joined) %in% geneset,]
samples2 <- na.omit(samples2)
reference = ref_normalized_joined
reference2 <- reference[rownames(reference) %in% geneset,]
reference2 <- na.omit(reference2)

celltypes = c("CD45-", "Microglia", "MDM", "Neutrophils", "CD4+ T-cells", "CD8+ T-cells")

reference3 <- matrix(nrow=2515, ncol=6)

cd45n = reference2[,c(1:11,44:48)]
mg = reference2[,c(12:21,49:52)]
mdm = reference2[,c(22:29,53:56)]
neu = reference2[,c(30:35,57:61)]
cd4 = reference2[,c(36:42,62:64)]
cd8 = reference2[,c(43,65:68)]

cd45n_med = apply(cd45n, 1, median)
mg_med = apply(mg, 1, median)
mdm_med = apply(mdm, 1, median)
neu_med = apply(neu, 1, median)
cd4_med = apply(cd4, 1, median)
cd8_med = apply(cd8, 1, median)

reference3[,1] = cd45n_med
reference3[,2] = mg_med
reference3[,3] = mdm_med
reference3[,4] = neu_med
reference3[,5] = cd4_med
reference3[,6] = cd8_med

rownames(reference3) = rownames(reference2)
colnames(reference3) = celltypes

# Function for regression analysis
regressionAnalysis = function(mixture, reference_cells, alpha = NULL) {
  # Dependencies 
  if (!require(glmnet)) {
    install.packages('glmnet')
  }
  
  library(glmnet)
  
  # Initialize results matrix
  res = matrix(nrow = length(colnames(reference_cells)),
               ncol = length(colnames(mixture)))
  colnames(res) = colnames(mixture)
  rownames(res) = colnames(reference_cells)
  
  # Set alpha
  if (is.null(alpha)) {
    ALPHA = 0.5
  } else {
    ALPHA = alpha
  }
  
  # Run regression analysis separateky for each sample
  for (s in 1:length(colnames(mixture))) {
    
    # Choose lambda using cross-validation
    LAMBDA = cv.glmnet(x = reference_cells, y = mixture[, s], alpha = ALPHA)$lambda.1se
    
    # Calculate coeffients
    res[, s] = coef(glmnet(x = reference_cells, y = mixture[, s], alpha = ALPHA), 
                    s = LAMBDA)[-1, 1]
  }
  
  return(res)
}

deconv_res <- regressionAnalysis(as.matrix(samples2), as.matrix(reference3), 0.25)


# IMMUNE CLUSTER ACTIVITY

# Load immune cluster genes from Luoto et al. 
cluster_all <- read.xlsx("193230_3_supp_4816818_pf5jrc.xlsx", startRow = 3)

dat <- matrix(ncol = 168, nrow = 8)
nam <- c("1 Macrophages and T cell response", 
         "2 Negative regulation of lymphocyte response", 
         "3 Leukocyte migration", 
         "4 Humoral response and lymphocytes",
         "5 Antigen presentation and interferon response", 
         "6 Leukocyte differentiation and chemotaxis",
         "7 Gamma delta T cells", 
         "8 Negative regulation of T-cell activation, PD-L1")
rownames(dat) <- nam

for (j in 1:8){
  cluster <- cluster_all[,j]
  cluster <- na.omit(cluster)
  
  indices <- match(as.character(cluster), trimws(as.character(rownames(quantile_counts_log2)), which="left"))
  indices <- indices[!is.na(indices)]
  cluster_counts <- quantile_counts_log2[indices,]
  varianceGenes = apply(cluster_counts, 1, var)
  # Level of the expression 
  filt = cluster_counts[which(varianceGenes > 0.05),]
  Correlations <- cor(t(as.matrix(filt)))
  
  # Count how many possitively correlating genes each gene has
  max_i <- 0
  max_pos <- 0
  for ( gene in 1:dim(Correlations)[1]) {
    row <- Correlations[gene,]
    n_pos <- length(row[which(row[] > 0)])
    if ( n_pos > max_pos ) {
      max_pos <- n_pos
      max_i <- gene
    }
  }
  
  # The the biggest group of positively correlating genes
  biggest <- Correlations[max_i,]
  positive <- biggest[which(biggest[] > 0)]
  length(positive)
  
  # Rescale counts to 0...1 

  genes <- trimws(as.character(names(positive)), which="left")
  countstrim <- quantile_counts_log2[match(genes, 
                                           trimws(as.character(rownames(quantile_counts_log2)), 
                                                  which="left")),]
  backup <- countstrim
  print(dim(countstrim))
  library("scales")
  for ( row in 1:dim(countstrim)[1]) {
    i <- order(countstrim[row,])
    countstrim[row,i[1:4]] <- 0
    countstrim[row,i[(length(i)-3):length(i)]] <- 1
    countstrim[row,i[5:(length(i)-4)]] <- rescale(countstrim[row,i[5:(length(i)-4)]])
    # 5 smallest to 0, 5 biggest to 1 and everything else between there
  }
  
  # Take median of these rescaled genes
  activity <- apply(countstrim, 2, median)
  
  
  dat[j,] <- activity
  colnames(dat) <- names(activity)
}

activity_tab = dat

#CORRELATIONS

#activity vs deconv pearson
cor_matrix = matrix(nrow = nrow(deconv_res), ncol = nrow(activity_tab))
rownames(cor_matrix) = rownames(deconv_res)
colnames(cor_matrix) = rownames(activity_tab)

for (i in 1:nrow(deconv_res)){
  for (j in 1:nrow(activity_tab)){
    cor_matrix[i,j] = cor(deconv_res[i,], activity_tab[j,], method = "pearson")
  }
}

cor_matrix_activity_deconv_pearson = cor_matrix

#activity vs deconv spearman
cor_matrix = matrix(nrow = nrow(deconv_res), ncol = nrow(activity_tab))
rownames(cor_matrix) = rownames(deconv_res)
colnames(cor_matrix) = rownames(activity_tab)

for (i in 1:nrow(deconv_res)){
  for (j in 1:nrow(activity_tab)){
    cor_matrix[i,j] = cor(deconv_res[i,], activity_tab[j,], method = "spearman")
  }
}

cor_matrix_activity_deconv_spearman = cor_matrix

#deconv pearson
cor_deconv_pearson = cor(t(deconv_res), method = "pearson")
#deconv spearman 
cor_deconv_spearman = cor(t(deconv_res), method = "spearman")

#activity pearson 
cor_activity_pearson = cor(t(activity_tab), method = "pearson")
#activity spearman 
cor_activity_spearman = cor(t(activity_tab), method = "spearman")



#EXPRESSION OF IMMUNE GENE SET

#Load immune gene set from Danaher et al. 
immune_genes_table <- read.xlsx("40425_2017_215_MOESM1_ESM.xlsx", sheet = 4)

immune_genes_table <- immune_genes_table[c(10,1:9,13:60,11:12),]
immune_genes_table <- immune_genes_table[-7,]
immune_genes <- as.vector(immune_genes_table[,1])

immune_samples <- all_gbm_samples[match(immune_genes, trimws(as.character(rownames(all_gbm_samples)), which ="left")),]
rownames(immune_samples) <- trimws(as.character(rownames(immune_samples)), which ="left")

idx <- rowSums( immune_samples >= 1 ) >= 8

immune_samples_filtered = immune_samples[idx,]


m = apply(immune_samples, 1, median)
median_norm_immune_samples = immune_samples - m



