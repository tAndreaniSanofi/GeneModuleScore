# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: MOFA model on RNAseq, Metabolomics, Proteomics, Lipidomics data

rm(list=ls())

library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(reticulate)
library(corrplot)
library(MOFA2)
library(basilisk)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(dplyr)
library(readxl)
library(xlsx)
library(reshape2)

#Start loading the multi omics data
synovium_transcriptome_visit3_rituximab <- read.table("R4RA_synovium_transcriptome_visit3_rituximab_Id_converted.txt",header=T)
blood_transcriptome_visit3_rituximab <- read.table("R4RA_blood_transcriptome_visit3_rituximab_Id_converted.txt",header=T)
lipidomics_visit3_rituximab <- read.table("R4RA_lipidomics_transcriptome_visit3_rituximab_Id_converted.txt",header=T)
metabolomics_visit3_rituximab <- read.csv("R4RA_Pos_Neg_Early_Ritux_Response_Metabolomics data_normalized_Renamed.csv",header=T)
proteomics_visit3_rituximab <- read.table("R4RA_proteomics_transcriptome_visit3_rituximab_Id_converted.txt",header=T)

#subset only the shared among all
shared_patients_across_modalities <- c("CAGLR4RA0932","QMULR4RA0294","LEEDR4RA0740","BASIR4RA1101","LOUVR4RA0956","PAVIR4RA1068","QMULR4RA0162","SOUTR4RA0665","LOUVR4RA0928","QMULR4RA0300","QMULR4RA0490","QMULR4RA0259","QMULR4RA0600","CARDR4RA0862","BASIR4RA0913","QMULR4RA0388","LOUVR4RA0743","NEWCR4RA0993","LOUVR4RA0926","QMULR4RA0078","LOUVR4RA0988","QMULR4RA0331",                                       "LOUVR4RA0628","LISBR4RA0897","QMULR4RA0486","LOUVR4RA0796","GUYSR4RA0860","LEEDR4RA0944","HOMER4RA0973","LOUVR4RA0691","BASIR4RA1143","QMULR4RA0130","WHIPR4RA0627","QMULR4RA0306","QMULR4RA0211","BASIR4RA0646","QMULR4RA0505","QMULR4RA0495","QMULR4RA0094","QMULR4RA0124","QMULR4RA0500","NEWCR4RA0958","NOVAR4RA1185","QMULR4RA0044",
"QMULR4RA0592","QMULR4RA0244","BARCR4RA1056","QMULR4RA0220","SOUTR4RA1080","LISBR4RA0720","CAGLR4RA0617","LISBR4RA1121","QMULR4RA0182","QMULR4RA0296","LOUVR4RA1129","QMULR4RA0511","LOUVR4RA0795","MANCR4RA1103","QMULR4RA0212","QMULR4RA0448")

#subsetting
synovium_transcriptome_visit3_rituximab <- synovium_transcriptome_visit3_rituximab[synovium_transcriptome_visit3_rituximab$PatientId %in% shared_patients_across_modalities, ]
blood_transcriptome_visit3_rituximab <- blood_transcriptome_visit3_rituximab[blood_transcriptome_visit3_rituximab$PatientId %in% shared_patients_across_modalities, ]
lipidomics_visit3_rituximab <- lipidomics_visit3_rituximab[lipidomics_visit3_rituximab$Patient_Id %in% shared_patients_across_modalities, ]
metabolomics_visit3_rituximab <- metabolomics_visit3_rituximab[metabolomics_visit3_rituximab$PatientId %in% shared_patients_across_modalities, ]
proteomics_visit3_rituximab <- proteomics_visit3_rituximab[proteomics_visit3_rituximab$Patient.I.D. %in% shared_patients_across_modalities, ]


# Prepare Transcriptome Synovium dataset --------------------------------------------------
rownames(synovium_transcriptome_visit3_rituximab) <- synovium_transcriptome_visit3_rituximab$PatientId
synovium_transcriptome_visit3_rituximab <- synovium_transcriptome_visit3_rituximab[,-c(1)]
synovium_transcriptome_visit3_rituximab <- t(synovium_transcriptome_visit3_rituximab)
dim(synovium_transcriptome_visit3_rituximab)
synovium_transcriptome_visit3_rituximab[, c(1:60)] <- sapply(synovium_transcriptome_visit3_rituximab[, c(1:60)], as.numeric)
synovMat <- synovium_transcriptome_visit3_rituximab
synovMat <- as.matrix(synovMat)
head(synovMat)

# Data distribution of RNAseq data synovium
boxplot(synovMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Synovium Transcriptome data")


# Prepare Transcriptome Blood dataset --------------------------------------------------
rownames(blood_transcriptome_visit3_rituximab) <- blood_transcriptome_visit3_rituximab$PatientId
blood_transcriptome_visit3_rituximab <- blood_transcriptome_visit3_rituximab[,-c(1)]
blood_transcriptome_visit3_rituximab <- t(blood_transcriptome_visit3_rituximab)
blood_transcriptome_visit3_rituximab[, c(1:60)] <- sapply(blood_transcriptome_visit3_rituximab[, c(1:60)], as.numeric)
bloodMat <- blood_transcriptome_visit3_rituximab
bloodMat <- as.matrix(bloodMat)
head(bloodMat)

# Data distribution of RNAseq data blood
boxplot(bloodMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Blood Transcriptome data")


# Prepare Lipidomics Serum dataset --------------------------------------------------
lipidomics_visit3_rituximab <- lipidomics_visit3_rituximab[,-c(2)]
rownames(lipidomics_visit3_rituximab) <- lipidomics_visit3_rituximab$Patient_Id
lipidomics_visit3_rituximab <- lipidomics_visit3_rituximab[,-c(1)]
lipidomics_visit3_rituximab <- t(lipidomics_visit3_rituximab)
lipidomics_visit3_rituximab[, c(1:60)] <- sapply(lipidomics_visit3_rituximab[, c(1:60)], as.numeric)
lipidomicsMat <- lipidomics_visit3_rituximab
lipidomicsMat <- as.matrix(lipidomicsMat)
head(lipidomicsMat)

# Data distribution of Lipidomics data
boxplot(lipidomicsMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Lipidomics Serum data")


# Prepare Metabolomics Serum dataset --------------------------------------------------
metabolomics_visit3_rituximab <- metabolomics_visit3_rituximab[,-c(2)]
rownames(metabolomics_visit3_rituximab) <- metabolomics_visit3_rituximab$PatientId
metabolomics_visit3_rituximab <- metabolomics_visit3_rituximab[,-c(1)]
metabolomics_visit3_rituximab <- t(metabolomics_visit3_rituximab)
metabolomics_visit3_rituximab[, c(1:60)] <- sapply(metabolomics_visit3_rituximab[, c(1:60)], as.numeric)
metaMat <- metabolomics_visit3_rituximab
metaMat <- as.matrix(metaMat)
dim(metaMat)

# Data distribution of Metabolomics data
boxplot(metaMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Metabolomics Serum data")


# Prepare Proteomics Serum dataset --------------------------------------------------
rownames(proteomics_visit3_rituximab) <- proteomics_visit3_rituximab$Patient.I.D.
proteomics_visit3_rituximab <- proteomics_visit3_rituximab[,-c(1)]
proteomics_visit3_rituximab <- t(proteomics_visit3_rituximab)
proteomics_visit3_rituximab[, c(1:60)] <- sapply(proteomics_visit3_rituximab[, c(1:60)], as.numeric)
proteoMat <- proteomics_visit3_rituximab
proteoMat <- as.matrix(proteoMat)
head(proteoMat)

# Data distribution of Proteomics data
boxplot(proteoMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Proteomics Serum data")



##Load Metada
metadata <-read_xlsx("/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/Pitzalis/scripts/Gene_Module/R4RA/R4RA_full_metadata_20221111.xlsx")
metadata <- as.data.frame(metadata)
metadata_visit3_rituximab <- subset(metadata, Visit == "3" & Randomized.medication == "Rituximab")
metadata_visit3_rituximab <- metadata_visit3_rituximab[metadata_visit3_rituximab$Patient.I.D. %in% shared_patients_across_modalities, ]


library("DESeq2")
dds_syn <- DESeqDataSetFromMatrix(countData = round(synovMat),
                                  colData = metadata_visit3_rituximab,
                                  design = ~ CDAI.response.status.V7
)


# Select top 500 most variable genes
exprMat_syn <- assay(dds_syn)
nTop <- 500
sds_syn <- genefilter::rowSds(exprMat_syn)
exprMat_syn <- exprMat_syn[order(sds_syn, decreasing = T)[1:nTop], ]
exprMat_syn[, c(1:60)] <- sapply(exprMat_syn[, c(1:60)], as.numeric)
str(exprMat_syn)
# Data distribution of metabolomics data
boxplot(exprMat_syn, outline = FALSE, col = "cornflowerblue", main = "Transformed RNAseq Synovium data")
exprMat_syn <- as.matrix(exprMat_syn)
str(exprMat_syn)





# Prepare Blood RNA-seq data --------------------------------------------

library("DESeq2")
dds_blood <- DESeqDataSetFromMatrix(countData = round(bloodMat),
                                    colData = metadata_visit3_rituximab,
                                    design = ~ CDAI.response.status.V7
)


# Select top 500 most variable genes
exprMat_blood <- assay(dds_blood)
nTop <- 500
sds_blood <- genefilter::rowSds(exprMat_blood)
exprMat_blood <- exprMat_blood[order(sds_blood, decreasing = T)[1:nTop], ]
exprMat_blood[, c(1:60)] <- sapply(exprMat_blood[, c(1:60)], as.numeric)
str(exprMat_blood)
# Data distribution of metabolomics data
boxplot(exprMat_blood, outline = FALSE, col = "cornflowerblue", main = "Transformed RNAseq Blood data")
exprMat_blood <- as.matrix(exprMat_blood)
str(exprMat_blood)




# Create the MOFA obejct --------------------------------------------------

# List of data matrix
mofaData <- list(mRNA_synovium = exprMat_syn, mRNA_Blood = exprMat_blood, Lipidomics = lipidomicsMat, Proteins = proteoMat, Metabolites = metaMat)
#mofaData <- list(mRNA_Blood = exprMat_blood, Lipidomics = lipidomicsMat, Proteins = proteoMat, Metabolites = metaMat)
lapply(mofaData, dim)
# Extract samples that appears in both assays
sampleList <- lapply(mofaData, colnames)
useSamples <- intersect(sampleList$mRNA_Blood, sampleList$Metabolites)


# Only keep samples that appears in both assays
f2 <- function(x) x[, useSamples]
mofaData <- lapply(mofaData, f2)

# Build MOFA object -------------------------------------------------------

# Create the MOFA object
MOFAobject <- create_mofa(mofaData)
MOFAobject

# Plot sample size of each modality
fig_data <- plot_data_overview(MOFAobject) +
  scale_fill_brewer(palette = "Set2")

ggsave("./plot/plot_data_overview_vst.png",
       plot = fig_data,
       width = 15, height = 10, units = "cm"
)

# Setup MOFA training parameters ------------------------------------------

# List data options
DataOptions <- get_default_data_options(MOFAobject)
DataOptions

# Define model options
ModelOptions <- get_default_model_options(MOFAobject)
ModelOptions$num_factors <- 10 # Number of factors
ModelOptions

# Define training options
TrainOptions <- get_default_training_options(MOFAobject)
TrainOptions$drop_factor_threshold <- 0.02
TrainOptions$convergence_mode <- "slow"
TrainOptions$seed <- 1234
TrainOptions$verbose <- TRUE
TrainOptions

# Prepare the MOFA object -------------------------------------------------

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = DataOptions,
                           model_options = ModelOptions,
                           training_options = TrainOptions,
                           
)

# Train the MOFA model ----------------------------------------------------

# Define file path for saving the trained model
outfile <- file.path("./result/model_Blood_Synovium_RNA-seq_vst.hdf5")

# Train the MOFA model
#go to the bash  and :
#1) create a conda environment e.g. conda create --name mofa_train python==3.7.5
#2) conda activate mofa_train
#3) pip install mofapy2
#4) conda deactivate
#now the tools are installed and you can train the model with run_mofa

library(basilisk)
useBasiliskEnv("/cloud-home/I0439277/.magellan/conda/envs/mofa_train")
MOFAobject <- run_mofa(MOFAobject, outfile, use_basilisk = FALSE)

# Variance explained by MOFA for each omic data ---------------------------
plot_factor_cor(MOFAobject)
# Load model
model <- load_model("./result/model_Blood_Synovium_RNA-seq_vst.hdf5")
model

# Calculate the variance explained (R2) per factor in each view
calculate_variance_explained(model)

# Plot variance explained
fig_v1 <- plot_variance_explained(model) +
  scale_fill_viridis_c(breaks = c(5, 10, 15, 20, 25))

ggsave("./plot/plot_var_per_factor_viridis_vst.png", plot = fig_v1, width = 15, height = 10, units = "cm")

fig_v2 <- plot_variance_explained(MOFAobject, plot_total = T)[[2]]
ggsave("./plot/plot_var_per_view_vst.png", plot = fig_v2, width = 15, height = 10, units = "cm")

# Sanity check: Factors should be uncorrelated ----------------------------

# Get factor data
Z <- get_factors(model)

# Compute correlation
r <- abs(cor(
  x = do.call(rbind, Z), y = do.call(rbind, Z),
  method = "pearson", use = "complete.obs"
))

# Corrplot
png(
  filename = "./plot/factor_corrPlot_vst.png",
  width = 600, height = 600, units = "px"
)
corrplot(r,
         method = "color", type = "upper", order = "hclust",
         col = rev(brewer.pal(n = 8, name = "RdYlBu")),
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, tl.cex = 1, # Text label color, rotation and size
         diag = FALSE
)
dev.off()

samples_shared <- MOFAobject@samples_metadata$sample
#devtools::install_github("koenderks/aRtsy")
#install.packages("gtable")
plot_factor_cor(MOFAobject)
plot_variance_explained(MOFAobject, max_r2=15)
plot_variance_explained(MOFAobject, plot_total = T)[[2]]


######################
##Interpretation output
######################
##Interpretation output
metadata <-read_xlsx("/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/Pitzalis/scripts/Gene_Module/R4RA/R4RA_full_metadata_20221111.xlsx")
metadata <- as.data.frame(metadata)
metadata_visit3_rituximab <- subset(metadata, Visit == "3" & Randomized.medication == "Rituximab")
metadata_visit3_rituximab <- metadata_visit3_rituximab[metadata_visit3_rituximab$Patient.I.D. %in% useSamples, ]
metadata_visit3_rituximab <- metadata_visit3_rituximab[,-c(1,2)]
colnames(metadata_visit3_rituximab)[1] <- 'sample'
samples_metadata(MOFAobject) <- metadata_visit3_rituximab
dim(MOFAobject@samples_metadata)
plot_data_overview(MOFAobject)
plot_factor_cor(MOFAobject)
plot_variance_explained(MOFAobject, max_r2=15)
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
colnames(MOFAobject@data$Metabolites$group1)
colnames(MOFAobject@data$mRNA_synovium$group1)
colnames(MOFAobject@data$mRNA_Blood$group1)==colnames(MOFAobject@data$mRNA_synovium$group1)
correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("CDAI.response.status.V7","Pathotype","CCP.Visit1","RF.Visit1","RF_ACPA_status.Visit1","CD3","CD20","CD68L","CD68SL","CD21","CD138","Tender.Joint.Count","Swollen.Joint.Count","VAS.Global","Creatinine","Joint.stiffness","Gender"), 
                                  plot="log_pval")

plot_factor(MOFAobject, 
            factor =4, 
            color_by = "Factor4"
)

plot_weights(MOFAobject,
             view = "Metabolites",
             factors = 4,
             nfeatures = 20,     # Top number of features to highlight
             scale = F,           # Scale weights from -1 to 1
)




plot_top_weights(MOFAobject,
                 view = "mRNA_synovium",
                 factor =4,
                 nfeatures = 20,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "mRNA_Blood",
                 factor =4,
                 nfeatures = 20,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "mRNA_Blood",
                 factor =5,
                 nfeatures = 20,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)


plot_top_weights(MOFAobject,
                 view = "Lipidomics",
                 factor =3,
                 nfeatures = 20,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)


plot_factor(MOFAobject, 
            factors = 3, 
            color_by = "Proteins",
            add_violin = TRUE,
            dodge = TRUE
)



plot_factor(MOFAobject, 
            factors = 3, 
            color_by = "CDAI.response.status.V7",
            add_violin = TRUE,
            dodge = TRUE
)

plot_factor(MOFAobject, 
            factors = 3, 
            color_by = "Creatinine",
            add_violin = TRUE,
            dodge = TRUE
)


plot_factor(MOFAobject, 
            factors = 3, 
            color_by = "RF_ACPA_status.Visit1",
            add_violin = TRUE,
            dodge = TRUE
)


plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =4,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
) + labs(y="RNA expression Blood")





# Positive CDAI
plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =4,  
                  features = 20,
                  sign = "positive",
                  color_by = "Gender"
) + labs(y="RNA expression Blood")


plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =5,  
                  features = 20,
                  sign = "positive",
                  color_by = "CDAI.response.status.V7"
)

plot_data_scatter(MOFAobject, 
                  view = "mRNA_synovium",
                  factor =3,  
                  features = 20,
                  sign = "positive",
                  color_by = "CDAI.response.status.V7"
) 


plot_data_scatter(MOFAobject, 
                  view = "mRNA_synovium",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)+ labs(y="RNA expression Synovium")


plot_data_scatter(MOFAobject, 
                  view = "Proteins",
                  factor =3,  
                  features = 20,
                  sign = "positive",
                  color_by = "CDAI.response.status.V7"
)


plot_data_scatter(MOFAobject, 
                  view = "Lipidomics",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
) + labs(y="Lipidomics Expression")

plot_data_scatter(MOFAobject, 
                  view = "Lipidomics",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)

plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =4,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
) + labs(y="Matabolites Expression")


plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =3,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)





### negative CDAI
plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =3,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)


plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)

plot_data_scatter(MOFAobject, 
                  view = "mRNA_synovium",
                  factor =3,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
) 


plot_data_scatter(MOFAobject, 
                  view = "mRNA_synovium",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)


plot_data_scatter(MOFAobject, 
                  view = "Proteins",
                  factor =3,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)


plot_data_scatter(MOFAobject, 
                  view = "Proteins",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)

plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =3,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
) 


plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)





# Positive RF_ACPA_status.Visit1
plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =3,  
                  features = 20,
                  sign = "positive",
                  color_by = "RF_ACPA_status.Visit1"
) 


plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =5,  
                  features = 20,
                  sign = "positive",
                  color_by = "RF_ACPA_status.Visit1"
)

plot_data_scatter(MOFAobject, 
                  view = "mRNA_synovium",
                  factor =3,  
                  features = 20,
                  sign = "positive",
                  color_by = "RF_ACPA_status.Visit1"
) 


plot_data_scatter(MOFAobject, 
                  view = "mRNA_synovium",
                  factor =5,  
                  features = 20,
                  sign = "positive",
                  color_by = "RF_ACPA_status.Visit1"
)


plot_data_scatter(MOFAobject, 
                  view = "Proteins",
                  factor =3,  
                  features = 20,
                  sign = "positive",
                  color_by = "RF_ACPA_status.Visit1"
)


plot_data_scatter(MOFAobject, 
                  view = "Proteins",
                  factor =5,  
                  features = 20,
                  sign = "positive",
                  color_by = "RF_ACPA_status.Visit1"
)

plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =3,  
                  features = 20,
                  sign = "positive",
                  color_by = "RF_ACPA_status.Visit1"
) 


plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =5,  
                  features = 20,
                  sign = "positive",
                  color_by = "RF_ACPA_status.Visit1"
)





### negative
plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =3,  
                  features = 20,
                  sign = "negative",
                  color_by = "RF_ACPA_status.Visit1"
)


plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "RF_ACPA_status.Visit1"
)

plot_data_scatter(MOFAobject, 
                  view = "mRNA_synovium",
                  factor =3,  
                  features = 20,
                  sign = "negative",
                  color_by = "RF_ACPA_status.Visit1"
) 


plot_data_scatter(MOFAobject, 
                  view = "mRNA_synovium",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "RF_ACPA_status.Visit1"
)


plot_data_scatter(MOFAobject, 
                  view = "Proteins",
                  factor =3,  
                  features = 20,
                  sign = "negative",
                  color_by = "RF_ACPA_status.Visit1"
)


plot_data_scatter(MOFAobject, 
                  view = "Proteins",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "RF_ACPA_status.Visit1"
)

plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =3,  
                  features = 20,
                  sign = "negative",
                  color_by = "RF_ACPA_status.Visit1"
) 


plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "RF_ACPA_status.Visit1"
)

