# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: MOFA model on RNAseq and metabolomics data
# Data: DESeq2 object after variance stabilizing transformation (VST)
# Data: Metabolon peak area data

# Reference: https://huber-group-embl.github.io/mofaCLL/analysisProcedure.html
# Updated tutorial: https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html

# Set environment variables -----------------------------------------------

#Sys.setenv(PATH = paste("/cloud-home/I0442220/.cache/basilisk/1.2.1/MOFA2-1.0.1/mofa_env",
#                        Sys.getenv()["PATH"],
#                        sep = ";"
#))

# Load library ------------------------------------------------------------


#Sys.setenv(PATH = paste("/cloud-home/I0439277/.cache/basilisk/1.2.1/MOFA2-1.0.1/mofa_env/",
#                        Sys.getenv()["PATH"],
#                        sep = ";"
#))

#https://stackoverflow.com/questions/75772757/r-for-data-science-chapter-2-2-2-already-getting-error-messages-before-ive-star
#unlink("/cloud-home/I0439277/R/x86_64-pc-linux-gnu-library/4.0.0/00LOCK-rlang", recursive = TRUE)
#remove.packages("rlang")
#remove.packages("ggplot2")
#install.packages("ggplot2")
#install.packages("rlang")

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

#Format Proteomics Data
setwd("cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/Pitzalis/R4RA/test_multiomics/Test/")
proteomics <- read.csv("R4RA_olinkdata_withmeta_Proteomics.csv")
proteomics <- subset(proteomics, Visit == "3" & Randomised.medication == "Tocilizumab")
proteomics1 <- proteomics[ ,c(16,4,13)]
test1 <- dcast(proteomics1,Patient.I.D.~OlinkID)
write.table(test1,"R4RA_olinkdata_withmeta_Proteomics_Toci_Visit3.txt",quote=F,sep = "\t",col.names = T,row.names = F)

proteomics_visit3_toci <- read.table("R4RA_olinkdata_withmeta_Proteomics_Toci_Visit3.txt",header = T)
proteomics_visit3_toci <- proteomics_visit3_toci[order(proteomics_visit3_toci$Patient.I.D.),]

#Metabolomics Neg Mode
metabolomics <- read.csv("../../tocilizumab/R4RA_Pos_Neg_Early_Toci_Response_Metabolomics data_normalized_Renamed.csv",header=T)
metabolomics_visit3_toci <- metabolomics[metabolomics$Patient_Id %in% proteomics_visit3_toci$Patient.I.D., ]
metabolomics_visit3_toci <- metabolomics_visit3_toci[order(metabolomics_visit3_toci$Patient_Id),]

#Lipidomics
lipidomics <- read.csv("../../tocilizumab/Early_Toci_Response_data_normalized_Renamed_Lipidomics.csv",header=T)
lipidomics_visit3_toci <- lipidomics[lipidomics$Patient_Id %in% proteomics_visit3_toci$Patient.I.D., ]
lipidomics_visit3_toci <- lipidomics_visit3_toci[order(lipidomics_visit3_toci$Patient_Id),]



######### Load Metadata
metadata <-read_xlsx("/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/Pitzalis/scripts/Gene_Module/R4RA/R4RA_full_metadata_20221111.xlsx")
metadata <- as.data.frame(metadata)
metadata_visit3_toci <- subset(metadata, Visit == "3" & Randomized.medication == "Tocilizumab")
patient_id <- metadata[,c(2,3,4)]
metadata_visit3_toci <- metadata_visit3_toci[metadata_visit3_toci$Patient.I.D. %in% proteomics_visit3_toci$Patient.I.D., ]
metadata_visit3_toci <- metadata_visit3_toci[order(metadata_visit3_toci$Patient.I.D.),]
Id <-  metadata_visit3_toci[,c(2,3)]
Id_conversion <- Id[Id$QiagenID_Synovium %in% synovium_transcriptome_t$Patient_Id, ]
Id_conversion <- Id_conversion[order(Id_conversion$QiagenID_Synovium),]


############## Load Gene expression Synovium
synovium_transcriptome <- read.table("Pitzalis_R4RA_Synovium_All_samples.RnaSeq_Genes.40Outliers_removed.Count_20221205.txt",header=T)
synovium_transcriptome <- data.frame(synovium_transcriptome, row.names = 1)
synovium_transcriptome_t <- as.data.frame(t(synovium_transcriptome))
synovium_transcriptome_t$Patient_Id <- rownames(synovium_transcriptome_t)
synovium_transcriptome_t <- synovium_transcriptome_t[synovium_transcriptome_t$Patient_Id %in% metadata_visit3_toci$QiagenID_Synovium, ]
synovium_transcriptome_t <- synovium_transcriptome_t[order(synovium_transcriptome_t$Patient_Id),]
synovium_transcriptome_t <- cbind(Id_conversion$Patient.I.D.,synovium_transcriptome_t)
synovium_transcriptome_t <- as.data.frame(synovium_transcriptome_t)
colnames(synovium_transcriptome_t)[1] = "PatientId"
synovium_transcriptome_t <- synovium_transcriptome_t[,-60701]
synovium_transcriptome_visit3_toci <- synovium_transcriptome_t

################ Load Gene expression Blood
blood_transcriptome <- read.table("Pitzalis_R4RA_Blood_All_samples.RnaSeq_Genes.7Outliers_removed.Count.vst_20221121.txt",header=T)
blood_transcriptome_t <- as.data.frame(t(blood_transcriptome))
names(blood_transcriptome_t) <- blood_transcriptome_t[1,]
blood_transcriptome_t <- blood_transcriptome_t[-1,]
blood_transcriptome_t$Patient_Id <- rownames(blood_transcriptome_t)
blood_transcriptome_t <- blood_transcriptome_t[order(blood_transcriptome_t$Patient_Id),]
blood_transcriptome_t <- blood_transcriptome_t[blood_transcriptome_t$Patient_Id %in% metadata_visit3_toci$QiagenID_Synovium, ]
Id_conversion_blood <- Id[Id$QiagenID_Synovium %in% blood_transcriptome_t$Patient_Id, ]
Id_conversion_blood <- Id_conversion_blood[order(Id_conversion_blood$QiagenID_Synovium),]
blood_transcriptome_t <- cbind(Id_conversion_blood$Patient.I.D.,blood_transcriptome_t)
blood_transcriptome_t <- as.data.frame(blood_transcriptome_t)
colnames(blood_transcriptome_t)[1] = "PatientId"
blood_transcriptome_t <- blood_transcriptome_t[,-60701]
blood_transcriptome_visit3_toci <- blood_transcriptome_t
head(blood_transcriptome_visit3_toci[1:3,1:3])

#select common patients
blood_transcriptome_visit3_toci$PatientId
dim(blood_transcriptome_visit3_toci)
synovium_transcriptome_visit3_toci$PatientId
dim(synovium_transcriptome_visit3_toci)
lipidomics_visit3_toci$Patient_Id
dim(lipidomics_visit3_toci)
metabolomics_visit3_rituximab$PatientId
dim(metabolomics_visit3_toci)
proteomics_visit3_toci$Patient.I.D.
dim(proteomics_visit3_toci)
metadata_visit3_toci$Patient.I.D.
dim(metadata_visit3_toci)

dim(blood_transcriptome_visit3_toci)
dim(synovium_transcriptome_visit3_toci)
dim(lipidomics_visit3_toci)
dim(metabolomics_visit3_toci)
dim(proteomics_visit3_toci)

metabolomics_visit3_toci$Patient_Id

getwd()
write.table(synovium_transcriptome_visit3_toci,"../../tocilizumab/R4RA_synovium_transcriptome_visit3_toci_Id_converted.txt",quote=F,sep = "\t",col.names = T,row.names = F)
write.table(blood_transcriptome_visit3_toci,"../../tocilizumab/R4RA_blood_transcriptome_visit3_toci_Id_converted.txt",quote=F,sep = "\t",col.names = T,row.names = F)
write.table(lipidomics_visit3_toci,"../../tocilizumab/R4RA_lipidomics_visit3_toci_Id_converted.txt",quote=F,sep = "\t",col.names = T,row.names = F)
write.table(metabolomics_visit3_toci,"../../tocilizumab/R4RA_metabolomics_visit3_toci_Id_converted.txt",quote=F,sep = "\t",col.names = T,row.names = F)
write.table(proteomics_visit3_toci,"../../tocilizumab/R4RA_proteomics_visit3_toci_Id_converted.txt",quote=F,sep = "\t",col.names = T,row.names = F)

# shared : c(NOVAR4RA1089 LOUVR4RA0815 SENDR4RA1124 NEWCR4RA0801 LOUVR4RA1038 QMULR4RA0400 QMULR4RA0260 QMULR4RA0160 QMULR4RA0535 BARCR4RA0820 LISBR4RA1201 LOUVR4RA0605 QMULR4RA0196 NEWCR4RA1014 QMULR4RA0474 QMULR4RA0123 MANCR4RA0881 CARDR4RA0807 QMULR4RA0042 LISBR4RA0802 QMULR4RA0565 QMULR4RA0221 QMULR4RA0567 NOVAR4RA1091 QMULR4RA0597 LISBR4RA0636 CAGLR4RA0667 BARCR4RA0616 SOUTR4RA1084 QMULR4RA0128 NEWCR4RA0715 BARCR4RA0699 LOUVR4RA0888 LOUVR4RA1009 QMULR4RA0048 GUYSR4RA1198 QMULR4RA0401 LOUVR4RA1179 LOUVR4RA1042 SOUTR4RA0677 QMULR4RA0102 NOVAR4RA0768 CARDR4RA1186 QMULR4RA0315 HOMER4RA0804 CARDR4RA0676 BASIR4RA1157 QMULR4RA0373 LOUVR4RA0889 LOUVR4RA1083 HOMER4RA0837 QMULR4RA0197 SOUTR4RA0921 LOUVR4RA0713 LEUVR4RA1018 QMULR4RA0142 LOUVR4RA0975 QMULR4RA0292 LOUVR4RA0915)

#subset only the shared among all
shared_patients_across_modalities <- c("NOVAR4RA1089","LOUVR4RA0815","SENDR4RA1124"," NEWCR4RA0801","LOUVR4RA1038","QMULR4RA0400","QMULR4RA0260","QMULR4RA0160","QMULR4RA0535","BARCR4RA0820","LISBR4RA1201","LOUVR4RA0605","QMULR4RA0196","NEWCR4RA1014","QMULR4RA0474","QMULR4RA0123","MANCR4RA0881","CARDR4RA0807","QMULR4RA0042","LISBR4RA0802","QMULR4RA0565","QMULR4RA0221","QMULR4RA0567","NOVAR4RA1091","QMULR4RA0597","LISBR4RA0636","CAGLR4RA0667","BARCR4RA0616","SOUTR4RA1084","QMULR4RA0128","NEWCR4RA0715","BARCR4RA0699","LOUVR4RA0888","LOUVR4RA1009","QMULR4RA0048","GUYSR4RA1198","QMULR4RA0401","LOUVR4RA1179","LOUVR4RA1042","SOUTR4RA0677","QMULR4RA0102","NOVAR4RA0768","CARDR4RA1186","QMULR4RA0315","HOMER4RA0804","CARDR4RA0676","BASIR4RA1157","QMULR4RA0373","LOUVR4RA0889","LOUVR4RA1083","HOMER4RA0837","QMULR4RA0197","SOUTR4RA0921","LOUVR4RA0713","LEUVR4RA1018","QMULR4RA0142","LOUVR4RA0975","QMULR4RA0292","LOUVR4RA0915")



#load modalities
synovium_transcriptome_visit3_toci <- read.table("../../tocilizumab/R4RA_synovium_transcriptome_visit3_toci_Id_converted.txt",header=T)
blood_transcriptome_visit3_toci <- read.table("../../tocilizumab/R4RA_blood_transcriptome_visit3_toci_Id_converted.txt",header=T)
lipidomics_visit3_toci <- read.table("../../tocilizumab/R4RA_lipidomics_visit3_toci_Id_converted.txt",header=T)
metabolomics_visit3_toci <- read.csv("../../tocilizumab/R4RA_Pos_Neg_Early_Toci_Response_Metabolomics data_normalized_Renamed.csv",header=T)
proteomics_visit3_toci <- read.table("../../tocilizumab/R4RA_proteomics_visit3_toci_Id_converted.txt",header=T)

#subsetting
synovium_transcriptome_visit3_toci <- synovium_transcriptome_visit3_toci[synovium_transcriptome_visit3_toci$PatientId %in% shared_patients_across_modalities, ]
blood_transcriptome_visit3_toci <- blood_transcriptome_visit3_toci[blood_transcriptome_visit3_toci$PatientId %in% shared_patients_across_modalities, ]
lipidomics_visit3_toci <- lipidomics_visit3_toci[lipidomics_visit3_toci$Patient_Id %in% shared_patients_across_modalities, ]
metabolomics_visit3_toci <- metabolomics_visit3_toci[metabolomics_visit3_toci$Patient_Id %in% shared_patients_across_modalities, ]
proteomics_visit3_toci <- proteomics_visit3_toci[proteomics_visit3_toci$Patient.I.D. %in% shared_patients_across_modalities, ]


# Prepare Transcriptome Synovium dataset --------------------------------------------------
rownames(synovium_transcriptome_visit3_toci) <- synovium_transcriptome_visit3_toci$PatientId
dim(blood_transcriptome_visit3_toci)
synovium_transcriptome_visit3_toci <- synovium_transcriptome_visit3_toci[,-c(1)]
synovium_transcriptome_visit3_toci <- t(synovium_transcriptome_visit3_toci)
dim(synovium_transcriptome_visit3_toci)
synovium_transcriptome_visit3_toci[, c(1:58)] <- sapply(synovium_transcriptome_visit3_toci[, c(1:58)], as.numeric)
synovMat <- synovium_transcriptome_visit3_toci
synovMat <- as.matrix(synovMat)
head(synovMat)

# Data distribution of RNAseq data
boxplot(synovMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Synovium Transcriptome data")


# Prepare Transcriptome Blood dataset --------------------------------------------------

rownames(blood_transcriptome_visit3_toci) <- blood_transcriptome_visit3_toci$PatientId
blood_transcriptome_visit3_toci <- blood_transcriptome_visit3_toci[,-c(1)]
blood_transcriptome_visit3_toci <- t(blood_transcriptome_visit3_toci)
blood_transcriptome_visit3_toci[,c(1:58)] <- sapply(blood_transcriptome_visit3_toci[, c(1:58)], as.numeric)
bloodMat <- blood_transcriptome_visit3_toci
bloodMat<-as.data.frame(bloodMat)
i <- 1:58           
bloodMat[ , i] <- apply(bloodMat[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))


bloodMat <-as.matrix(bloodMat)
head(bloodMat)

# Data distribution of RNAseq data
boxplot(bloodMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Blood Transcriptome data")


# Prepare Lipidomics Serum dataset --------------------------------------------------
rownames(lipidomics_visit3_toci) <- lipidomics_visit3_toci$Patient_Id
head(lipidomics_visit3_toci)
lipidomics_visit3_toci <- lipidomics_visit3_toci[,-c(1)]
lipidomics_visit3_toci <- t(lipidomics_visit3_toci)
dim(lipidomics_visit3_toci)
lipidomics_visit3_toci[, c(1:58)] <- sapply(lipidomics_visit3_toci[, c(1:58)], as.numeric)
lipidomicsMat <- lipidomics_visit3_toci
lipidomicsMat <- as.matrix(lipidomicsMat)
dim(lipidomicsMat)

# Data distribution of RNAseq data
boxplot(lipidomicsMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Lipidomics Serum data")


# Prepare Metabolomics Serum dataset --------------------------------------------------
rownames(metabolomics_visit3_toci) <- metabolomics_visit3_toci$Patient_Id
metabolomics_visit3_toci <- metabolomics_visit3_toci[,-c(1)]
metabolomics_visit3_toci <- t(metabolomics_visit3_toci)
metabolomics_visit3_toci[, c(1:58)] <- sapply(metabolomics_visit3_toci[, c(1:58)], as.numeric)
metaMat <- metabolomics_visit3_toci
metaMat <- as.matrix(metaMat)
dim(metabolomics_visit3_toci)
head(metaMat)
# Data distribution of RNAseq data
boxplot(metaMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Metabolomics Serum data")


# Prepare Proteomics Serum dataset --------------------------------------------------
rownames(proteomics_visit3_toci) <- proteomics_visit3_toci$Patient.I.D.
proteomics_visit3_toci <- proteomics_visit3_toci[,-c(1)]
proteomics_visit3_toci <- t(proteomics_visit3_toci)
proteomics_visit3_toci[, c(1:58)] <- sapply(proteomics_visit3_toci[, c(1:58)], as.numeric)
proteoMat <- proteomics_visit3_toci
proteoMat <- as.matrix(proteoMat)
head(proteoMat)

# Data distribution of RNAseq data
boxplot(proteoMat, outline = FALSE, col = "cornflowerblue", main = "Transformed Proteomics Serum data")



##Load Metada
metadata <-read_xlsx("/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/Pitzalis/scripts/Gene_Module/R4RA/R4RA_full_metadata_20221111.xlsx")
metadata <- as.data.frame(metadata)
metadata_visit3_toci <- subset(metadata, Visit == "3" & Randomized.medication == "Tocilizumab")
metadata_visit3_toci <- metadata_visit3_toci[metadata_visit3_toci$Patient.I.D. %in% shared_patients_across_modalities, ]


library("DESeq2")
dds_syn <- DESeqDataSetFromMatrix(countData = round(synovMat),
                                  colData = metadata_visit3_toci,
                                  design = ~ CDAI.response.status.V7
)


# Select top 500 most variable genes
exprMat_syn <- assay(dds_syn)
nTop <- 500
sds_syn <- genefilter::rowSds(exprMat_syn)
exprMat_syn <- exprMat_syn[order(sds_syn, decreasing = T)[1:nTop], ]
exprMat_syn[, c(1:58)] <- sapply(exprMat_syn[, c(1:58)], as.numeric)
str(exprMat_syn)
# Data distribution of metabolomics data
boxplot(exprMat_syn, outline = FALSE, col = "cornflowerblue", main = "Transformed RNAseq Synovium data")
exprMat_syn <- as.matrix(exprMat_syn)
str(exprMat_syn)





# Prepare Blood RNA-seq data --------------------------------------------

library("DESeq2")
dds_blood <- DESeqDataSetFromMatrix(countData = round(bloodMat),
                                    colData = metadata_visit3_toci,
                                    design = ~ CDAI.response.status.V7
)


# Select top 500 most variable genes
exprMat_blood <- assay(dds_blood)
nTop <- 500
sds_blood <- genefilter::rowSds(exprMat_blood)
exprMat_blood <- exprMat_blood[order(sds_blood, decreasing = T)[1:nTop], ]
exprMat_blood[, c(1:58)] <- sapply(exprMat_blood[, c(1:58)], as.numeric)
str(exprMat_blood)
# Data distribution of metabolomics data
boxplot(exprMat_blood, outline = FALSE, col = "cornflowerblue", main = "Transformed RNAseq Blood data")
exprMat_blood <- as.matrix(exprMat_blood)
str(exprMat_blood)




# Create the MOFA obejct --------------------------------------------------

# List of data matrix
#mofaData <- list(mRNA_synovium = exprMat_syn, mRNA_Blood = exprMat_blood, Lipidomics = lipidomicsMat, Proteins = proteoMat, Metabolites = metaMat)
mofaData <- list(mRNA_Blood = exprMat_blood, Lipidomics = lipidomicsMat, Proteins = proteoMat, Metabolites = metaMat)
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
metadata_visit3_toci <- subset(metadata, Visit == "3" & Randomized.medication == "Tocilizumab")
metadata_visit3_toci <- metadata_visit3_toci[metadata_visit3_toci$Patient.I.D. %in% useSamples, ]
metadata_visit3_toci <- metadata_visit3_toci[,-c(1,2)]
colnames(metadata_visit3_toci)[1] <- 'sample'
samples_metadata(MOFAobject) <- metadata_visit3_toci
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
            factor =5, 
            color_by = "Factor5"
)

plot_weights(MOFAobject,
             view = "Lipidomics",
             factors = 5,
             nfeatures = 20,     # Top number of features to highlight
             scale = F,           # Scale weights from -1 to 1
)




plot_top_weights(MOFAobject,
                 view = "Lipidomics",
                 factor =5,
                 nfeatures = 20,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "mRNA_Blood",
                 factor =3,
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
                  view = "Lipidomics",
                  factor =3,  
                  features = 20,
                  sign = "positive",
                  color_by = "CDAI.response.status.V7"
) + labs(y="Lipidomics")





# Positive CDAI
plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =5,  
                  features = 20,
                  sign = "positive",
                  color_by = "CDAI.response.status.V7"
) + labs(y="RNA expression Blood")


plot_data_scatter(MOFAobject, 
                  view = "mRNA_Blood",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
) + labs(y="RNA expression Blood")

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
)


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
                  sign = "positive",
                  color_by = "CDAI.response.status.V7"
)+ labs(y="Lipidomics Expression")

plot_data_scatter(MOFAobject, 
                  view = "Lipidomics",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)+ labs(y="Lipidomics Expression")

plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =5,  
                  features = 20,
                  sign = "positive",
                  color_by = "CDAI.response.status.V7"
)+ labs(y="Metabolomics Expression")

plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =5,  
                  features = 20,
                  sign = "negative",
                  color_by = "CDAI.response.status.V7"
)+ labs(y="Metabolomics Expression")

plot_data_scatter(MOFAobject, 
                  view = "Metabolites",
                  factor =3,  
                  features = 20,
                  sign = "positive",
                  color_by = "CDAI.response.status.V7"
) 


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

