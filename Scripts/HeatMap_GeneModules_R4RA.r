rm(list=ls())
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(dplyr)
library(xtable)
library(readxl)

#Load the gene modules matrix
mat <- read.delim("/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/Pitzalis/scripts/Gene_Module/R4RA/Master_Table_Parsed.tsv")
mat_sort <- mat[order(mat$Id),]
samples <- mat_sort$Id
dim(mat)

#Load the metadata
metadata <-read_xlsx("/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/Pitzalis/scripts/Gene_Module/R4RA/R4RA_full_metadata_20221111.xlsx")
metadata <- as.data.frame(metadata)
str(metadata)

#Quality check that metadata patient Id and gene module matrix patient Id are in agreement
R4RA <- metadata[metadata$QiagenID_Synovium %in% samples, ]
R4RA_sort <- R4RA[order(R4RA$QiagenID_Synovium),]
R4RA_sort$QiagenID_Synovium==mat_sort$Id

#subset the metadata by Visit (3 or 7) first and by Drug (Rituximab or Tocilizumab)
metadata_visit3_visit7 <- subset(metadata,Visit=="3" | Visit=="7") 
metadata_visit3_visit7_rituximab_tocilizumab <- subset(metadata_visit3_visit7, Randomized.medication=="Rituximab"| Randomized.medication=="Tocilizumab")
length(unique(metadata_visit3_visit7_rituximab_tocilizumab$QiagenID_Synovium))


#remove from the metadata all the samples that are not outliers (by selecting only samples that are in the gene module matrix)
metadata_visit3_visit7_rituximab_tocilizumab_no_outliers <- metadata_visit3_visit7_rituximab_tocilizumab[metadata_visit3_visit7_rituximab_tocilizumab$QiagenID_Synovium %in% samples, ]
metadata_visit3_visit7_rituximab_tocilizumab_no_outliers <- metadata_visit3_visit7_rituximab_tocilizumab_no_outliers[,-c(1)]
metadata_visit3_visit7_rituximab_tocilizumab_no_outliers <- as.data.frame(metadata_visit3_visit7_rituximab_tocilizumab_no_outliers)

#subset the metadata with the clinical and patient information to plot
metadata <- R4RA_sort[,c("QiagenID_Synovium","Visit","Randomized.medication","CDAI.response.status.V7","Pathotype","CD3","CD20","CD68L","CD68SL","CD138","DAS28.ESR","DAS28.CRP",
                         "VAS.Global","VAS.physician","RF.Visit1","CCP.Visit1","DAS28.ESR.status","DAS28.ESR.status.V7","DAS28.CRP.status","DAS28.CRP.status.V7","BMI","Tender.Joint.Count",
                         "Swollen.Joint.Count","Neutrophils","Age.Visit1")]

#obatain samples' metadata in different conditions
metadata_visit3 <- subset(metadata,Visit=="3")
metadata_visit3_rituximab <- subset(metadata_visit3,Randomized.medication=="Rituximab")
metadata_visit3_tocilizumab <- subset(metadata_visit3,Randomized.medication=="Tocilizumab")
metadata_visit7 <- subset(metadata,Visit=="7")
metadata_visit7_rituximab <- subset(metadata_visit7,Randomized.medication=="Rituximab")
metadata_visit7_tocilizumab <- subset(metadata_visit7,Randomized.medication=="Tocilizumab")
metadata_rituximab <- subset(metadata,Randomized.medication=="Rituximab")
metadata_tocilizumab <- subset(metadata,Randomized.medication=="Tocilizumab")

#exatract only the patient Id for each combination
samples_vist3 <- metadata_visit3$QiagenID_Synovium
samples_vist7 <- metadata_visit7$QiagenID_Synovium
samples_rituximab <- metadata_rituximab$QiagenID_Synovium
samples_tocilizumab <- metadata_tocilizumab$QiagenID_Synovium
samples_vist3_tocilizumab <- metadata_visit3_tocilizumab$QiagenID_Synovium
samples_vist7_tocilizumab <- metadata_visit7_tocilizumab$QiagenID_Synovium
samples_vist3_rituximab <- metadata_visit3_rituximab$QiagenID_Synovium
samples_vist7_rituximab <- metadata_visit7_rituximab$QiagenID_Synovium


#obtain module score for each treatment and group
R4RA_gene_modules <- data.frame(mat_sort, row.names = 1)
mat_sort

mat_sort_visit3 <- mat_sort[mat_sort$Id %in% samples_vist3,]
R4RA_gene_modules_visit3 <- data.frame(mat_sort_visit3, row.names = 1)

mat_sort_visit3_rituximab <- mat_sort[mat_sort$Id %in% samples_vist3_rituximab,]
R4RA_gene_modules_visit3_rituximab  <- data.frame(mat_sort_visit3_rituximab, row.names = 1)

mat_sort_visit3_tocilizumab <- mat_sort[mat_sort$Id %in% samples_vist3_tocilizumab,]
R4RA_gene_modules_visit3_tocilizumab  <- data.frame(mat_sort_visit3_tocilizumab, row.names = 1)

mat_sort_visit7 <- mat_sort[mat_sort$Id %in% samples_vist7,]
R4RA_gene_modules_visit7 <- data.frame(mat_sort_visit7, row.names = 1)

mat_sort_visit7_rituximab <- mat_sort[mat_sort$Id %in% samples_vist7_rituximab,]
R4RA_gene_modules_visit7_rituximab  <- data.frame(mat_sort_visit7_rituximab, row.names = 1)

mat_sort_visit7_tocilizumab <- mat_sort[mat_sort$Id %in% samples_vist7_tocilizumab,]
R4RA_gene_modules_visit7_tocilizumab  <- data.frame(mat_sort_visit7_tocilizumab, row.names = 1)

mat_sort_tocilizumab <- mat_sort[mat_sort$Id %in% samples_tocilizumab,]
R4RA_gene_modules_tocilizumab <- data.frame(mat_sort_tocilizumab, row.names = 1)

mat_sort_rituximab <- mat_sort[mat_sort$Id %in% samples_rituximab,]
R4RA_gene_modules_rituximab <- data.frame(mat_sort_rituximab, row.names = 1)



#obtain the matrices in the format for the plot (by transposing them)
mat <- t(R4RA_gene_modules)
mat_visit3 <- t(as.matrix(R4RA_gene_modules_visit3))
mat_visit3_rituximab <- t(as.matrix(R4RA_gene_modules_visit3_rituximab))
mat_visit3_tocilizumab <- t(as.matrix(R4RA_gene_modules_visit3_tocilizumab))
mat_visit7 <- t(as.matrix(R4RA_gene_modules_visit7))
mat_visit7_rituximab <- t(as.matrix(R4RA_gene_modules_visit7_rituximab))
mat_visit7_tocilizumab <- t(as.matrix(R4RA_gene_modules_visit7_tocilizumab))
mat_rituximab <- t(as.matrix(R4RA_gene_modules_rituximab))
mat_tocilizumab <- t(as.matrix(R4RA_gene_modules_tocilizumab))


#create an object for the variables to include in the plots (one for each type of plot)
ann <- data.frame(
  Pathotype = as.character(metadata$Pathotype),
  CD3 = as.numeric(metadata$CD3),
  CD20 = as.numeric(metadata$CD20),
  CD68L = as.numeric(metadata$CD68L),
  CD68SL = as.numeric(metadata$CD68SL),
  CD138 = as.numeric(metadata$CD138),
  Visit = as.character(metadata$Visit),
  Randomized.medication = as.character(metadata$Randomized.medication),
  CDAI.response.status.V7 =  as.character(metadata$CDAI.response.status.V7),
  VAS.Global = as.numeric(metadata$VAS.Global),
  VAS.physician = as.numeric(metadata$VAS.physician),
  RF.Visit1 = as.numeric(metadata$RF.Visit1),
  CCP.Visit1 = as.numeric(metadata$CCP.Visit1),
  DAS28.CRP = as.numeric(metadata$DAS28.CRP),
  DAS28.ESR = as.numeric(metadata$DAS28.ESR),
  DAS28.CRP.status = as.character(metadata$DAS28.CRP.status),
  DAS28.CRP.status.V7 = as.character(metadata$DAS28.CRP.status.V7),
  DAS28.ESR.status = as.character(metadata$DAS28.ESR.status),
  DAS28.ESR.status.V7 = as.character(metadata$DAS28.ESR.status.V7),
  BMI = as.numeric(metadata$BMI),
  Tender.Joint.Count = as.numeric(metadata$Tender.Joint.Count),
  Swollen.Joint.Count = as.numeric(metadata$Swollen.Joint.Count),
  Neutrophils = as.numeric(metadata$Neutrophils),
  Age.Visit1 = as.numeric(metadata$Age.Visit1),
  stringsAsFactors = FALSE)


#visit 3
ann_visit3 <- data.frame(
  Pathotype = as.character(metadata_visit3$Pathotype),
  CD3 = as.numeric(metadata_visit3$CD3),
  CD20 = as.numeric(metadata_visit3$CD20),
  CD68L = as.numeric(metadata_visit3$CD68L),
  CD68SL = as.numeric(metadata_visit3$CD68SL),
  CD138 = as.numeric(metadata_visit3$CD138),
  Visit = as.character(metadata_visit3$Visit),
  Randomized.medication = as.character(metadata_visit3$Randomized.medication),
  CDAI.response.status.V7 =  as.character(metadata_visit3$CDAI.response.status.V7),
  VAS.Global = as.numeric(metadata_visit3$VAS.Global),
  VAS.physician = as.numeric(metadata_visit3$VAS.physician),
  RF.Visit1 = as.numeric(metadata_visit3$RF.Visit1),
  CCP.Visit1 = as.numeric(metadata_visit3$CCP.Visit1),
  DAS28.CRP = as.numeric(metadata_visit3$DAS28.CRP),
  DAS28.ESR = as.numeric(metadata_visit3$DAS28.ESR),
  DAS28.CRP.status = as.character(metadata_visit3$DAS28.CRP.status),
  DAS28.CRP.status.V7 = as.character(metadata_visit3$DAS28.CRP.status.V7),
  DAS28.ESR.status = as.character(metadata_visit3$DAS28.ESR.status),
  DAS28.ESR.status.V7 = as.character(metadata_visit3$DAS28.ESR.status.V7),
  BMI = as.numeric(metadata_visit3$BMI),
  Tender.Joint.Count = as.numeric(metadata_visit3$Tender.Joint.Count),
  Swollen.Joint.Count = as.numeric(metadata_visit3$Swollen.Joint.Count),
  Neutrophils = as.numeric(metadata_visit3$Neutrophils),
  Age.Visit1 = as.numeric(metadata_visit3$Age.Visit1),
  stringsAsFactors = FALSE)

#visit 7
ann_visit7 <- data.frame(
  Pathotype = as.character(metadata_visit7$Pathotype),
  CD3 = as.numeric(metadata_visit7$CD3),
  CD20 = as.numeric(metadata_visit7$CD20),
  CD68L = as.numeric(metadata_visit7$CD68L),
  CD68SL = as.numeric(metadata_visit7$CD68SL),
  CD138 = as.numeric(metadata_visit7$CD138),
  Visit = as.character(metadata_visit7$Visit),
  Randomized.medication = as.character(metadata_visit7$Randomized.medication),
  CDAI.response.status.V7 =  as.character(metadata_visit7$CDAI.response.status.V7),
  VAS.Global = as.numeric(metadata_visit7$VAS.Global),
  VAS.physician = as.numeric(metadata_visit7$VAS.physician),
  RF.Visit1 = as.numeric(metadata_visit7$RF.Visit1),
  CCP.Visit1 = as.numeric(metadata_visit7$CCP.Visit1),
  DAS28.CRP = as.numeric(metadata_visit7$DAS28.CRP),
  DAS28.ESR = as.numeric(metadata_visit7$DAS28.ESR),
  DAS28.CRP.status = as.character(metadata_visit7$DAS28.CRP.status),
  DAS28.CRP.status.V7 = as.character(metadata_visit7$DAS28.CRP.status.V7),
  DAS28.ESR.status = as.character(metadata_visit7$DAS28.ESR.status),
  DAS28.ESR.status.V7 = as.character(metadata_visit7$DAS28.ESR.status.V7),
  BMI = as.numeric(metadata_visit7$BMI),
  Tender.Joint.Count = as.numeric(metadata_visit7$Tender.Joint.Count),
  Swollen.Joint.Count = as.numeric(metadata_visit7$Swollen.Joint.Count),
  Neutrophils = as.numeric(metadata_visit7$Neutrophils),
  Age.Visit1 = as.numeric(metadata_visit7$Age.Visit1),
  stringsAsFactors = FALSE)

#Rituximab
ann_rituximab <- data.frame(
  Pathotype = as.character(metadata_rituximab$Pathotype),
  CD3 = as.numeric(metadata_rituximab$CD3),
  CD20 = as.numeric(metadata_rituximab$CD20),
  CD68L = as.numeric(metadata_rituximab$CD68L),
  CD68SL = as.numeric(metadata_rituximab$CD68SL),
  CD138 = as.numeric(metadata_rituximab$CD138),
  Visit = as.character(metadata_rituximab$Visit),
  Randomized.medication = as.character(metadata_rituximab$Randomized.medication),
  CDAI.response.status.V7 =  as.character(metadata_rituximab$CDAI.response.status.V7),
  VAS.Global = as.numeric(metadata_rituximab$VAS.Global),
  VAS.physician = as.numeric(metadata_rituximab$VAS.physician),
  RF.Visit1 = as.numeric(metadata_rituximab$RF.Visit1),
  CCP.Visit1 = as.numeric(metadata_rituximab$CCP.Visit1),
  DAS28.CRP = as.numeric(metadata_rituximab$DAS28.CRP),
  DAS28.ESR = as.numeric(metadata_rituximab$DAS28.ESR),
  DAS28.CRP.status = as.character(metadata_rituximab$DAS28.CRP.status),
  DAS28.CRP.status.V7 = as.character(metadata_rituximab$DAS28.CRP.status.V7),
  DAS28.ESR.status = as.character(metadata_rituximab$DAS28.ESR.status),
  DAS28.ESR.status.V7 = as.character(metadata_rituximab$DAS28.ESR.status.V7),
  BMI = as.numeric(metadata_rituximab$BMI),
  Tender.Joint.Count = as.numeric(metadata_rituximab$Tender.Joint.Count),
  Swollen.Joint.Count = as.numeric(metadata_rituximab$Swollen.Joint.Count),
  Neutrophils = as.numeric(metadata_rituximab$Neutrophils),
  Age.Visit1 = as.numeric(metadata_rituximab$Age.Visit1),
  stringsAsFactors = FALSE)

#Rituximab visit 3
ann_rituximab_visit3 <- data.frame(
  Pathotype = as.character(metadata_visit3_rituximab$Pathotype),
  CD3 = as.numeric(metadata_visit3_rituximab$CD3),
  CD20 = as.numeric(metadata_visit3_rituximab$CD20),
  CD68L = as.numeric(metadata_visit3_rituximab$CD68L),
  CD68SL = as.numeric(metadata_visit3_rituximab$CD68SL),
  CD138 = as.numeric(metadata_visit3_rituximab$CD138),
  Visit = as.character(metadata_visit3_rituximab$Visit),
  Randomized.medication = as.character(metadata_visit3_rituximab$Randomized.medication),
  CDAI.response.status.V7 =  as.character(metadata_visit3_rituximab$CDAI.response.status.V7),
  VAS.Global = as.numeric(metadata_visit3_rituximab$VAS.Global),
  VAS.physician = as.numeric(metadata_visit3_rituximab$VAS.physician),
  RF.Visit1 = as.numeric(metadata_visit3_rituximab$RF.Visit1),
  CCP.Visit1 = as.numeric(metadata_visit3_rituximab$CCP.Visit1),
  DAS28.CRP = as.numeric(metadata_visit3_rituximab$DAS28.CRP),
  DAS28.ESR = as.numeric(metadata_visit3_rituximab$DAS28.ESR),
  DAS28.CRP.status = as.character(metadata_visit3_rituximab$DAS28.CRP.status),
  DAS28.CRP.status.V7 = as.character(metadata_visit3_rituximab$DAS28.CRP.status.V7),
  DAS28.ESR.status = as.character(metadata_visit3_rituximab$DAS28.ESR.status),
  DAS28.ESR.status.V7 = as.character(metadata_visit3_rituximab$DAS28.ESR.status.V7),
  BMI = as.numeric(metadata_visit3_rituximab$BMI),
  Tender.Joint.Count = as.numeric(metadata_visit3_rituximab$Tender.Joint.Count),
  Swollen.Joint.Count = as.numeric(metadata_visit3_rituximab$Swollen.Joint.Count),
  Neutrophils = as.numeric(metadata_visit3_rituximab$Neutrophils),
  Age.Visit1 = as.numeric(metadata_visit3_rituximab$Age.Visit1),
  stringsAsFactors = FALSE)


#Rituximab visit 7
ann_rituximab_visit7 <- data.frame(
  Pathotype = as.character(metadata_visit7_rituximab$Pathotype),
  CD3 = as.numeric(metadata_visit7_rituximab$CD3),
  CD20 = as.numeric(metadata_visit7_rituximab$CD20),
  CD68L = as.numeric(metadata_visit7_rituximab$CD68L),
  CD68SL = as.numeric(metadata_visit7_rituximab$CD68SL),
  CD138 = as.numeric(metadata_visit7_rituximab$CD138),
  Visit = as.character(metadata_visit7_rituximab$Visit),
  Randomized.medication = as.character(metadata_visit7_rituximab$Randomized.medication),
  CDAI.response.status.V7 =  as.character(metadata_visit7_rituximab$CDAI.response.status.V7),
  VAS.Global = as.numeric(metadata_visit7_rituximab$VAS.Global),
  VAS.physician = as.numeric(metadata_visit7_rituximab$VAS.physician),
  RF.Visit1 = as.numeric(metadata_visit7_rituximab$RF.Visit1),
  CCP.Visit1 = as.numeric(metadata_visit7_rituximab$CCP.Visit1),
  DAS28.CRP = as.numeric(metadata_visit7_rituximab$DAS28.CRP),
  DAS28.ESR = as.numeric(metadata_visit7_rituximab$DAS28.ESR),
  DAS28.CRP.status = as.character(metadata_visit7_rituximab$DAS28.CRP.status),
  DAS28.CRP.status.V7 = as.character(metadata_visit7_rituximab$DAS28.CRP.status.V7),
  DAS28.ESR.status = as.character(metadata_visit7_rituximab$DAS28.ESR.status),
  DAS28.ESR.status.V7 = as.character(metadata_visit7_rituximab$DAS28.ESR.status.V7),
  BMI = as.numeric(metadata_visit7_rituximab$BMI),
  Tender.Joint.Count = as.numeric(metadata_visit7_rituximab$Tender.Joint.Count),
  Swollen.Joint.Count = as.numeric(metadata_visit7_rituximab$Swollen.Joint.Count),
  Neutrophils = as.numeric(metadata_visit7_rituximab$Neutrophils),
  Age.Visit1 = as.numeric(metadata_visit7_rituximab$Age.Visit1),
  stringsAsFactors = FALSE)

#Tocilizumab
ann_tocilizumab <- data.frame(
  Pathotype = as.character(metadata_tocilizumab$Pathotype),
  CD3 = as.numeric(metadata_tocilizumab$CD3),
  CD20 = as.numeric(metadata_tocilizumab$CD20),
  CD68L = as.numeric(metadata_tocilizumab$CD68L),
  CD68SL = as.numeric(metadata_tocilizumab$CD68SL),
  CD138 = as.numeric(metadata_tocilizumab$CD138),
  Visit = as.character(metadata_tocilizumab$Visit),
  Randomized.medication = as.character(metadata_tocilizumab$Randomized.medication),
  CDAI.response.status.V7 =  as.character(metadata_tocilizumab$CDAI.response.status.V7),
  VAS.Global = as.numeric(metadata_tocilizumab$VAS.Global),
  VAS.physician = as.numeric(metadata_tocilizumab$VAS.physician),
  RF.Visit1 = as.numeric(metadata_tocilizumab$RF.Visit1),
  CCP.Visit1 = as.numeric(metadata_tocilizumab$CCP.Visit1),
  DAS28.CRP = as.numeric(metadata_tocilizumab$DAS28.CRP),
  DAS28.ESR = as.numeric(metadata_tocilizumab$DAS28.ESR),
  DAS28.CRP.status = as.character(metadata_tocilizumab$DAS28.CRP.status),
  DAS28.CRP.status.V7 = as.character(metadata_tocilizumab$DAS28.CRP.status.V7),
  DAS28.ESR.status = as.character(metadata_tocilizumab$DAS28.ESR.status),
  DAS28.ESR.status.V7 = as.character(metadata_tocilizumab$DAS28.ESR.status.V7),
  BMI = as.numeric(metadata_tocilizumab$BMI),
  Tender.Joint.Count = as.numeric(metadata_tocilizumab$Tender.Joint.Count),
  Swollen.Joint.Count = as.numeric(metadata_tocilizumab$Swollen.Joint.Count),
  Neutrophils = as.numeric(metadata_tocilizumab$Neutrophils),
  Age.Visit1 = as.numeric(metadata_tocilizumab$Age.Visit1),
  stringsAsFactors = FALSE)

#Tocilizumab visit3
ann_tocilizumab_visit3 <- data.frame(
  Pathotype = as.character(metadata_visit3_tocilizumab$Pathotype),
  CD3 = as.numeric(metadata_visit3_tocilizumab$CD3),
  CD20 = as.numeric(metadata_visit3_tocilizumab$CD20),
  CD68L = as.numeric(metadata_visit3_tocilizumab$CD68L),
  CD68SL = as.numeric(metadata_visit3_tocilizumab$CD68SL),
  CD138 = as.numeric(metadata_visit3_tocilizumab$CD138),
  Visit = as.character(metadata_visit3_tocilizumab$Visit),
  Randomized.medication = as.character(metadata_visit3_tocilizumab$Randomized.medication),
  CDAI.response.status.V7 =  as.character(metadata_visit3_tocilizumab$CDAI.response.status.V7),
  VAS.Global = as.numeric(metadata_visit3_tocilizumab$VAS.Global),
  VAS.physician = as.numeric(metadata_visit3_tocilizumab$VAS.physician),
  RF.Visit1 = as.numeric(metadata_visit3_tocilizumab$RF.Visit1),
  CCP.Visit1 = as.numeric(metadata_visit3_tocilizumab$CCP.Visit1),
  DAS28.CRP = as.numeric(metadata_visit3_tocilizumab$DAS28.CRP),
  DAS28.ESR = as.numeric(metadata_visit3_tocilizumab$DAS28.ESR),
  DAS28.CRP.status = as.character(metadata_visit3_tocilizumab$DAS28.CRP.status),
  DAS28.CRP.status.V7 = as.character(metadata_visit3_tocilizumab$DAS28.CRP.status.V7),
  DAS28.ESR.status = as.character(metadata_visit3_tocilizumab$DAS28.ESR.status),
  DAS28.ESR.status.V7 = as.character(metadata_visit3_tocilizumab$DAS28.ESR.status.V7),
  BMI = as.numeric(metadata_visit3_tocilizumab$BMI),
  Tender.Joint.Count = as.numeric(metadata_visit3_tocilizumab$Tender.Joint.Count),
  Swollen.Joint.Count = as.numeric(metadata_visit3_tocilizumab$Swollen.Joint.Count),
  Neutrophils = as.numeric(metadata_visit3_tocilizumab$Neutrophils),
  Age.Visit1 = as.numeric(metadata_visit3_tocilizumab$Age.Visit1),
  stringsAsFactors = FALSE)

#Tocilizumab visit7
ann_tocilizumab_visit7 <- data.frame(
  Pathotype = as.character(metadata_visit7_tocilizumab$Pathotype),
  CD3 = as.numeric(metadata_visit7_tocilizumab$CD3),
  CD20 = as.numeric(metadata_visit7_tocilizumab$CD20),
  CD68L = as.numeric(metadata_visit7_tocilizumab$CD68L),
  CD68SL = as.numeric(metadata_visit7_tocilizumab$CD68SL),
  CD138 = as.numeric(metadata_visit7_tocilizumab$CD138),
  Visit = as.character(metadata_visit7_tocilizumab$Visit),
  Randomized.medication = as.character(metadata_visit7_tocilizumab$Randomized.medication),
  CDAI.response.status.V7 =  as.character(metadata_visit7_tocilizumab$CDAI.response.status.V7),
  VAS.Global = as.numeric(metadata_visit7_tocilizumab$VAS.Global),
  VAS.physician = as.numeric(metadata_visit7_tocilizumab$VAS.physician),
  RF.Visit1 = as.numeric(metadata_visit7_tocilizumab$RF.Visit1),
  CCP.Visit1 = as.numeric(metadata_visit7_tocilizumab$CCP.Visit1),
  DAS28.CRP = as.numeric(metadata_visit7_tocilizumab$DAS28.CRP),
  DAS28.ESR = as.numeric(metadata_visit7_tocilizumab$DAS28.ESR),
  DAS28.CRP.status = as.character(metadata_visit7_tocilizumab$DAS28.CRP.status),
  DAS28.CRP.status.V7 = as.character(metadata_visit7_tocilizumab$DAS28.CRP.status.V7),
  DAS28.ESR.status = as.character(metadata_visit7_tocilizumab$DAS28.ESR.status),
  DAS28.ESR.status.V7 = as.character(metadata_visit7_tocilizumab$DAS28.ESR.status.V7),
  BMI = as.numeric(metadata_visit7_tocilizumab$BMI),
  Tender.Joint.Count = as.numeric(metadata_visit7_tocilizumab$Tender.Joint.Count),
  Swollen.Joint.Count = as.numeric(metadata_visit7_tocilizumab$Swollen.Joint.Count),
  Neutrophils = as.numeric(metadata_visit7_tocilizumab$Neutrophils),
  Age.Visit1 = as.numeric(metadata_visit7_tocilizumab$Age.Visit1),
  stringsAsFactors = FALSE)

#colors
Age.Visit1_colour = colorRamp2(c(0,81), c("white", "navy"))
Neutrophils_colour = colorRamp2(c(1,14), c("white", "green"))
SwollenJointCount_colour = colorRamp2(c(0,26), c("white", "pink"))
TenderJointCount_colour = colorRamp2(c(0,28), c("white", "purple"))
BMI_colour = colorRamp2(c(13,52), c("white","powderblue"))
DAS28_ESR_colour = colorRamp2(c(2,8), c("white", "cyan3"))
DAS28_CRP_colour = colorRamp2(c(2,9), c("white", "cyan3"))
CCP.Visit1_colour = colorRamp2(c(0,640), c("white", "coral"))
RF.Visit1_colour = colorRamp2(c(0,4405), c("white", "red"))
VAS.physician_colour = colorRamp2(c(2,100), c("white", "purple"))
VAS.Global_colour = colorRamp2(c(1,100), c("white", "purple"))


# create the colour mapping
colours <- list(
  Pathotype = c('Lymphoid' = 'blue', 'Myeloid' = 'red', 'Fibroid' = 'green3', 'TBC' = 'grey', 'Ungraded' = 'grey'),
  CD3 = c('0' = '#F7FCF5', '1' = '#C7E9C0', '2' = '#74C476', '3' = '#238B45', '4' = '#00441B'),
  CD20 = c('0' = '#F7FBFF', '1' = '#C6DBEF', '2' = '#6BAED6', '3' = '#2171B5', '4' = '#08306B'),
  CD68L = c('0' = '#FFF5F0', '1' = '#FCBBA1', '2' = '#FB6A4A', '3' = '#CB181D', '4' = '#67000D'),
  CD68SL = c('0' = '#FFF5EB', '1' = '#FDD0A2', '2' = '#FD8D3C', '3' = '#D94801', '4' = '#7F2704'),
  CD138 = c('0' = '#FFF5EB', '1' = '#FDD0A2', '2' = '#FD8D3C', '3' = '#D94801', '4' = '#7F2704','NA' = 'grey'),
  Visit = c('3'= 'red','7' = 'green'),
  Randomized.medication=c('Rituximab'='blue','Tocilizumab'='orange'),
  CDAI.response.status.V7=c('Responder' = 'turquoise4','Non.Responder' = 'black', 'NA' = 'grey'),
  VAS.Global = VAS.Global_colour,
  VAS.physician = VAS.physician_colour,
  RF.Visit1 = RF.Visit1_colour,
  CCP.Visit1 = CCP.Visit1_colour,
  DAS28.CRP = DAS28_CRP_colour,
  DAS28.CRP.status = c('Remission' = '#b0e2ff','Moderate.DA' = '#9ecbe5', 'Low.DA' = '#8cb4cc', 'High.DA' = 'black'),
  DAS28.CRP.status.V7 = c('Remission' = '#b0e2ff','Moderate.DA' = '#9ecbe5', 'Low.DA' = '#8cb4cc', 'High.DA' = 'black'),
  DAS28.ESR = DAS28_ESR_colour,
  DAS28.ESR.status = c('Remission' = '#b0e2ff','Moderate.DA' = '#9ecbe5', 'Low.DA' = '#8cb4cc', 'High.DA' = 'black'),
  DAS28.ESR.status.V7 = c('Remission' = '#b0e2ff','Moderate.DA' = '#9ecbe5', 'Low.DA' = '#8cb4cc', 'High.DA' = 'black'),
  BMI = BMI_colour,
  Tender.Joint.Count = TenderJointCount_colour,
  Swollen.Joint.Count = SwollenJointCount_colour,
  Neutrophils=Neutrophils_colour,
  Age.Visit1 = Age.Visit1_colour
)




colAnn_tocilizumab_visit7 <- HeatmapAnnotation(
  df = ann_tocilizumab_visit7,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'grey', # default colour for any NA values in the annotation data-frame, 'ann'
  col = colours,
  annotation_name_side = 'left',
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Pathotype = list(
      nrow = 4, # number of rows across which the legend will be arranged
      by_row = FALSE,
      title = 'Pathotype',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    CD3 = list(
      nrow = 5,
      by_row = FALSE,
      title = 'CD3',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    CD20 = list(
      nrow = 5,
      by_row = FALSE,
      title = 'CD20',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    CD68L = list(
      nrow = 5,
      by_row = FALSE,
      title = 'CD68L',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    CD68SL = list(
      nrow = 5,
      by_row = FALSE,
      title = 'CD68SL',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    CD138 = list(
      nrow = 5,
      by_row = FALSE,
      title = 'CD138',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    Visit = list(
      nrow = 4, # number of rows across which the legend will be arranged
      by_row = FALSE,
      title = 'Visit',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    Randomized.medication = list(
      nrow = 4, # number of rows across which the legend will be arranged
      by_row = FALSE,
      title = 'Randomized.medication',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    CDAI.response.status.V7= list(
      nrow = 4, # number of rows across which the legend will be arranged
      by_row = FALSE,
      title = 'CDAI.response.status.V7',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    VAS.Global =list(
      nrow = 4, # number of rows across which the legend will be arranged
      by_row = FALSE,
      title = 'VAS.Global',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    VAS.physician = list(
      nrow = 4, # number of rows across which the legend will be arranged
      by_row = FALSE,
      title = 'VAS.physician',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    RF.Visit1 = list(
      nrow = 4, # number of rows across which the legend will be arranged
      by_row = FALSE,
      title = 'RF.Visit1',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    CCP.Visit1 = list(
      nrow = 4, # number of rows across which the legend will be arranged
      by_row = FALSE,
      title = 'CCP.Visit1',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    DAS28.CRP = list(
      nrow = 5,
      by_row = FALSE,
      title = 'DAS28.CRP',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    DAS28.CRP.status = list(
      nrow = 5,
      by_row = FALSE,
      title = 'DAS28.CRP.status',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    DAS28.CRP.status.V7 = list(
      nrow = 5,
      by_row = FALSE,
      title = 'DAS28.CRP.status.V7',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    DAS28.ESR = list(
      nrow = 5,
      by_row = FALSE,
      title = 'DAS28.ESR',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    DAS28.ESR.status = list(
      nrow = 5,
      by_row = FALSE,
      title = 'DAS28.ESR.status',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    DAS28.ESR.status.V7 = list(
      nrow = 5,
      by_row = FALSE,
      title = ' DAS28.ESR.status.V7',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    BMI = list(
      nrow = 5,
      by_row = FALSE,
      title = 'BMI',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    Tender.Joint.Count = list(
      nrow = 5,
      by_row = FALSE,
      title = 'Tender.Joint.Count',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    Swollen.Joint.Count = list(
      nrow = 5,
      by_row = FALSE,
      title = 'Swollen.Joint.Count',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    Neutrophils = list(
      nrow = 5,
      by_row = FALSE,
      title = 'Neutrophils',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold')),
    Age.Visit1 = list(
      nrow = 5,
      by_row = FALSE,
      title = 'Age.Visit1',
      title_position = 'lefttop-rot',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 8, fontface = 'bold'),
      labels_gp = gpar(fontsize = 8, fontface = 'bold'))))

heat_visit3 <- mat_visit3
heat_visit7 <- mat_visit7
heat_rituximab <- mat_rituximab
heat_rituximab_visit3 <- mat_visit3_rituximab
heat_rituximab_visit7 <- mat_visit7_rituximab
heat_tocilizumab <- mat_tocilizumab
heat_tocilizumab_visit3 <- mat_visit3_tocilizumab
heat_tocilizumab_visit7 <- mat_visit7_tocilizumab

hmap <- Heatmap(heat_tocilizumab_visit7,
                top_annotation = colAnn_tocilizumab_visit7,
                show_row_names = TRUE,
                name = 'Module Score',
                row_dend_reorder = FALSE,
                column_dend_reorder = FALSE,
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 8, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),
                )
draw(hmap)

