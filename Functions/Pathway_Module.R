##
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  all_pathways_and_genes <- args[1]
  pathway_of_interest <- args[2]
  annotation <- args[3]
  matrix <- args[4]
  output <- args[5]
 
 #Import tables with genes x cell types
 master_table <- read.table(file=all_pathways_and_genes,header = T,sep="\t")
 master_table <- master_table [, c(2,3)]

 #Subset the genes for each cell type
 pathway <- master_table[which(master_table$Pathway==pathway_of_interest),]
 print(pathway)
 #Subset Ensembl gene using the entrez genes of Cell Reports Paper
 Annotation <- read.csv(file=annotation,header=T,sep="\t")
 head(Annotation)
 Annotation <- Annotation [, c(1,4)]
 Annotation$EntrezID <- as.character(Annotation$EntrezID)
 s <- strsplit(Annotation$EntrezID, split = ";")
 Annotation_adapted <- data.frame(GeneID = rep(Annotation$GeneID, sapply(s, length)),  EntrezID= unlist(s))
 pathway_ensembl <- Annotation_adapted[Annotation_adapted$EntrezID %in% pathway$Gene, ]

 #Import the expression Matrix for RA_PEAC2_PEAC2 and and compute the gene module using the RLE (=rlog) values as per the manuscript
 RA_PEAC2_PEAC2_table <- read.table(file=matrix,header = T,sep="\t")
 row.names(RA_PEAC2_PEAC2_table) <- RA_PEAC2_PEAC2_table$ID
 RA_PEAC2_PEAC2_table[1] <- NULL
 # B.mod is list of the genes in your module
 pathway_ensembl_mod <- pathway_ensembl$Gene
 # gene exp. scaling
 # Gene_exp is expression matrix. Cols are samples, rows are genes.
 Gene_scaled <- scale(t(RA_PEAC2_PEAC2_table))
 Gene_scaled <- t(Gene_scaled)[colnames(Gene_scaled) %in% pathway_ensembl_mod, ]
 Gene_scaled <- as.data.frame(Gene_scaled)

 # calculate module scores
 Gene_module_scores <- colMeans(Gene_scaled,na.rm = T)
 Gene_module_scores_vst <- as.data.frame(Gene_module_scores)
 row.names(Gene_module_scores_vst)
 Id <- row.names(Gene_module_scores_vst)
 Gene_module_scores_vst <- as.data.frame(cbind(Id,Gene_module_scores_vst$Gene_module_scores))
 colnames(Gene_module_scores_vst)[2] <- pathway_of_interest 
 write.table(Gene_module_scores_vst,file=output,quote=F,row.names = F,col.names = T,sep="\t")
}
main()
