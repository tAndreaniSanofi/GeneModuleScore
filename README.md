# GMS
Calculation of Gene Module Scores from gene expression matrices to characterize the cell type(s) in different patients.

### Input:  
The script/function to run  
Gene Expression Matrix  
Gene annotation  
Cell type of interest (from a list of selected)  


### Output:  
A file where rows are patients and column is the cell type. Each value represents the gene module score.


### Description  
Example how to run he code:  
```
Rscript Gene_Module.R cell_types_and_genes.txt "Plasmacell" Human.B38_OmicsoftGenCode.V33.Genes.annotation3.txt   RA_Synovium.PEAC.Baseline.RNASeq.Genes.Outliers_removed_Batch_corrected.Count.vst_20220503.txt Plasmacell_Gene_module.txt  
```
## Exaplanation of the 4 input  
Rscript is used to run an R script from unix command line while "Gene_Module.R" is the script to run  
The 5 main parameters that the script takes are:  
1) cell_types_and_genes.txt: file in which for each cell type there is a list of gene specific;  
2) "Plasmacell": Name of the cell to compute the gene module score, this is a full list (please use the exact word(s) using "":  
	- Basophils
	- CD14+ Monocytes
	- CD14+CD16+ Monocytes
	- CD14+CD16- Monocytes
	- CD14-CD16+ Monocytes
	- CD19+ B Cells
	- CD34+ Progenitors
	- CD4+ T Cells
	- CD4+CD25+CD45RA+ naive regulatory T cells
	- CD4+CD25+CD45RA- memory regulatory T cells
	- CD4+CD25-CD45RA+ naive conventional T cells
	- CD4+CD25-CD45RA- memory conventional T cells
	- CD8+ T Cells
	- Dendritic Cells - plasmacytoid
	- Endothelial Cells - Lymphatic
	- Endothelial Cells - Microvascular
	- Eosinophils
	- Fibroblast - skin
	- Mast cell
	- Natural Killer Cells
	- Neutrophils
	- Plasmacell
	- Synoviocyte
	- gamma delta positive T cells

3) Human.B38_OmicsoftGenCode.V33.Genes.annotation3.txt: annotation of the gene's name  
4) RA_Synovium.PEAC.Baseline.RNASeq.Genes.Outliers_removed_Batch_corrected.Count.vst_20220503.txt: Matrix with the gene expression values  
5) Plasmacell_Gene_module.txt: Output with the gene module score for the cell type under investigation for each patient. In the example "Plasmacell" but can be changed according to the cell type in the parameter 2  





