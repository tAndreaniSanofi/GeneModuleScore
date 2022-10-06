### GMS
Calculation of Gene Module Scores from gene expression matrices to characterize the cell type(s) in different patients.

### Input (check the folder input):  
The script/function to run  
Gene Expression Matrix (bring the one that you have from your samples)  
Gene annotation  
Cell type of interest (from a list of selected)  


### Output:  
A file where rows are patients and column is the cell type. Each value represents the gene module score.


### Description  
Example how to run the code on the bash command line:  
```
Rscript Gene_Module.R cell_types_and_genes.txt "Plasmacell" Human.B38_OmicsoftGenCode.V33.Genes.annotation3.txt Gene_Expression_Matrix.txt Plasmacell_Gene_module.txt  
```
## Exaplanation of the 4 input and 1 output  
Rscript is used to run an R script from unix command line while "Gene_Module.R" is the script to run  
The 5 main parameters that the script takes are:  
1) cell_types_and_genes.txt: file in which for each cell type there is a list of gene specific;  
2) "Plasmacell": Name of the cell to compute the gene module score, this is a full list (please use the exact word(s) using "":  
	-  Basophils
	-  CD14+ Monocytes
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
4) Gene_Expression_Matrix.txt: Matrix with the gene expression values  (row genes, columns Patient Id)
5) Plasmacell_Gene_module.txt: Output with the gene module score for the cell type under investigation for each patient. In the example "Plasmacell" but can be changed according to the cell type in the parameter 2  

## Metabolic Pathways (list from here: https://www.kegg.jp/kegg/pathway.html#metabolism). Only Human related pathway are used for the Module score.    
Name of the metabolic pathway to compute the module score, this is a full list (please use the exact word(s) using "":  
	- Carbon fixation in photosynthetic organisms  
	- Carbon fixation pathways in prokaryotes  
	- Methane metabolism  
	- Nitrogen metabolism  
	- Sulfur metabolism  
	- Fatty acid biosynthesis  
	- Fatty acid elongation  
	- Fatty acid degradation  
	- Cutin, suberine and wax biosynthesis  
	- Steroid biosynthesis  
	- Primary bile acid biosynthesis  
	- Secondary bile acid biosynthesis  
	- Steroid hormone biosynthesis  
	- Glycerolipid metabolism  
	- Glycerophospholipid metabolism  
	- Ether lipid metabolism  
	- Sphingolipid metabolism  
	- Arachidonic acid metabolism  
	- Linoleic acid metabolism  
	- alpha-Linolenic acid metabolism  
	- Biosynthesis of unsaturated fatty acids  
	- Purine metabolism  
	- Pyrimidine metabolism  
	- Alanine, aspartate and glutamate metabolism  
	- Glycine, serine and threonine metabolism  
	- Cysteine and methionine metabolism  
	- Valine, leucine and isoleucine degradation  
	- Valine, leucine and isoleucine biosynthesis  
	- Lysine biosynthesis  
	- Lysine degradation  
	- Arginine biosynthesis  
	- Arginine and proline metabolism  
	- Histidine metabolism  
	- Tyrosine metabolism  
	- Phenylalanine metabolism  
	- Tryptophan metabolism  
	- Phenylalanine, tyrosine and tryptophan biosynthesis  
	- beta-Alanine metabolism  
	- Taurine and hypotaurine metabolism  
	- Phosphonate and phosphinate metabolism  
	- Selenocompound metabolism  
	- Cyanoamino acid metabolism  
	- D-Amino acid metabolism (here)  
	- Glutathione metabolism  
	- N-Glycan biosynthesis  
	- Mannose type O-glycan biosynthesis  
	- Various types of N-glycan biosynthesis  
	- Mucin type O-glycan biosynthesis  
	- Mannose type O-glycan biosynthesis  
	- Other types of O-glycan biosynthesis  
	- Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate  
	- Glycosaminoglycan biosynthesis - heparan sulfate / heparin  
	- Glycosaminoglycan biosynthesis - keratan sulfate  
	- Glycosaminoglycan degradation  
	- Glycosylphosphatidylinositol (GPI)-anchor biosynthesis  
	- Glycosphingolipid biosynthesis - lacto and neolacto series  
	- Glycosphingolipid biosynthesis - globo and isoglobo series  
	- Glycosphingolipid biosynthesis - ganglio series  
	- Other glycan degradation  
	- Lipopolysaccharide biosynthesis  
	- O-Antigen repeat unit biosynthesis  
	- O-Antigen nucleotide sugar biosynthesis  
	- Peptidoglycan biosynthesis  
	- Teichoic acid biosynthesis  
	- Lipoarabinomannan (LAM) biosynthesis  
	- Arabinogalactan biosynthesis - Mycobacterium  
	- Exopolysaccharide biosynthesis New!  
	- Thiamine metabolism  
	- Riboflavin metabolism  
	- Vitamin B6 metabolism  
	- Nicotinate and nicotinamide metabolism  
	- Pantothenate and CoA biosynthesis  
	- Biotin metabolism  
	- Lipoic acid metabolism  
	- Folate biosynthesis  
	- One carbon pool by folate  
	- Retinol metabolism  
	- Porphyrin metabolism  
	- Ubiquinone and other terpenoid-quinone biosynthesis  
	- Terpenoid backbone biosynthesis  
	- Monoterpenoid biosynthesis  
	- Sesquiterpenoid and triterpenoid biosynthesis  
	- Diterpenoid biosynthesis  
	- Carotenoid biosynthesis  
	- Brassinosteroid biosynthesis  
	- Insect hormone biosynthesis  
	- Zeatin biosynthesis  
	- Limonene and pinene degradation  
	- Geraniol degradation  
	- Type I polyketide structures  
	- Biosynthesis of 12-, 14- and 16-membered macrolides  
	- Biosynthesis of ansamycins  
	- Biosynthesis of enediyne antibiotics  
	- Biosynthesis of type II polyketide backbone  
	- Biosynthesis of type II polyketide products  
	- Tetracycline biosynthesis  
	- Polyketide sugar unit biosynthesis  
	- Nonribosomal peptide structures  
	- Biosynthesis of siderophore group nonribosomal peptides  
	- Biosynthesis of vancomycin group antibiotics  
	- Phenylpropanoid biosynthesis  
	- Stilbenoid, diarylheptanoid and gingerol biosynthesis  
	- Flavonoid biosynthesis  
	- Flavone and flavonol biosynthesis  
	- Anthocyanin biosynthesis  
	- Isoflavonoid biosynthesis  
	- Indole alkaloid biosynthesis  
	- Indole diterpene alkaloid biosynthesis  
	- Isoquinoline alkaloid biosynthesis  
	- Tropane, piperidine and pyridine alkaloid biosynthesis  
	- Biosynthesis of various alkaloids  
	- Caffeine metabolism  
	- Betalain biosynthesis  
	- Glucosinolate biosynthesis  
	- Benzoxazinoid biosynthesis  
	- Penicillin and cephalosporin biosynthesis  
	- Carbapenem biosynthesis  
	- Monobactam biosynthesis  
	- Clavulanic acid biosynthesis  
	- Streptomycin biosynthesis  
	- Neomycin, kanamycin and gentamicin biosynthesis  
	- Acarbose and validamycin biosynthesis  
	- Novobiocin biosynthesis  
	- Staurosporine biosynthesis  
	- Phenazine biosynthesis  
	- Prodigiosin biosynthesis  
	- Aflatoxin biosynthesis  
	- Biosynthesis of various antibiotics  
	- Biosynthesis of various plant secondary metabolites  
	- Biosynthesis of various other secondary metabolites  
	- Benzoate degradation  
	- Aminobenzoate degradation  
	- Fluorobenzoate degradation  
	- Chloroalkane and chloroalkene degradation  
	- Chlorocyclohexane and chlorobenzene degradation  
	- Toluene degradation  
	- Xylene degradation  
	- Nitrotoluene degradation  
	- Ethylbenzene degradation  
	- Styrene degradation  
	- Atrazine degradation  
	- Caprolactam degradation  
	- Bisphenol degradation  
	- Dioxin degradation  
	- Naphthalene degradation  
	- Polycyclic aromatic hydrocarbon degradation  
	- Furfural degradation  
	- Steroid degradation  
	- Metabolism of xenobiotics by cytochrome P450  
	- Drug metabolism - cytochrome P450  
	- Drug metabolism - other enzymes  
	- Overview of biosynthetic pathways  
	- Biosynthesis of plant secondary metabolites  
	- Biosynthesis of phenylpropanoids  
	- Biosynthesis of terpenoids and steroids  
	- Biosynthesis of alkaloids derived from shikimate pathway  
	- Biosynthesis of alkaloids derived from ornithine, lysine and nicotinic acid  
	- Biosynthesis of alkaloids derived from histidine and purine  
