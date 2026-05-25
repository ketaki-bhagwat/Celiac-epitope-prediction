# Celiac-epitope-prediction
Predicts peptides from gliadin that bind to HLA alleles HLA-DQ2 and/or HLA-DQ8

Celiac disease is a chronic auto-immune disorder that manifests through triggering a strong immune response against gluten, a protein storage complex found in wheat. Gluten is composed of two polypeptides, gliadin and glutenin. Immune response against gluten is orchestrated by T-cells which recognise peptides from gliadin when presented on HLA class II alleles HLA-DQ8 and/or HLA-DQ2. Predicting gliadin peptides that bind to DQ2 or DQ8 can be useful to know which amino acids are involved in the peptide-MHC interactions and if their modifications can produce a wheat strain that can be tolerated by celiac disease patients. 
	This code breaks up a protein of interest (in this case, gliadin) into 15-mer peptides. It then sends the peptides in batches to IEDB's MHCII prediction tool for predicting their binding to HLA-DQ2 or HLADQ8. It returns a result that records the core peptide from every peptide sequence that binds to DQ2 and or DQ8, along with the score and rank which provide an insight into the strength of binding. It also identifies strong binders as those with a rank less than 10 and provides their list.

	Most of the peptides identified as strong binders by this code bind to HLA-DQ2, which is consistent with published literature on celiac disease. For your protein of interest, replace the uniprot ID with your protein of interest while calling the protein_fetcher function.

Dependencies - pandas,requests,io,IEDB API
