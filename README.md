# diet_rnaseq

Gene expression analyses 
This repo documents work done to analyze our gene expression analysis of cactus mice fed either a standard diet or a low fat diet. We parameterized aligning and counting RNA seq reads using STAR and htseq-count. For details on this, please see the ReadCount directory. We used a number of scripts we wrote for this purpose in the parameterization process, notably scripts to index the genome, align reads, and count expression as well as one that just aligns, and counts expression. We have additionally deposited scripts to extract mapping rates and other statistics for each sample and scripts to automate summary statistics on these mapping statistics within this directory.

Following this we conducted differential expression analyses using the R package DESeq2. We have a detailed R markdown files of our differential expression analyses. In addition to doing analyses, this script produces a number of the figures in the publication. Relevant data for this script (e.g., gene counts for each sample, genome annotation documents used) are all in the data directory.

Additionally, we used the WGCNA package in R to examine how genes are coexpressed and their correlation to treatment groups.

The scripts are split by tissue, analysis, and one big one. We suggest you open the diet_final.Rproj to open all the required files. To run the things we recommend starting with the rnaseq_dehy.Rmd file. We take advantage of Rs use of global variables in these analyses and this script WILL fail due to missing variables. This is when you will switch to running each tissue_GO.Rmd file and each WGCNA_tissue.Rmd file. Once you have done this head back to the main rnaseq_dehy.Rmd file where the code failed due to missing variables and start up again to finish it off. Is this good coding practice, nope and I learned while writing these scripts that this was a bad idea.
