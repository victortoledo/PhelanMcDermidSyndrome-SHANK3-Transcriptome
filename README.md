# PhelanMcDermidSyndrome-SHANK3-Transcriptome
This repository contains the scripts used to analyze the transcriptomic profile of iPSC-derived neurons harboring mutations in the SHANK3 gene, associated with the Phelan-McDermid Syndrome

Codes are ordered according to the steps of the pipeline:

0trimmomatic.sh -> FASTQ files were trimmed using Trimmomatic

1alignSTAR.sh -> alignment to the human transcriptome was performed with STAR

2quantifyRSEM -> aligned reads were quantified by RSEM

3samplesPreprocessing.Rmd -> all samples were analyzed in order to identify clustering patterns, outliers and potential contamination

4diffExpression.Rmd -> counts were summarized to genes and differential expression analysis was performed

5wgcna.R -> co-expression patterns were identified using weighted gene co-expression network analysis (WGCNA)



If you have any question, feel free to e-mail me at victor.toledo.btos@gmail.com
