# Genome-wide cline analysis identifies new locus contributing to a barrier to gene flow across an _Antirrhinum_ hybrid zone

This repository documents the analyses scripts associated with the manuscript "Genome-wide cline analysis identifies new locus contributing to a barrier to gene flow across an _Antirrhinum_ hybrid zone" by David L Field, Sean Stankowski, Taylor Reiter, Jitka Polechova, Desmond Bradley, Daniel Richardson, Annabel Whibley, Arka Pal, Daria Shipilina, Louis Boell, Melinda Pickup, Yongbiao Xue, Enrico Coen, Nicholas Barton. https://www.biorxiv.org/content/10.1101/2025.02.17.638607v1.full


## Table of Contents
- [Input Data sets](#input-data-sets)
- [Workflow](#workflow)
  - [(i) Data, functions, parameters](#i-data-functions-parameters)
  - [(ii) Genome scan sliding window](#ii-genome-scan-sliding-window)
  - [(iii) FastClines geographic](#iii-fastclines-geographic)
  - [(iv) Post processing sliding window and FastClines outputs](#iv-post-processing-sliding-window-and-fastclines-outputs)
  - [(v) RNAseq Differential Gene Expression (DGE)](#v-rnaseq-differential-gene-expression-dge)
  - [(vi) Genotype-Phenotype associations](#vi-genotype-phenotype-associations)
  - [(vii) MLE Cline fitting](#vii-mle-cline-fitting)
  - [(viii) FastClines simulations](#viii-fastclines-simulations)

## Scripts

Analyses and figures for main text and supplementary materials are conducted in a series of R markdown scripts contained in main folder of _https://github.com/dfield007/genome_wide_clines_

- FastClines_Amajus_main.Rmd
- FastClines_Amajus_RNAseq.Rmd
- FastClines_Amajus_GenoPheno.Rmd
- FastClines_Amajus_MLEclines.Rmd
- FastClines_Amajus_simulations.Rmd

and also a number of Python scripts called within the markdowns 

- SlidingWindows.py
- FastClines.py

further custom script functions called from R markdowns contained in `genome_wide_clines/myRfunctions`.

Most of the main analyses and generation of figures are conducted in the script `FastClines_Amajus_main.Rmd`

## Workflow pipeline summary

### (i) Data, functions, parameters

Input datasets:

- Antirrhinum Reference genome (online link)
- Antirrhinum Annotated genome: GWHBJVT00000000.gff (GFF annotations have been uploaded to NCBI WGS under accession number SUB15081867)
- Whole Genome PoolSeq dataset: (SRA under accession number SUB15081578)
- Whole Genome PoolSeq dataset population meta details: `Amajus_PoolSeq_popDetails.txt` (online link)
- RNAseq dataset: (SRA under accession number SUB15081578)
- RNAseq candidate genes for Chromosome 5: `Amajus_candidate_genes_Chr5_RNAseq_Rimport.csv` (Dryad link)
- Hybrid zone plant genotypes and spatial dataset: `` (Dryad link)
- Hybrid zone plant colour phenotype dataset: `Amajus_FlowerPhenoMeasures.csv` (Dryad link)

Functions: see `Source-my-functions` code chunk in `FastClines_Amajus_main.Rmd` to call all functions from `genome_wide_clines/myRfunctions`

### (ii) Genome scan sliding window

- Genome scan sliding windows (Fst, Dxy, pi) performed in `FastClines_Amajus_main.Rmd` 
- splits original sync files into separate chromosome chunks using BASH scripts
- Calls python script `SlidingWindows_v1.4.py` available here: _https://github.com/dfield007/slidingWindows_

### (iii) FastClines geographic

- Geographic cline analyses (cline centre and width) performed in `FastClines_Amajus_main.Rmd` 
- Calls python script `FastClines.py`available here: _https://github.com/dfield007/fastClines_

### (iv) post processing sliding window and FastClines outputs

- Performed by `FastClines_Amajus_main.Rmd` (unless stated otherwise)
- imports individual SlidingWindow outputs and concatenates chromosome outputs
- combine fastClines with SlidngWindows for permutation tests (_python code_)
- combining flower colour gene locations with sliding window and clines data
- Fig 1. map of hybrid zone sampling
- Fig 3a. genome wide cline width by cline centre
- Fig 4. summary box plots of cline width in different genomic regions
- Fig 5. Genome scans for cline properties, frequency and diversity and differentiation sliding window analyses
- Fig 6. 
- Fig 8a,b. Genome scans of Chromosome 5 and RUBIA region (Fst, Dxy, pi, cline width)
- Fig S5 - Fig S13

### (v) RNAseq Differential Gene Expression (DGE)

- main bioinformatics pipeline performed by `FastClines_Amajus_RNAseq.Rmd`
- post analyses and figures in`FastClines_Amajus_main.Rmd`
- Fig 8c. Differential gene expression across RUBIA region

### (vi) Genotype-Phenotype associations

- main analyses run by `FastClines_Amajus_GenoPheno.Rmd`
- Fig 8d. ROS/RUB interaction
- Fig S15 - Fig S17, Table S5

### (vii) MLE Cline fitting

- main analyses run by `FastClines_Amajus_MLEclines.Rmd`
- Fig 8d. ROS/RUB interaction
- Fig S3, Fig S4

### (viii) FastClines simulations

- main analyses run by `FastClines_Amajus_MLEclines.Rmd`
- creates figures for supplementary Fig S1 and Fig S2
