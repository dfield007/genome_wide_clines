# Genome-wide cline analysis identifies new locus contributing to a barrier to gene flow across an _Antirrhinum_ hybrid zone

This repository documents the analyses scripts associated with the manuscript "Genome-wide cline analysis identifies new locus contributing to a barrier to gene flow across an _Antirrhinum_ hybrid zone" by David L Field, Sean Stankowski, Taylor Reiter, Jitka Polechova, Desmond Bradley, Daniel Richardson, Arka Pal, Daria Shipilina, Louis Boell, Melinda Pickup, Yongbiao Xue, Enrico Coen, Nicholas Barton. https://www.biorxiv.org/content/10.1101/2025.02.17.638607v1.full

## Table of Contents
- [Input Data sets](#Input Data sets)
- [Workflow](#methods)
  - [(i) Data, functions, parameters](#i-Data, functions, parameters)
  - [(ii) Genome scan sliding window](#ii-Genome scan sliding window)
  - [(iii) FastClines geographic ](#iii-FastClines geographic cline estimates)
  - [(iv) post processing sliding window and FastClines outputs ](#iv-post processing sliding window and FastClines outputs)
  - [(v) RNAseq DGE ](#v-RNAseq DGE)
  - [(vi) Genotype-Phenotype associations ](#vi-Genotype-Phenotype associations)
  - [(vii) Slow Clines analyses ](#vii-Slow Clines analyses)

# Input Data sets

- Antirrhinum Reference genome (online link)
- Antirrhinum Annotated genome (online link) GWHBJVT00000000.gff
- Whole Genome PoolSeq dataset (online link)
- RNAseq dataset (online link)
- Hybrid zone plant genotypes and spatial dataset (Dryad link)
- Hybrid zone plant colour phenotype dataset (Dryad link)

# Workflow

Analyses and figures for main text and supplementary materials are conducted in a series of R markdown scripts:

- WGS_FastClines_Antirrhinum.Rmd
- RNAseq_Antirrhinum.Rmd
- GenoPheno_Antirrhinum.Rmd
- SlowClines_Antirrhinum.Rmd

and also a number of Python scripts called within the markdowns:

- SlidingWindows.py
- FastClines.py

Most of the main analyses and generation of figures are conducted in the  script R/WGS_FastClines_Antirrhinum.Rmd

The workflow pipeline summary

# (i) Data, functions, parameters

Input datasets:

- Antirrhinum Reference genome 
- Antirrhinum Annotated genome GWHBJVT00000000.gff
- Whole Genome PoolSeq dataset
- RNAseq dataset
- Hybrid zone plant genotypes and spatial dataset
- Hybrid zone plant colour phenotype dataset

Functions: see `Source-my-functions` code chunk

- WGS_FastClines_Antirrhinum.Rmd
- Genome scan sliding window analyses (Fst, Dxy, pi)
- FastClines analysis (cline centre and width)
