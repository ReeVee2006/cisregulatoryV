# Cis-regulatory region variants

---------------------------------

## About
Associated scripts and input data

This github contains the scripts and data tables required for generation of the information and figures in Villani et al. Regional-specific calibration enables application of computational evidence for clinical classification of 5’ cis-regulatory variants in Mendelian disease

## Usage requirements
Script use requires:

R version 4.2.3

R packages; tidyverse 2.0.0, GenomicRanges 1.50.2, ggpubr 0.6.0, cowplot 1.1.1, grid 4.2.3, vcfR 1.14.0, rstatix 0.7.2

Data provided in the data package has been generated using the following packages:

htslib 1.9, bedops 2.4.35, VEP 99.2, SpliceAI 1.3.1, CADD GRCh38-v1.6, FATHMM-MKL MKL, FATHMM-XF XF, REMM V0.4, VEP cmd line tool 99.2, VEP online tool version 109, Biomart Ensembl release 110 - July 2023, gnomAD v3.0 VCF, phyloP100way 2015.05.11, GERP 99.2, Linsight 14 April 2020, Eigen 1.0.1, ClinVar VCF variant files clinvar_20221105, SpliceAI 1.3.1

Additionally data has been sourced from the following online resources:

EPDnew; https://epd.expasy.org/epd

ncVarDB; https://github.com/Gardner-BinfLab/ncVarDB

Alsheikh AJ, Wollenhaupt S, King EA, Reeb J, Ghosh S, Stolzenburg LR, Tamim S, Lazar J, Davis JW, Jacob HJ. The landscape of GWAS validation; systematic review identifying 309 validated non-coding variants across 130 human diseases. BMC Med Genomics. 2022 Apr 1;15(1):74. doi: 10.1186/s12920-022-01216-w. PMID: 35365203; PMCID: PMC8973751.

## Scripts provided
Script_Fig2_ncvardb; Analysis of non-coding variants in ncVarDB impact prediction scores and molecular consequence

Script_Fig3_cont_A; Extraction of data from genomAD VCF for correlation analysis

Script_Fig3_cont_B; Generation of control variant correlation figure from matched regions

Script_cisreg FILTER; Analysis of cis-regulatory region variants, including filter process for collecting variants to create references sets

Script_Fig4_calibration_A; Collation of impact prediction scores and categories for cis-regulatory region variants, allocation of categories

Script_Fig4_calibration_B; Clinical performance evaluation of cis-regulatory region variants based on categories previously allocated

Script_Fig5_heatmap; Heatmap comparison of clinical performance between impact prediction scores, 

Script_Fig5_score_eval; Calculation of category performances for prediction tools, CADD REMM combination.

Script_Fig6_EPDnew; Analysis of reference set location in EPDnew and evaluation of EPDnew incluson in impact prediction

Script_Supp_Fig_dis; Summary features of collected disease-associated cis-regulatory region variants

Script_Supp_Fig_remm; Analysis of reference set location in promoter features

Script_Supp_Fig_promoter; Analysis of reference set location in promoter features

Script_Supp_evalution; Calculation of prediction tool category performance components. 
