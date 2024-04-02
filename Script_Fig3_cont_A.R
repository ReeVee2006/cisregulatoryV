######################################################
##Analysing the correlation between allele frequency (maxAF) and conservation in regulatory regions
#process - R code takes VEP annotated gnomAD control variants (grc38) to determine correlation and graph bins of maxAF 
#input - VEP annotated variants in regulatory ROIs from gnomAD files (grc38), maxAF calculated from designated pops, and used as
#Version 1 - created 03.03.23 

#######################################################
#set up environment

#install required packages("tidyverse")("GenomicRanges")("vcfR")("rstatix)#

library(tidyverse)
library(vcfR)
library(GenomicRanges)
library(rstatix)

#set wd

#load data - VEP output 
#full control variant maxAF file
vepcons = read.vcfR("ts23_afcons_merge.vcf.gz")

#---------------------------------------------
#extracting annotations

vep_info2 = data.frame(vepcons@fix)
vep_info2$maxAF = gsub(".*maxAF=([^;]+)[;].*", "\\1", vep_info2$INFO)

vep_info = vcfR2tidy(vepcons, info_only = TRUE)$fix
vep_info_2 = vep_info %>% 
  mutate(CSQ = strsplit(as.character(CSQ),",")) %>% unnest(CSQ)
n_col = max(stringr::str_count(vep_info_2$CSQ, "|" )) +1
vep_info_3 = cbind(vep_info_2,(stringi::stri_split_fixed(vep_info_2$CSQ, "|", simplify = TRUE, omit_empty = NA)))

#check on content, count number of distinct variants
vep_info_3$ID %>% as.factor() %>% n_distinct

#extracting info column names
vep_header = data.frame(vepcons@meta)
#identification of the CSQ INFO column
infoline = str_which(vep_header$vepcons.meta, "CSQ")

#
infocsq1 = print(gsub("##","", vep_header[infoline,]))
infocsq = as.vector(stringi::stri_split_fixed(infocsq1, "|", simplify = TRUE))
#replace the first b value with simply "Allele' currently includes the whole start of the INFO column

#how to automate the calculation of the start and end of 'new column names)
startcol = ncol(vep_info)+1
endcol = ncol(vep_info_3)
colnames(vep_info_3)[startcol:endcol] = c(infocsq)
colnames(vep_info_3)[colnames(vep_info_3) == "phyloP100\">" ] ="c_phyloP100"

#clean up
colnames(vep_info_3)
cons = vep_info_3 %>% dplyr::select(1:7,"Conservation","c_phyloP100") %>% unique()
afcons = cbind(cons,vep_info2$maxAF)
colnames(afcons)[colnames(afcons) == "vep_info2$maxAF" ] ="maxAF"

#clear environment
#rm(endcol, startcol, n_col, veptxt, infocsq,infocsq1,vep_info, vep_info_2, vepannot)
#rm(vep_info_3, vep_info2, vep_info_2)

#---------------------------------------
##########################
#correlation of AF v conservation

sink("230323_afvcons.txt", append = TRUE)

#prepare the data for analysis and graphing
colnames(afcons)
afcons %>% n_distinct()
summary(afcons)

#select only the rows with filter category PASS
nrow(afcons[afcons$maxAF == "PASS", ] )
afcons = afcons[afcons$FILTER =="PASS", ]
#count number of variants with maxAF not calcucated ie "."
nrow(afcons[afcons$maxAF == ".", ] )
#exclude variants with maxAF not calcucated ie "."
afcons = afcons[afcons$maxAF != ".", ]
#remove rows with NA in maxAF column, rows with phylop or cons NA will remain, and excluded in plots
afcons = afcons[!is.na(afcons$maxAF),]
#convert to numeric class
afcons$maxAF = as.numeric(afcons$maxAF)
afcons$Conservation = as.numeric(afcons$Conservation)
afcons$c_phyloP100 = as.numeric(afcons$c_phyloP100)
summary(afcons)

#set wd for output of all of the results, make a date folder to save specific results event

#
binvalues = c(0, 0.00001,0.00002,0.0001,0.001,0.01,0.05,1)
binlabels = c("0-0.00001","0.00001-0.00002","0.00002-0.0001", "0.0001-0.001", "0.001-0.01", "0.01-0.05", "0.05-1")
afcons$bin  = cut(afcons$maxAF, breaks = binvalues, include.lowest = TRUE, labels = binlabels)
afcons %>% group_by(bin) %>% summary()

write.table(afcons, "temp/testset_vepafcons.txt",sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

cor.test(afcons$maxAF,afcons$c_phyloP100, method = 'spearman')
cor.test(afcons$maxAF,afcons$Conservation, method = 'spearman')

print("binned analysis")

print("summary stats maxAF")
group_by(afcons, bin) %>%
  summarise(
    count = n(),
    mean = mean(maxAF, na.rm = TRUE),
    sd = sd(maxAF, na.rm = TRUE),
    median = median(maxAF, na.rm = TRUE),
    IQR = IQR(maxAF, na.rm = TRUE)
  )

print("summary stats phylop")
group_by(afcons, bin) %>%
  summarise(
    count = n(),
    mean = mean(c_phyloP100, na.rm = TRUE),
    sd = sd(c_phyloP100, na.rm = TRUE),
    median = median(c_phyloP100, na.rm = TRUE),
    IQR = IQR(c_phyloP100, na.rm = TRUE)
  )

print("summary stats Conservation")
group_by(afcons, bin) %>%
  summarise(
    count = n(),
    mean = mean(Conservation, na.rm = TRUE),
    sd = sd(Conservation, na.rm = TRUE),
    median = median(Conservation, na.rm = TRUE),
    IQR = IQR(Conservation, na.rm = TRUE)
  )

ungroup(afcons)

kruskal.test(c_phyloP100 ~ bin, data = afcons)
kruskal.test(Conservation ~ bin, data = afcons)

dev(off)
