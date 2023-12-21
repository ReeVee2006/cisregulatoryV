##cis reg variant - Fig 7
##promoter analysis - CADD web GUI sub-annotations
##R code complemented by base annotation datasets in corresponding Data folder - including source files for project and sources files for annotation, create a folder for the specific annotation set which is directd within specific folder in code 
##CADD annotation analysis - input files via web/GUI based annotation  - GRCh38 v1.6
##21.09.23 RV

#install required packages
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(cowplot)

#--------------------------------------------
#setwd 
setwd("D:/cisreg_manuscript/SUPP")

## Load datasets to annotate. Cisreg variants, control and disease to combine dataset and create into VCF file
#load disease-associated variant set
cisregDraw = read.table("Data/Supplemental_DisVar.txt", header= TRUE, sep = "\t")
#load control variant set
cisregCraw = read.table("Data/Supplemental_ContVar.txt", header= TRUE, sep = "\t")

#current annotation dataset
cisregAll38_vcf = read.table("temp/cisreg_all38_grc38_vcf.txt", sep = "\t", header = TRUE)

#load annotation data
cisregANNOT_c = read.table("temp/cisreg_testset_ANNOTcategories.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
colnames(cisregANNOT_c)

#----------------------------------------------------------------------------------------------
#using the testset annotation to look at the dataset features
CADD = unique(subset(cisregANNOT_c, select = c("variant38","INFO","CADD_minDistTSS","CADD_minDistTSE","CADD_GC","CADD_CpG","CADD_tOverlapMotifs",
                                      "CADD_RemapOverlapTF","CADD_EnsembleRegulatoryFeature","CADD_EncodeDNase.sum","CADD_EncodeDNase.max",
                                     "CADD_RawScore","CADD_PHRED")))
CADD$variant38 %>% n_distinct()
CADD$INFO =  str_replace(CADD$INFO, "Cont", "Control")
CADD$INFO =  str_replace(CADD$INFO, "Dis", "Disease")
CADD$INFO = as.factor(CADD$INFO)
summary(CADD$INFO)

#summary stats of the CADD Phred by group
CADD %>% group_by(INFO) %>% 
  summarize(min = min(CADD_PHRED),
            q1 = quantile(CADD_PHRED, 0.25),
            median = median(CADD_PHRED),
            mean = mean(CADD_PHRED),
            q3 = quantile(CADD_PHRED, 0.75),
            max = max(CADD_PHRED))

#graphing CADD score control versus disease CADD_PHRED - WITH PEJAVER MISSENSE THRESHOLDS
cadd = ggplot(CADD, aes(INFO, CADD_PHRED))
cadd + geom_boxplot() + geom_hline(yintercept = 25.3) + geom_hline(yintercept = 22.7)

sum(CADD$CADD_PHRED >= 22.7)
sum(CADD$CADD_PHRED >= 22.7 & CADD$INFO == "Dis")
sum(CADD$CADD_PHRED >= 22.7 & CADD$INFO == "Dis")/448
sum(CADD$CADD_PHRED <= 22.7)
sum(CADD$CADD_PHRED <= 22.7& CADD$INFO == "Dis")/448

#graphing C v Disease for CADD component scores
CADD$exp.group = as.factor(CADD$INFO)

####GC-------------------------------------------------------
sum(is.na(CADD$CADD_GC))
CADD %>% group_by(INFO) %>% 
  summarize(min = min(CADD_GC),
            q1 = quantile(CADD_GC, 0.25),
            median = median(CADD_GC),
            mean = mean(CADD_GC),
            q3 = quantile(CADD_GC, 0.75),
            max = max(CADD_GC), 
            n = sum(!is.na(CADD_GC)))

#CADD_GC content
gc = ggplot(CADD, aes(INFO,CADD_GC, fill = INFO))+ 
  geom_boxplot() + labs(x = "", y = "GC content")+
  theme(legend.position="none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#619CFF"), na.value= "#999999")

####CpG--------------------------------------------------------
sum(is.na(CADD$CADD_CpG))
CADD %>% group_by(INFO) %>% 
  summarize(min = min(CADD_CpG),
            q1 = quantile(CADD_CpG, 0.25),
            median = median(CADD_CpG),
            mean = mean(CADD_CpG),
            q3 = quantile(CADD_CpG, 0.75),
            max = max(CADD_CpG), 
            n = sum(!is.na(CADD_CpG)))

cpg = ggplot(CADD, aes(INFO,CADD_CpG, fill = INFO))+ 
  geom_boxplot() + labs(x = "", y = "CpG content")+
  theme(legend.position="none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#619CFF"), na.value= "#999999")

####RemapOverlapTF-------------------------------
#summary stats
sum(is.na(CADD$CADD_RemapOverlapTF))

#check idnetified NAs, adjusting so no overlapTF means 0 overlap, convert NA -> 0 to create adj column
CADD$CADD_RemapOverlapTFadj = CADD$CADD_RemapOverlapTF
CADD$CADD_RemapOverlapTFadj[is.na(CADD$CADD_RemapOverlapTFadj)] <- 0
summary(CADD$CADD_RemapOverlapTF)
summary(as.numeric(CADD$CADD_RemapOverlapTFadj))

CADD %>% group_by(INFO) %>% 
  summarize(min = min(CADD_RemapOverlapTFadj),
            q1 = quantile(CADD_RemapOverlapTFadj, 0.25),
            median = median(CADD_RemapOverlapTFadj),
            mean = mean(CADD_RemapOverlapTFadj),
            q3 = quantile(CADD_RemapOverlapTFadj, 0.75),
            max = max(CADD_RemapOverlapTFadj), 
            n = sum(!is.na(CADD_RemapOverlapTFadj)))

#CADD_RemapOverlapTF
ggplot(CADD, aes(INFO,CADD_RemapOverlapTFadj, fill = INFO))+ 
  geom_boxplot() + labs(x = "", y = "RemapOverlapTF adj")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#619CFF"), na.value= "#999999")

remaptf = ggplot(CADD, aes(INFO,CADD_RemapOverlapTFadj, fill = INFO))+ 
  geom_boxplot() + labs(x = "", y = "RemapOverlapTF(log10)")+ scale_y_log10()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#619CFF"), na.value= "#999999")

##EncodeDNase.max
#summary stats
sum(is.na(CADD$CADD_EncodeDNase.max))

#adjusting so no overlapTF means 0 overlap, convert NA -> 0 to create adj column
CADD$CADD_EncodeDNase.maxadj = CADD$CADD_EncodeDNase.max
CADD$CADD_EncodeDNase.maxadj[is.na(CADD$CADD_EncodeDNase.maxadj)] <- 0
summary(CADD$CADD_EncodeDNase.max)
summary(as.numeric(CADD$CADD_EncodeDNase.maxadj))

CADD %>% group_by(INFO) %>% 
  summarize(min = min(CADD_EncodeDNase.maxadj),
            q1 = quantile(CADD_EncodeDNase.maxadj, 0.25),
            median = median(CADD_EncodeDNase.maxadj),
            mean = mean(CADD_EncodeDNase.maxadj),
            q3 = quantile(CADD_EncodeDNase.maxadj, 0.75),
            max = max(CADD_EncodeDNase.maxadj), 
            n = sum(!is.na(CADD_EncodeDNase.maxadj)))

ggplot(CADD, aes(INFO,CADD_EncodeDNase.maxadj, fill = INFO))+ 
  geom_boxplot() + labs(x = "", y = "EncodeDNase.max")+ 
  theme(legend.title = element_blank(), )+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#619CFF"), na.value= "#999999")

#use the adjusted endoce DNAse max (with NAs converted to 0 ie DNase max = 0)
dnasemax = ggplot(CADD, aes(INFO,CADD_EncodeDNase.maxadj, fill = INFO))+ 
  geom_boxplot() + labs(x = "", y = "EncodeDNase max(log10)") + scale_y_log10() +
  theme(legend.title = element_blank(),legend.position="none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#619CFF"), na.value= "#999999")

#####EnsembleRegulatoryFeature------------------------------------------
#CADD$EnsembleRegulatoryFeature = replace(CADD$EnsembleRegulatoryFeature,is.na(CADD$EnsembleRegulatoryFeature),0)
#CADD$EnsembleRegulatoryFeature[is.na(CADD$EnsembleRegulatoryFeature)] = 0
CADD$CADD_EnsembleRegulatoryFeature = as.factor(CADD$CADD_EnsembleRegulatoryFeature)
summary(CADD$CADD_EnsembleRegulatoryFeature)
ggplot(CADD, aes(CADD_EnsembleRegulatoryFeature, fill = INFO))+ 
  geom_bar(position = "fill") + labs(x = "reg feature ", y = "proportion(sourced from CADD)") + 
  theme(axis.text.x = element_text(angle = 45),legend.title = element_blank())

#graph of how many variants overlap with ANY Ensembl reg feature
CADD$EnsRegFeatSUM = ifelse(CADD$CADD_EnsembleRegulatoryFeature == "<NA>","NA","reg feature overlap")
CADD$EnsRegFeatSUM[is.na(CADD$EnsRegFeatSUM)] <- "no overlap"
CADD$EnsRegFeatSUM = as.factor(CADD$EnsRegFeatSUM)
summary(CADD$EnsRegFeatSUM)

##extracting proportions of Reg feature overlap, relative to parent exp.group type - to graph
EnsRF_df = CADD %>% group_by(INFO) %>% count(EnsRegFeatSUM)
EnsRF_df$prop_exp.group = ifelse(EnsRF_df$INFO == "Control", EnsRF_df$n/9251, EnsRF_df$n/441)
EnsRF_df$fillgrp = c(1:4)
print(EnsRF_df)


#EnsRegFeatSUM
CADD  %>% count(INFO, EnsRegFeatSUM)

ggplot(CADD, aes(INFO, fill = EnsRegFeatSUM))+ 
  geom_bar(position = "fill") + labs(y = "Ensembl regulatory region overlap")

ggplot(EnsRF_df, aes(INFO, prop_exp.group, fill = EnsRegFeatSUM))+ 
  geom_col(position = "dodge") + labs(x = "", y = "regulatory feature overlap (proportion)") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#666666","#999999","#9999CC"), na.value= "#333333")

#graphing the reg feature overlap with control and disease separately - proportion
EnsRF_df_dis = EnsRF_df[EnsRF_df$INFO == "Disease", ]
EnsRF_df_cont = EnsRF_df[EnsRF_df$INFO == "Control", ]
EnsRF_df_overlap = EnsRF_df[EnsRF_df$EnsRegFeatSUM == "reg feature overlap", ]

ensoverlap = ggplot(EnsRF_df_overlap, aes(INFO, prop_exp.group, fill = INFO))+ 
  geom_col(colour="#333333") + labs(x = "", y = "reg. feature overlap (proportion)") +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value= "#333333")

legend = get_legend(ggplot(EnsRF_df_overlap, aes(INFO, prop_exp.group, fill = INFO))+ 
                             geom_col(colour="#333333") + labs(x = "", y = "regulatory feature overlap (proportion)"+
                             theme(legend.title = element_blank())+
                             scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value= "#333333")))

plot_grid(plot_grid(gc, cpg,
          remaptf,dnasemax,
          align = "h", axis = "b", rel_widths = c(1, 1),
          labels = c("A","B", "C", "D"),
          ncol = 2, nrow = 2),
          plot_grid(enscont, ensdis,
                    align = "h", axis = "b", rel_widths = c(1, 1),
                    labels = c("E","F"),
                    ncol = 3, nrow = 1),
          rel_heights = c(2, 1),
          ncol = 1, nrow = 2)
ggsave(path = "D:/cisreg_manuscript/figures", filename = "Fig_promoter_full.png", width = 9, height = 12, device='png', dpi=600)

theme_set(theme_gray(base_size = 14))

#final figure

plot_grid(plot_grid(gc, cpg,remaptf,dnasemax, ensoverlap,
                    align = "h", axis = "b", rel_widths = c(1, 1),
                    labels = c("A","B", "C", "D","E"),
                    ncol = 2, nrow = 3))

ggsave(path = "D:/cisreg_manuscript/figures", filename = "Fig7.png", width = 9, height = 12, device='png', dpi=600)

###statistics###############################

#test for matched variance - F test

CADDcont = CADD[CADD$INFO == "Control",]
CADDdis = CADD[CADD$INFO == "Disease",]

#test for population mean difference - 
#statistical test approach determined via density plots 
#QQplots just to check
#CADD_GCcont
qqnorm(CADDcont$CADD_GC, pch = 1, frame = FALSE)
qqline(CADDcont$CADD_GC, col = "steelblue", lwd = 2)
#CADD_GCdis
qqnorm(CADDdis$CADD_GC, pch = 1, frame = FALSE)
qqline(CADDdis$CADD_GC, col = "steelblue", lwd = 2)
#CpCADD_GCont
qqnorm(CADDcont$CADD_CpG, pch = 1, frame = FALSE)
qqline(CADDcont$CADD_CpG, col = "steelblue", lwd = 2)
#CpCADD_GCont
qqnorm(CADDdis$CADD_CpG, pch = 1, frame = FALSE)
qqline(CADDdis$CADD_CpG, col = "steelblue", lwd = 2)
##CADD_RemapOverlapTFcont
qqnorm(CADDcont$CADD_RemapOverlapTFadj, pch = 1, frame = FALSE)
qqline(CADDcont$CADD_RemapOverlapTFadj, col = "steelblue", lwd = 2)
#CADD_RemapOverlapTFcont
qqnorm(CADDdis$CADD_RemapOverlapTFadj, pch = 1, frame = FALSE)
qqline(CADDdis$CADD_RemapOverlapTFadj, col = "steelblue", lwd = 2)
##CADD_EncodeDNase.max
qqnorm(CADDcont$CADD_EncodeDNase.max, pch = 1, frame = FALSE)
qqline(CADDcont$CADD_EncodeDNase.max, col = "steelblue", lwd = 2)
#CADD_EncodeDNase.max
qqnorm(CADDdis$CADD_EncodeDNase.max, pch = 1, frame = FALSE)
qqline(CADDdis$CADD_EncodeDNase.max, col = "steelblue", lwd = 2)

#large sample number means shapiro-wilk not valid, scores however clearly not parametric
#non-parametric  wilcox test used (MannU)
#CADD_GC
wilcox.test( formula = CADD_GC ~ INFO, data = CADD)
wilcox.test( formula = CADD_CpG ~ INFO, data = CADD)
wilcox.test( formula = CADD_RemapOverlapTFadj ~ INFO, data = CADD)
wilcox.test( formula = CADD_EncodeDNase.max ~ INFO, data = CADD)

#chi square test to determine if variants associated with EnsemblRegulatoryRegion feature
ERFnew_df = data.frame(CADD$INFO,CADD$EnsRegFeatSUM)
ERFnew_df = table(CADD$INFO,CADD$EnsRegFeatSUM)
print(ERFnew_df)

print(chisq.test(ERFnew_df))
