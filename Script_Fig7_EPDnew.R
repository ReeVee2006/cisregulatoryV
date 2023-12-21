##cis reg variants - 
##fig 8 - EPDnew analysis
##input datasets in data folder - 
##visualisation and analysis - 21.11.23 RV

library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(grid)
library(cowplot)

#-setup-------------------------------------------
#setwd 
setwd("D:/cisreg_manuscript/SUPP")

## Load datasets
#load disease-associated variant set
cisregDraw = read.table("Data/Supplemental_DisVar.txt", header= TRUE, sep = "\t")
#load control variant set
cisregCraw = read.table("Data/Supplemental_ContVar.txt", header= TRUE, sep = "\t")

#upload annotation data
cisregANNOT_c = read.table("temp/cisreg_testset_ANNOTcategories.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
colnames(cisregANNOT_c)

#Input data summary check------------------------
summary(cisregANNOT_c)

cisregANNOT_c$INFO =  str_replace(cisregANNOT_c$INFO, "Cont", "Control")
cisregANNOT_c$INFO =  str_replace(cisregANNOT_c$INFO, "Dis", "Disease")
cisregANNOT_c$INFO = as.factor(cisregANNOT_c$INFO)
summary(cisregANNOT_c$INFO)

#count the total number of variants in the 
n_distinct(cisregANNOT_c$ID)
#Control n in test set (all variants)
sum(cisregANNOT_c$INFO=="Control")

#disease n in test set (all variants)
sum(cisregANNOT_c$INFO=="Disease")

#-EPDnew-region-annotation----------------------------------------------------------
##Annotation of dataset with EPDnew/promoter region location overlap
##EPDnew set up - HG38
#load EPDnew .bed file with all the EPDnew promoter locations
EPDnew = read.table("Data/Hs_EPDnew.bed", sep = " ", header= FALSE)

#create GRanges object from cisreg  variant list for input to annotations
#converting dataframe into GRanges object - from vcf like file
grcisregANNOT =  makeGRangesFromDataFrame(cisregANNOT_c,
                                          keep.extra.columns=FALSE,
                                          ignore.strand=TRUE,
                                          seqinfo=NULL,
                                          seqnames.field=c("grchr_38"),
                                          start.field=c("POS"),
                                          end.field=c("end.field"),
                                          starts.in.df.are.0based=FALSE)

#convert EPDnew regions file to GRanges object
grEPD <- makeGRangesFromDataFrame(EPDnew,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("V1"),
                                  start.field="V2",
                                  end.field=c("V3"),
                                  strand.field="V6",
                                  starts.in.df.are.0based=TRUE)

#overlap EPDnew and test set to see if test set are located within the BED file
findOverlaps(grEPD, grcisregANNOT, ignore.strand=TRUE)

#overlap EPDnew and cisregdisease est set to see if test set are located within the BED file
findOverlaps(grcisregANNOT, grEPD, ignore.strand=TRUE)
testVinEPD_gr = subsetByOverlaps(grcisregANNOT, grEPD, ignore.strand=TRUE)
testVinEPD_df = as(testVinEPD_gr, "data.frame")
head(testVinEPD_df)

#save annotation file with EPDnew
#write.table(testVinEPD_df, "Cis Reg Annot/Annotation testing/testVinEPD_df_010223.txt", sep = "\t", col.names = TRUE, quote = FALSE)

#ANNOTATION EPDnew - annotate variant file with in EPDnew or not
testVinEPD_df$CHROM = testVinEPD_df$seqnames
testVinEPD_df$CHROM = gsub("chr","", as.factor(testVinEPD_df$CHROM))
testVinEPD_df$chr_pos = paste(testVinEPD_df$CHROM, testVinEPD_df$start, sep = ":")
cisregANNOT_c$EPDnew <- ifelse(cisregANNOT_c$location38 %in% testVinEPD_df$chr_pos, "EPDnew overlap", "no overlap")
as.factor(cisregANNOT_c$EPDnew) %>% summary()

#-EPDnew basic data summary--------------------------------------------------
#
#number in EPDnew
summary(as.factor(cisregANNOT_c$EPDnew))
cisregANNOT_c  %>% count(INFO, EPDnew)

sum(cisregANNOT_c$INFO=="Disease" & cisregANNOT_c$EPDnew == "EPDnew_overlap")
sum(cisregANNOT_c$INFO=="Control" & cisregANNOT_c$EPDnew == "EPDnew_overlap")

###-EPDnew-diseaseVcontrol-summary-------------------------------------

#EPDnew
ggplot(cisregANNOT_c, aes(INFO, fill = EPDnew))+ 
  geom_bar(position = "fill") + labs(y = "# in EPDnew regions")

ggplot(cisregANNOT_c, aes(INFO, fill = EPDnew))+ 
  geom_bar(position = "dodge") + 
  labs(y = "# in EPDnew regions", x = "")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","#999999","#9999CC"), na.value= "#999999")

#calculate and graph proportion
EPDnew_df = cisregANNOT_c %>% group_by(INFO) %>% count(EPDnew)
EPDnew_df$prop_exp.group = ifelse(EPDnew_df$INFO == "Control", EPDnew_df$n/(sum(cisregANNOT_c$INFO=="Control")), EPDnew_df$n/(sum(cisregANNOT_c$INFO=="Disease")))
print(EPDnew_df)

#graph both groups in one plot
ggplot(EPDnew_df, aes(INFO, prop_exp.group, fill = EPDnew))+ 
  geom_col(position = "dodge") + labs(x = "", y = "EPDnew overlap (proportion)") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#666666","#999999","#9999CC"), na.value= "#999999")

#graphing the reg feature overlap with control and disease separately
EPDnew_df_dis = EPDnew_df[EPDnew_df$INFO == "Disease", ] 
EPDnew_df_cont = EPDnew_df[EPDnew_df$INFO == "Control", ]

#ggplot(cisregANNOT_ccont, aes(INFO, fill = EPDnew))+ 
  geom_bar(position = "dodge",colour="#333333") + labs(x = "", y = "# EPDnew region overlap")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","darkslategray","#9999CC"), na.value="#999999")

#ggplot(cisregANNOT_cdis, aes(INFO, fill = EPDnew))+ 
  geom_bar(position = "dodge",colour="#333333") + labs(x = "", y = "# EPDnew region overlap")+
  theme(legend.title = element_blank()) +
  scale_fill_manual(values=c("#F8766D","brown4","#9999CC"), na.value= "#999999")

epdc_prop = ggplot(EPDnew_df_cont, aes(INFO, prop_exp.group, fill = EPDnew))+ 
  geom_col(position = "dodge",colour="#333333") + labs(x = "", y = "EPDnew overlap (proportion)")+
  theme(legend.title = element_blank(),legend.position="none") +
  scale_fill_manual(values=c("#00BFC4","darkslategray","#9999CC"), na.value="#999999")

epdd_prop = ggplot(EPDnew_df_dis, aes(INFO, prop_exp.group,fill = EPDnew))+ 
  geom_col(position = "dodge",colour="#333333") + labs(x = "", y = "EPDnew overlap (proportion)")+
  theme(legend.title = element_blank(),legend.position="none") +
  scale_fill_manual(values=c("#F8766D","brown4","#9999CC"), na.value= "#999999")

ggarrange(caddc, caddd, ncol = 2, nrow = 1)

#EPD summary

EPDsumm = cisregANNOT_c %>% group_by(INFO, EPDnew) %>%
  summarise (count = n()) %>% mutate (perc= count/sum(count))

ggplot(EPDsumm, aes(x = factor(INFO), y = perc, fill = factor(EPDnew))) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) + 
  labs(x = "", y = "EPDnew region overlap (proportion)")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("brown4","#999999","#9999CC"), na.value= "#999999")

ggplot(EPDsumm, aes(x = factor(INFO), y = perc, fill = factor(EPDnew))) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7, color = "black") + 
  labs(x = "", y = "EPDnew region overlap (proportion)")+
  theme(legend.title = element_blank(), text = element_text(size = 14),axis.text.x = element_text(size = 16)) + 
  scale_fill_manual(values = c("brown4","lightgray"))
  
#chi square test to determine if variants associated with EPDnew
EPDnew_df = data.frame(cisregANNOT_c$INFO,cisregANNOT_c$EPDnew)
EPDnew_df = table(cisregANNOT_c$INFO,cisregANNOT_c$EPDnew)
print(EPDnew_df)

print(chisq.test(EPDnew_df))

###LIKELIHOOD_RATIO_CALCULATIONS#############################

# create a summary table for your diagnostic  evaluations
colnames(cisregANNOT_c)
results = cisregANNOT_c %>% select(.,"ID","INFO","EPDnew",
                                   "CADD_RawScore","CADD_PHRED","CADD_scorebin","CADD_testbin",
                                   "REMM_grc38","REMM_score","REMM_scorebin","REMM_testbin")

#basic values for stats calculations
#Control n in test set (all variants)
sum(results$INFO=="Cont")
#disease n in test set (all variants)
sum(results$INFO=="Dis")

#How many variants in test set scored
#Control n
cont_n = sum(results$INFO=="Cont")

#disease n
dis_n = sum(results$INFO=="Dis")

#test n
all_n = nrow(results)

##test1-LR-calculation-EPDnew------------------------------------

#calculate LR for binary EPDnew colocalisation
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test1 = results %>% select(., "INFO","EPDnew")
table(df_test1)
summary(df_test1)

#Calculate n True Pos (Event pos)
test1_TP = df_test1 %>% filter(INFO=="Dis" & EPDnew == "EPDnew_overlap") %>% count()
print(test1_TP)

#calculate n True neg (no event neg)
test1_TN = df_test1 %>% filter(INFO=="Cont" & EPDnew == "NA") %>% count()
print(test1_TN)

#calculate false positive (no event pos)
test1_FP = df_test1 %>% filter(INFO=="Cont" & EPDnew == "EPDnew_overlap") %>% count()
print(test1_FP)

#calculate false negative (event neg)
test1_FN = df_test1 %>% filter(INFO=="Dis" & EPDnew == "NA") %>% count()
print(test1_FN)

# calculate sensitivity
test1_sens = test1_TP/(test1_TP + test1_FN)
print(test1_sens)

# cal specificity
test1_spec = test1_TN/(test1_TN + test1_FP)
print(test1_spec)

# calculate the LR for pos +
test1_LRpos = (test1_TP/dis_n)/(test1_FP/cont_n)
print(test1_LRpos)

#std error of LR+
test1_seLRpos = sqrt(((1/test1_TP)-(1/dis_n))+((1/test1_FP)-(1/cont_n)))

#calc 95% CI for the LR+
#var(logLR) for the LR +
test1_varLRpos = (1/test1_TN-1/cont_n)+(1/test1_FN-1/dis_n)
print(test1_varLRpos)

#calc low CI for the LR+
exp(log(test1_LRpos) - 1.96 * sqrt(test1_varLRpos))

#calc high CI for the LR+
exp(log(test1_LRpos) + 1.96 * sqrt(test1_varLRpos))

# calculate the LR for negative -
#lr_neg(TP, FP, FN, TN)
#(1-sens)/spec
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
test1_LRneg = (test1_FN/dis_n)/(test1_TN/cont_n)
print(test1_LRneg)

#std error of LR
test1_seLRneg = sqrt(((1/test1_FN)-(1/dis_n))+((1/test1_TN)-(1/cont_n)))

#var(logLR) for the LR -
test1_varLRneg = (1/test1_TN-1/cont_n)+(1/test1_FN-1/dis_n)
print(test1_varLRneg)

#calc low CI for the LR-
exp(log(test1_LRneg) - 1.96 * sqrt(test1_varLRneg))

#calc high CI for the LR-
exp(log(test1_LRneg) + 1.96 * sqrt(test1_varLRneg))

##test2-EPDnew-negative-CADDcalibration------------------------
#establishing the groups

results_test2 = results[results$EPDnew == "no_overlap", ]
summary(as.factor(results_test2$EPDnew))

#test2_group1 = EPDnew +ve - LR calculated above
#test2_group2 = CADD and EPDnew negative
#test2_group3 = CADD and EPDnew negative 
#test2_group4 = CADD and EPDnew negative

#creating new groups
colnames(results_test2)
results_test2$test2_scorebin = case_when(
  results_test2$CADD_PHRED >= 10 ~ "test2_CADDhighEPneg",
  between(results_test2$CADD_PHRED,8,10) ~ "test2_CADDintEPneg",
  results_test2$CADD_PHRED <= 8 ~ "test2_CADDlowEPneg"
)

#need to convert this to a table
summary(as.factor(results_test2$test2_scorebin))
summary(as.factor(results_test2$CADD_scorebin))

results_test2 %>% select("INFO","test2_scorebin") %>% table()
results_test2 %>% select("INFO","EPDnew") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test2 = results_test2[!is.na(results_test2$test2_scorebin), ] %>% select(., "INFO","test2_scorebin")

table(df_test2)
summary(df_test2)
#Control n in test set (all variants)
sum(df_test2$INFO=="Cont")
#disease n in test set (all variants)
sum(df_test2$INFO=="Dis")

#How many variants in test set scored
#Control n
test2_cont_n = sum(df_test2$INFO=="Cont")

#disease n
test2_dis_n = sum(df_test2$INFO=="Dis")

#test n
test2_all_n = nrow(df_test2)

#test1-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(INFO=="Dis" & test2_scorebin == "test2_CADDhighEPneg") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(INFO=="Cont" & test2_scorebin == "test2_CADDlowEPneg") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(INFO=="Cont" & test2_scorebin == "test2_CADDhighEPneg") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(INFO=="Dis" & test2_scorebin == "test2_CADDlowEPneg") %>% count()
print(test2_FN)


#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(INFO=="Dis" & test2_scorebin == "test2_CADDintEPneg") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(INFO=="Cont" & test2_scorebin == "test2_CADDintEPneg") %>% count()

# calculate sensitivity
test2_sens = test2_TP/(test2_TP + test2_FN)
print(test2_sens)

# cal specificity
test2_spec = test2_TN/(test2_TN + test2_FP)
print(test2_spec)

#
# calculate the LR for pos +
#lr_pos(TP, FP, FN, TN) - cmd using package
#sens/(1- spec) - interpretation of haldane fitting with wiggins
# LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
test2_LRpos = (test2_TP/test2_dis_n)/(test2_FP/test2_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR -
test2_varLRpos = (1/test2_TN-1/test2_cont_n)+(1/test2_FN-1/test2_dis_n)
print(test2_varLRpos)

#calc low CI for the LR-
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR-
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

# calculate the LR for negative -
#lr_neg(TP, FP, FN, TN)
#(1-sens)/spec
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
test2_LRneg = (test2_FN/test2_dis_n)/(test2_TN/test2_cont_n)
print(test2_LRneg)

#std error of LR
test2_seLRneg = sqrt(((1/test2_FN)-(1/test2_dis_n))+((1/test2_TN)-(1/test2_cont_n)))

#var(logLR) for the LR -
test2_varLRneg = (1/test2_TN-1/test2_cont_n)+(1/test2_FN-1/test2_dis_n)
print(test2_varLRneg)

#calc low CI for the LR-
exp(log(test2_LRneg) - 1.96 * sqrt(test2_varLRneg))

#calc high CI for the LR-
exp(log(test2_LRneg) + 1.96 * sqrt(test2_varLRneg))

#
# calculate the LR for no call
test2_LRnocall = (test2_UNP/test2_dis_n)/(test2_UNN/test2_cont_n)
print(test2_LRnocall)

#var(logLR) for the LR no call
test2_varLRnocall = (1/test2_UNN-1/test2_cont_n)+(1/test2_UNP-1/test2_dis_n)
print(test2_varLRnocall)

#calc low CI for the LR no call
exp(log(test2_LRnocall) - 1.96 * sqrt(test2_varLRnocall))

#calc high CI for the LR no call
exp(log(test2_LRnocall) + 1.96 * sqrt(test2_varLRnocall))

##test3-EPDnew-negative-REMMcalibration############################
#TEST 3
#establishing the groups

results_test3 = results[results$EPDnew == "NA", ]
summary(as.factor(results_test3$EPDnew))

#test3_group1 = EPDnew +ve - LR calculated above
#test3_group2 = REMM and EPDnew negative
#test3_group3 = REMM and EPDnew negative 
#test3_group4 = REMM and EPDnew negative

#creating new groups
colnames(results_test3)
results_test3$test3_scorebin = case_when(
  results_test3$REMM_score >= 0.86 ~ "test3_REMMhighEPneg",
  between(results_test3$REMM_score,0.8,0.86) ~ "test3_REMMintEPneg",
  results_test3$REMM_score <= 0.8 ~ "test3_REMMlowEPneg"
)

#need to convert this to a table
summary(as.factor(results_test3$test3_scorebin))
summary(as.factor(results_test3$REMM_scorebin))

results_test3 %>% select("INFO","test3_scorebin") %>% table()
results_test3 %>% select("INFO","EPDnew") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test3 = results_test3[!is.na(results_test3$test3_scorebin), ] %>% select(., "INFO","test3_scorebin")

table(df_test3)
summary(df_test3)
#Control n in test set (all variants)
sum(df_test3$INFO=="Cont")
#disease n in test set (all variants)
sum(df_test3$INFO=="Dis")

#How many variants in test set scored
#Control n
test3_cont_n = sum(df_test3$INFO=="Cont")

#disease n
test3_dis_n = sum(df_test3$INFO=="Dis")

#test n
test3_all_n = nrow(df_test3)

##test3-LR-calculation-REMM-EPDnew-negative calculation----------------------------------------------

#Calculate n True Pos (Event pos)
test3_TP = df_test3 %>% filter(INFO=="Dis" & test3_scorebin == "test3_REMMhighEPneg") %>% count()
print(test3_TP)

#calculate n True neg (no event neg)
test3_TN = df_test3 %>% filter(INFO=="Cont" & test3_scorebin == "test3_REMMlowEPneg") %>% count()
print(test3_TN)

#calculate false positive (no event pos)
test3_FP = df_test3 %>% filter(INFO=="Cont" & test3_scorebin == "test3_REMMhighEPneg") %>% count()
print(test3_FP)

#calculate false negative (event neg)
test3_FN = df_test3 %>% filter(INFO=="Dis" & test3_scorebin == "test3_REMMlowEPneg") %>% count()
print(test3_FN)


#Calculate Uncertain (event no call)
test3_UNP = df_test3 %>% filter(INFO=="Dis" & test3_scorebin == "test3_REMMintEPneg") %>% count()

#Calculate Uncertain (no event event no call )
test3_UNN = df_test3 %>% filter(INFO=="Cont" & test3_scorebin == "test3_REMMintEPneg") %>% count()

# calculate sensitivity
test3_sens = test3_TP/(test3_TP + test3_FN)
print(test3_sens)

# cal specificity
test3_spec = test3_TN/(test3_TN + test3_FP)
print(test3_spec)

# calculate the LR for pos
test3_LRpos = (test3_TP/test3_dis_n)/(test3_FP/test3_cont_n)
print(test3_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR -
test3_varLRpos = (1/test3_TN-1/test3_cont_n)+(1/test3_FN-1/test3_dis_n)
print(test3_varLRpos)

#calc low CI for the LR-
exp(log(test3_LRpos) - 1.96 * sqrt(test3_varLRpos))

#calc high CI for the LR-
exp(log(test3_LRpos) + 1.96 * sqrt(test3_varLRpos))

# calculate the LR for negative -
#lr_neg(TP, FP, FN, TN)
#(1-sens)/spec
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
test3_LRneg = (test3_FN/test3_dis_n)/(test3_TN/test3_cont_n)
print(test3_LRneg)

#std error of LR
test3_seLRneg = sqrt(((1/test3_FN)-(1/test3_dis_n))+((1/test3_TN)-(1/test3_cont_n)))

#var(logLR) for the LR -
test3_varLRneg = (1/test3_TN-1/test3_cont_n)+(1/test3_FN-1/test3_dis_n)
print(test3_varLRneg)

#calc low CI for the LR-
exp(log(test3_LRneg) - 1.96 * sqrt(test3_varLRneg))

#calc high CI for the LR-
exp(log(test3_LRneg) + 1.96 * sqrt(test3_varLRneg))

#
# calculate the LR for no call
test3_LRnocall = (test3_UNP/test3_dis_n)/(test3_UNN/test3_cont_n)
print(test3_LRnocall)

#var(logLR) for the LR no call
test3_varLRnocall = (1/test3_UNN-1/test3_cont_n)+(1/test3_UNP-1/test3_dis_n)
print(test3_varLRnocall)

#calc low CI for the LR no call
exp(log(test3_LRnocall) - 1.96 * sqrt(test3_varLRnocall))

#calc high CI for the LR no call
exp(log(test3_LRnocall) + 1.96 * sqrt(test3_varLRnocall))



######test4-LR-calculation-combinedREMM+CADD--------------------------------
#test4-groupings
#REMM >0.86 and CADD >10 
#REMM >0.86 or CADD >10
#REMM 0.8-0.86 AND CADD 8-10
#REMM >0.8 or CADD >8 
#REMM >0.8 and CADD >8


#Calculate n True Pos (Event pos)
TP = df_test4 %>% filter(INFO=="Dis" & test4_scorebin == "test4_group2") %>% count()
print(TP)

#calculate n True neg (no event neg)
TN = df_test4 %>% filter(INFO=="Cont" & test4_scorebin == "test4_group5") %>% count()
print(TN)

#calculate false positive (no event pos)
FP = df_test4 %>% filter(INFO=="Cont" & test4_scorebin == "test4_group2") %>% count()
print(FP)

#calculate false negative (event neg)
FN = df_test4 %>% filter(INFO=="Dis" & test4_scorebin == "test4_group5") %>% count()
print(FN)

#Calculate Uncertain (event no call)
UNP = df_test4 %>% filter(INFO=="Dis" & test4_scorebin == "test4_group3") %>% count()
print(UNP)

#Calculate Uncertain (no event event no call )
UNN = df_test4 %>% filter(INFO=="Cont" & test4_scorebin == "test4_group3") %>% count()
print(UNN)

# calculate sensitivity
sens = TP/(TP + FN)
print(sens)

# cal specificity
spec = TN/(TN + FP)
print(spec)

#
# calculate the LR for pos +
#lr_pos(TP, FP, FN, TN) - cmd using package
#sens/(1- spec) - interpretation of haldane fitting with wiggins
# LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
LRpos = (TP/dis_n)/(FP/cont_n)
print(LRpos)

#std error of LR+
seLRpos = sqrt(((1/TP)-(1/dis_n))+((1/FP)-(1/cont_n)))
#calc 95% CI for the LR+

#var(logLR) for the LR -
varLRpos = (1/TN-1/cont_n)+(1/FN-1/dis_n)
print(varLRpos)

#calc low CI for the LR-
exp(log(LRpos) - 1.96 * sqrt(varLRpos))

#calc high CI for the LR-
exp(log(LRpos) + 1.96 * sqrt(varLRpos))

# calculate the LR for negative -
#lr_neg(TP, FP, FN, TN)
#(1-sens)/spec
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
LRneg = (FN/dis_n)/(TN/cont_n)
print(LRneg)

#std error of LR
seLRneg = sqrt(((1/FN)-(1/dis_n))+((1/TN)-(1/cont_n)))

#var(logLR) for the LR -
varLRneg = (1/TN-1/cont_n)+(1/FN-1/dis_n)
print(varLRneg)

#calc low CI for the LR-
exp(log(LRneg) - 1.96 * sqrt(varLRneg))

#calc high CI for the LR-
exp(log(LRneg) + 1.96 * sqrt(varLRneg))

#
# calculate the LR for no call
LRnocall = (UNP/dis_n)/(UNN/cont_n)
print(LRnocall)

#std error of LR
#calc 95% CI for the LR

#var(logLR) for the LR -
varLRnocall = (1/UNN-1/cont_n)+(1/UNP-1/dis_n)
print(varLRnocall)

#calc low CI for the LR-
exp(log(LRnocall) - 1.96 * sqrt(varLRnocall))

#calc high CI for the LR-
exp(log(LRnocall) + 1.96 * sqrt(varLRnocall))

############COMBINING SCORES######

##########-checking-methods-to-combine------------------
#combine scores mathematically, see if calibration is better

###reduced-categories####test_comb1###################
#both high
#CADD>10 and REMM >0.86
#one high
#CADD>10 and REMM 0.8-0.86
#CADD8-10 and REMM >0.86
#both int
#CADD8-10 and REMM 0.8-0.86
#one low
#CADD8-10 and REMM <0.8
#CADD<8 and REMM 0.8-0.86
#both low
#CADD<8 and REMM <0.8
#conflicting
#CADD<8 and REMM >0.86
#CADD>10 and REMM <0.8

colnames(results)
results$test_comb1 = case_when(
  results$CADD_PHRED >= 10 & results$REMM_score >= 0.86 ~ "2high",
  results$CADD_PHRED >= 10 & between(results$REMM_score,0.80,0.86) ~ "1high",
  results$CADD_PHRED >= 10 & results$REMM_score <= 0.80 ~ "conflict",
  between(results$CADD_PHRED,8,10) & results$REMM_score >= 0.86 ~ "1high",
  between(results$CADD_PHRED,8,10) & between(results$REMM_score,0.80,0.86) ~ "int",
  between(results$CADD_PHRED,8,10) & results$REMM_score <= 0.80 ~ "1low",
  results$CADD_PHRED <= 8 & results$REMM_score >= 0.86 ~ "conflict",
  results$CADD_PHRED <= 8 & between(results$REMM_score,0.80,0.86) ~ "1low",
  results$CADD_PHRED <= 8 & results$REMM_score <= 0.80 ~ "2low",
)
summary(as.factor(results$test_comb1))

#need to convert this to a table
summary(as.factor(results$test_comb1))
summary(as.factor(results$CADD_scorebin))

results %>% select("INFO","test_comb1") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_comb1 = results[!is.na(results$test_comb1), ] %>% select(., "INFO","test_comb1")
summary(as.factor(df_comb1$test_comb1))

table(df_comb1)
summary(df_comb1)
#Control n in test set (all variants)
sum(df_comb1$INFO=="Cont")
#disease n in test set (all variants)
sum(df_comb1$INFO=="Dis")

#How many variants in test set scored
#Control n
comb1_cont_n = sum(df_comb1$INFO=="Cont")

#disease n
comb1_dis_n = sum(df_comb1$INFO=="Dis")

#test n
comb1_all_n = nrow(df_comb1)

#test1-evaluation 

#######LR-for-2high--------
#Calculate n True Pos (Event pos)
comb1_TP = df_comb1 %>% filter(INFO=="Dis" & test_comb1 == "2high") %>% count()
print(comb1_TP)

#calculate n True neg (no event neg)
comb1_TN = df_comb1 %>% filter(INFO=="Cont" & test_comb1 == "2low") %>% count()
print(comb1_TN)

#calculate false positive (no event pos)
comb1_FP = df_comb1 %>% filter(INFO=="Cont" & test_comb1 == "2high") %>% count()
print(comb1_FP)

#calculate false negative (event neg)
comb1_FN = df_comb1 %>% filter(INFO=="Dis" & test_comb1 == "2low") %>% count()
print(comb1_FN)

#calculate the variants with both scores intermediate
comb1_UNN = df_comb1 %>% filter(INFO=="Cont" & test_comb1 == "int") %>% count()
print(comb1_UNN)

comb1_UNP = df_comb1 %>% filter(INFO=="Dis" & test_comb1 == "int") %>% count()
print(comb1_UNP)

# calculate the LR for pos +
# LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
comb1_LRpos = (comb1_TP/comb1_dis_n)/(comb1_FP/comb1_cont_n)
print(comb1_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR -
comb1_varLRpos = (1/comb1_TN-1/comb1_cont_n)+(1/comb1_FN-1/comb1_dis_n)
print(comb1_varLRpos)

#calc low CI for the LR-
exp(log(comb1_LRpos) - 1.96 * sqrt(comb1_varLRpos))

#calc high CI for the LR-
exp(log(comb1_LRpos) + 1.96 * sqrt(comb1_varLRpos))

# calculate the LR for negative
#(1-sens)/spec
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
comb1_LRneg = (comb1_FN/comb1_dis_n)/(comb1_TN/comb1_cont_n)
print(comb1_LRneg)

#std error of LR
comb1_seLRneg = sqrt(((1/comb1_FN)-(1/comb1_dis_n))+((1/comb1_TN)-(1/comb1_cont_n)))

#var(logLR) for the LR -
comb1_varLRneg = (1/comb1_TN-1/comb1_cont_n)+(1/comb1_FN-1/comb1_dis_n)
print(comb1_varLRneg)

#calc low CI for the LR-
exp(log(comb1_LRneg) - 1.96 * sqrt(comb1_varLRneg))

#calc high CI for the LR-
exp(log(comb1_LRneg) + 1.96 * sqrt(comb1_varLRneg))

#
# calculate the LR for no call
comb1_LRnocall = (comb1_UNP/comb1_dis_n)/(comb1_UNN/comb1_cont_n)
print(comb1_LRnocall)

#var(logLR) for the LR no call
comb1_varLRnocall = (1/comb1_UNN-1/comb1_cont_n)+(1/comb1_UNP-1/comb1_dis_n)
print(comb1_varLRnocall)

#calc low CI for the LR no call
exp(log(comb1_LRnocall) - 1.96 * sqrt(comb1_varLRnocall))

#calc high CI for the LR no call
exp(log(comb1_LRnocall) + 1.96 * sqrt(comb1_varLRnocall))

##test2-LR-calculation-for-sngle-pos and sngl neg------------------------------------------

#Calculate n True Pos (Event pos)
comb1_TP = df_comb1 %>% filter(INFO=="Dis" & test_comb1 == "1high") %>% count()
print(comb1_TP)

#calculate n True neg (no event neg)
comb1_TN = df_comb1 %>% filter(INFO=="Cont" & test_comb1 == "1low") %>% count()
print(comb1_TN)

#calculate false positive (no event pos)
comb1_FP = df_comb1 %>% filter(INFO=="Cont" & test_comb1 == "1high") %>% count()
print(comb1_FP)

#calculate false negative (event neg)
comb1_FN = df_comb1 %>% filter(INFO=="Dis" & test_comb1 == "1low") %>% count()
print(comb1_FN)

# calculate the LR for pos +
# LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
comb1_LRpos = (comb1_TP/comb1_dis_n)/(comb1_FP/comb1_cont_n)
print(comb1_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR -
comb1_varLRpos = (1/comb1_TN-1/comb1_cont_n)+(1/comb1_FN-1/comb1_dis_n)
print(comb1_varLRpos)

#calc low CI for the LR-
exp(log(comb1_LRpos) - 1.96 * sqrt(comb1_varLRpos))

#calc high CI for the LR-
exp(log(comb1_LRpos) + 1.96 * sqrt(comb1_varLRpos))

# calculate the LR for negative -
#lr_neg(TP, FP, FN, TN)
#(1-sens)/spec
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
comb1_LRneg = (comb1_FN/comb1_dis_n)/(comb1_TN/comb1_cont_n)
print(comb1_LRneg)

#std error of LR
comb1_seLRneg = sqrt(((1/comb1_FN)-(1/comb1_dis_n))+((1/comb1_TN)-(1/comb1_cont_n)))

#var(logLR) for the LR -
comb1_varLRneg = (1/comb1_TN-1/comb1_cont_n)+(1/comb1_FN-1/comb1_dis_n)
print(comb1_varLRneg)

#calc low CI for the LR-
exp(log(comb1_LRneg) - 1.96 * sqrt(comb1_varLRneg))

#calc high CI for the LR-
exp(log(comb1_LRneg) + 1.96 * sqrt(comb1_varLRneg))

###--calculate the LR for the conflicting variants---------------


#Calculate Uncertain (event conflict)
comb1_ConCont = df_comb1 %>% filter(INFO=="Cont" & test_comb1 == "conflict") %>% count()
print(comb1_ConCont)

#Calculate Uncertain (no event conflict)
comb1_ConDis = df_comb1 %>% filter(INFO=="Dis" & test_comb1 == "conflict") %>% count()
print(comb1_ConDis)

#
# calculate the LR for conflict
comb1_LRnocall = (comb1_ConDis/comb1_dis_n)/(comb1_ConCont/comb1_cont_n)
print(comb1_LRnocall)

#var(logLR) for the LR conflict
comb1_varLRnocall = (1/comb1_ConCont-1/comb1_cont_n)+(1/comb1_ConDis-1/comb1_dis_n)
print(comb1_varLRnocall)

#calc low CI for the LR conflict
exp(log(comb1_LRnocall) - 1.96 * sqrt(comb1_varLRnocall))

#calc high CI for the LR conflict
exp(log(comb1_LRnocall) + 1.96 * sqrt(comb1_varLRnocall))


########repeat calculations but with EPDnew negative variants ----------

#need to convert this to a table
summary(as.factor(results$test_comb1))
summary(as.factor(results$CADD_scorebin))

results %>% filter(EPDnew == "NA") %>% select("INFO","test_comb1") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_comb1_2 = results[!is.na(results$test_comb1), ] %>% filter(EPDnew == "NA") %>% select(., "INFO","test_comb1")
summary(as.factor(df_comb1_2$test_comb1))

table(df_comb1_2)
summary(df_comb1_2)
#Control n in test set (all variants)
sum(df_comb1_2$INFO=="Cont")
#disease n in test set (all variants)
sum(df_comb1_2$INFO=="Dis")

#How many variants in test set scored
#Control n
comb1_cont_n = sum(df_comb1_2$INFO=="Cont")

#disease n
comb1_dis_n = sum(df_comb1_2$INFO=="Dis")

#test n
comb1_all_n = nrow(df_comb1_2)

#test1-evaluation 

# calculate sensitivity
comb1_sens = comb1_TP/(comb1_TP + comb1_FN)
print(comb1_sens)

# cal specificity
comb1_spec = comb1_TN/(comb1_TN + comb1_FP)
print(comb1_spec)

#######LR-for-2high--------
#Calculate n True Pos (Event pos)
comb1_TP = df_comb1_2 %>% filter(INFO=="Dis" & test_comb1 == "2high") %>% count()
print(comb1_TP)

#calculate n True neg (no event neg)
comb1_TN = df_comb1_2 %>% filter(INFO=="Cont" & test_comb1 == "2low") %>% count()
print(comb1_TN)

#calculate false positive (no event pos)
comb1_FP = df_comb1_2 %>% filter(INFO=="Cont" & test_comb1 == "2high") %>% count()
print(comb1_FP)

#calculate false negative (event neg)
comb1_FN = df_comb1_2 %>% filter(INFO=="Dis" & test_comb1 == "2low") %>% count()
print(comb1_FN)

#
# calculate the LR for pos +
#lr_pos(TP, FP, FN, TN) - cmd using package
#sens/(1- spec) - interpretation of haldane fitting with wiggins
# LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
comb1_LRpos = (comb1_TP/comb1_dis_n)/(comb1_FP/comb1_cont_n)
print(comb1_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR -
comb1_varLRpos = (1/comb1_TN-1/comb1_cont_n)+(1/comb1_FN-1/comb1_dis_n)
print(comb1_varLRpos)

#calc low CI for the LR-
exp(log(comb1_LRpos) - 1.96 * sqrt(comb1_varLRpos))

#calc high CI for the LR-
exp(log(comb1_LRpos) + 1.96 * sqrt(comb1_varLRpos))

# calculate the LR for negative -
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
comb1_LRneg = (comb1_FN/comb1_dis_n)/(comb1_TN/comb1_cont_n)
print(comb1_LRneg)

#std error of LR
comb1_seLRneg = sqrt(((1/comb1_FN)-(1/comb1_dis_n))+((1/comb1_TN)-(1/comb1_cont_n)))

#var(logLR) for the LR -
comb1_varLRneg = (1/comb1_TN-1/comb1_cont_n)+(1/comb1_FN-1/comb1_dis_n)
print(comb1_varLRneg)

#calc low CI for the LR-
exp(log(comb1_LRneg) - 1.96 * sqrt(comb1_varLRneg))

#calc high CI for the LR-
exp(log(comb1_LRneg) + 1.96 * sqrt(comb1_varLRneg))

#calculate the LR for intermediate variants
#Calculate Uncertain (event no call)
UNP = df_comb1_2 %>% filter(INFO=="Dis" & test_comb1 == "int") %>% count()
print(UNP)

#Calculate Uncertain (no event event no call )
UNN = df_comb1_2 %>% filter(INFO=="Cont" & test_comb1 == "int") %>% count()
print(UNN)

# calculate the LR for no call
comb1_LRnocall = (comb1_UNP/comb1_dis_n)/(comb1_UNN/comb1_cont_n)
print(comb1_LRnocall)

#var(logLR) for the LR no call
comb1_varLRnocall = (1/comb1_UNN-1/comb1_cont_n)+(1/comb1_UNP-1/comb1_dis_n)
print(comb1_varLRnocall)

#calc low CI for the LR no call
exp(log(comb1_LRnocall) - 1.96 * sqrt(comb1_varLRnocall))

#calc high CI for the LR no call
exp(log(comb1_LRnocall) + 1.96 * sqrt(comb1_varLRnocall))

##test2-LR-calculation-for-sngle-pos and sngl neg------------------------------------------

#Calculate n True Pos (Event pos)
comb1_TP = df_comb1_2 %>% filter(INFO=="Dis" & test_comb1 == "1high") %>% count()
print(comb1_TP)

#calculate n True neg (no event neg)
comb1_TN = df_comb1_2 %>% filter(INFO=="Cont" & test_comb1 == "1low") %>% count()
print(comb1_TN)

#calculate false positive (no event pos)
comb1_FP = df_comb1_2 %>% filter(INFO=="Cont" & test_comb1 == "1high") %>% count()
print(comb1_FP)

#calculate false negative (event neg)
comb1_FN = df_comb1_2 %>% filter(INFO=="Dis" & test_comb1 == "1low") %>% count()
print(comb1_FN)

# calculate the LR for pos +
# LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
comb1_LRpos = (comb1_TP/comb1_dis_n)/(comb1_FP/comb1_cont_n)
print(comb1_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR -
comb1_varLRpos = (1/comb1_TN-1/comb1_cont_n)+(1/comb1_FN-1/comb1_dis_n)
print(comb1_varLRpos)

#calc low CI for the LR-
exp(log(comb1_LRpos) - 1.96 * sqrt(comb1_varLRpos))

#calc high CI for the LR-
exp(log(comb1_LRpos) + 1.96 * sqrt(comb1_varLRpos))

# calculate the LR for negative -
#lr_neg(TP, FP, FN, TN)
#(1-sens)/spec
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
comb1_LRneg = (comb1_FN/comb1_dis_n)/(comb1_TN/comb1_cont_n)
print(comb1_LRneg)

#var(logLR) for the LR -
comb1_varLRneg = (1/comb1_TN-1/comb1_cont_n)+(1/comb1_FN-1/comb1_dis_n)
print(comb1_varLRneg)

#calc low CI for the LR-
exp(log(comb1_LRneg) - 1.96 * sqrt(comb1_varLRneg))

#calc high CI for the LR-
exp(log(comb1_LRneg) + 1.96 * sqrt(comb1_varLRneg))

###--calculate the LR for the conflicting variants---------------


#Calculate Uncertain (event conflict)
comb1_ConCont = df_comb1_2 %>% filter(INFO=="Cont" & test_comb1 == "conflict") %>% count()
print(comb1_ConCont)

#Calculate Uncertain (no event conflict)
comb1_ConDis = df_comb1_2 %>% filter(INFO=="Dis" & test_comb1 == "conflict") %>% count()
print(comb1_ConDis)

#
# calculate the LR for conflict
comb1_LRnocall = (comb1_ConDis/comb1_dis_n)/(comb1_ConCont/comb1_cont_n)
print(comb1_LRnocall)

#var(logLR) for the LR conflict
comb1_varLRnocall = (1/comb1_ConCont-1/comb1_cont_n)+(1/comb1_ConDis-1/comb1_dis_n)
print(comb1_varLRnocall)

#calc low CI for the LR conflict
exp(log(comb1_LRnocall) - 1.96 * sqrt(comb1_varLRnocall))

#calc high CI for the LR conflict
exp(log(comb1_LRnocall) + 1.96 * sqrt(comb1_varLRnocall))


######################ALL LR GRAPH#############################
#2 data tables as saved in combined score evaluation
test = read.table("Data/EpdCaddRemm3.txt", header= TRUE, sep = "\t")
colnames(test)
test$group_no = as.numeric(test$group_no)

#graphtest1
ggplot(data = test, aes(x=group2, y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip()

#graphtest2 - prelim all info together
ggplot(data = test, aes(x=group2, y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10")+
  geom_hline(yintercept = 0.23) + 
  geom_hline(yintercept = 0.48) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 2.1) +
  geom_hline(yintercept = 4.3)

#grphtest4
#convert the table to three separate
colnames(test)
test_epd = test[test$score == "EPD", ]
test_epd$group2 = factor(test_epd$group2,levels = c(test_epd$group2))
print(test_epd)
test_CADD = test[test$score == "CADD", ]
test_CADD$group2 = factor(test_CADD$group2,levels = c(test_CADD$group2))
print(test_CADD)
test_REMM = test[test$score == "REMM", ]
test_REMM$group2 = factor(test_REMM$group2,levels = c(test_REMM$group2))
print(test_REMM)
test_CADDREMM = test[test$score == "CADD_REMM", ]
test_CADDREMM$group2 = factor(test_CADDREMM$group2,levels = c(test_CADDREMM$group2))
print(test_CADDREMM)
test_EpdCADDREMM = test[test$score == "EPD_CADD_REMM", ]
test_EpdCADDREMM$group2 = factor(test_EpdCADDREMM$group2,levels = c(test_EpdCADDREMM$group2))
print(test_EpdCADDREMM)


epdplot = ggplot(data = test_epd, aes(x=group2, y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "EPDn", y = "Likelihood ratio log10 (95% CI)")+
  geom_hline(yintercept = 0.23,  colour = "darkgreen") +
  annotation_custom(grobTree(textGrob("moderate", x=0.08,  y=0.05, hjust=0,gp=gpar(col="darkgreen", fontsize=8))))+
  geom_hline(yintercept = 0.48, colour = "green") +
  annotation_custom(grobTree(textGrob("supporting", x=0.20,  y=0.05, hjust=0,gp=gpar(col="chartreuse3", fontsize=8))))+
  geom_hline(yintercept = 1, colour = "gold") +
  geom_hline(yintercept = 2.1, colour = "orange") +
  annotation_custom(grobTree(textGrob("supporting", x=0.57,  y=0.05, hjust=0,gp=gpar(col="orange", fontsize=8))))+
  geom_hline(yintercept = 4.3, colour = "red")+
  annotation_custom(grobTree(textGrob("moderate", x=0.7,  y=0.05, hjust=0,gp=gpar(col="red", fontsize=8))))+
  geom_hline(yintercept = 18.7, colour = "brown4")+
  annotation_custom(grobTree(textGrob("strong", x=0.94,  y=0.05, hjust=0,gp=gpar(col="brown4", fontsize=8))))

caddplot = ggplot(data = test_CADD, aes(x=group2, y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "EPDn/CADD", y = "Likelihood ratio log10 (95% CI)")+
  geom_hline(yintercept = 0.23,  colour = "darkgreen") +
  annotation_custom(grobTree(textGrob("moderate", x=0.08,  y=0.05, hjust=0,gp=gpar(col="darkgreen", fontsize=8))))+
  geom_hline(yintercept = 0.48, colour = "green") +
  annotation_custom(grobTree(textGrob("supporting", x=0.20,  y=0.05, hjust=0,gp=gpar(col="chartreuse3", fontsize=8))))+
  geom_hline(yintercept = 1, colour = "gold") +
  geom_hline(yintercept = 2.1, colour = "orange") +
  annotation_custom(grobTree(textGrob("supporting", x=0.57,  y=0.05, hjust=0,gp=gpar(col="orange", fontsize=8))))+
  geom_hline(yintercept = 4.3, colour = "red")+
  annotation_custom(grobTree(textGrob("moderate", x=0.7,  y=0.05, hjust=0,gp=gpar(col="red", fontsize=8))))+
  geom_hline(yintercept = 18.7, colour = "brown4")+
  annotation_custom(grobTree(textGrob("strong", x=0.94,  y=0.05, hjust=0,gp=gpar(col="brown4", fontsize=8))))

remmplot = ggplot(data = test_REMM, aes(x=group2, y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "EPDn/REMM", y = "Likelihood ratio log10 (95% CI)")+
  geom_hline(yintercept = 0.23,  colour = "darkgreen") +
  annotation_custom(grobTree(textGrob("moderate", x=0.08,  y=0.05, hjust=0,gp=gpar(col="darkgreen", fontsize=8))))+
  geom_hline(yintercept = 0.48, colour = "green") +
  annotation_custom(grobTree(textGrob("supporting", x=0.20,  y=0.05, hjust=0,gp=gpar(col="chartreuse3", fontsize=8))))+
  geom_hline(yintercept = 1, colour = "gold") +
  geom_hline(yintercept = 2.1, colour = "orange") +
  annotation_custom(grobTree(textGrob("supporting", x=0.57,  y=0.05, hjust=0,gp=gpar(col="orange", fontsize=8))))+
  geom_hline(yintercept = 4.3, colour = "red")+
  annotation_custom(grobTree(textGrob("moderate", x=0.7,  y=0.05, hjust=0,gp=gpar(col="red", fontsize=8))))+
  geom_hline(yintercept = 18.7, colour = "brown4")+
  annotation_custom(grobTree(textGrob("strong", x=0.94,  y=0.05, hjust=0,gp=gpar(col="brown4", fontsize=8))))

caddremmplot = ggplot(data = test_CADDREMM, aes(x=group2, y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "CADD/REMM", y = "Likelihood ratio log10 (95% CI)")+
  geom_hline(yintercept = 0.23,  colour = "darkgreen") +
  annotation_custom(grobTree(textGrob("moderate", x=0.08,  y=0.05, hjust=0,gp=gpar(col="darkgreen", fontsize=8))))+
  geom_hline(yintercept = 0.48, colour = "green") +
  annotation_custom(grobTree(textGrob("supporting", x=0.20,  y=0.05, hjust=0,gp=gpar(col="chartreuse3", fontsize=8))))+
  geom_hline(yintercept = 1, colour = "gold") +
  geom_hline(yintercept = 2.1, colour = "orange") +
  annotation_custom(grobTree(textGrob("supporting", x=0.57,  y=0.05, hjust=0,gp=gpar(col="orange", fontsize=8))))+
  geom_hline(yintercept = 4.3, colour = "red")+
  annotation_custom(grobTree(textGrob("moderate", x=0.7,  y=0.05, hjust=0,gp=gpar(col="red", fontsize=8))))+
  geom_hline(yintercept = 18.7, colour = "brown4")+
  annotation_custom(grobTree(textGrob("strong", x=0.94,  y=0.05, hjust=0,gp=gpar(col="brown4", fontsize=8))))

epdcaddremmplot = ggplot(data = test_EpdCADDREMM, aes(x= group2, y =LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange(stat = "identity")+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "EPDn/CADD/REMM", y = "Likelihood ratio log10 (95% CI)")+
  geom_hline(yintercept = 0.23,  colour = "darkgreen") +
  annotation_custom(grobTree(textGrob("moderate", x=0.08,  y=0.05, hjust=0,gp=gpar(col="darkgreen", fontsize=8))))+
  geom_hline(yintercept = 0.48, colour = "green") +
  annotation_custom(grobTree(textGrob("supporting", x=0.20,  y=0.05, hjust=0,gp=gpar(col="chartreuse3", fontsize=8))))+
  geom_hline(yintercept = 1, colour = "gold") +
  geom_hline(yintercept = 2.1, colour = "orange") +
  annotation_custom(grobTree(textGrob("supporting", x=0.57,  y=0.05, hjust=0,gp=gpar(col="orange", fontsize=8))))+
  geom_hline(yintercept = 4.3, colour = "red")+
  annotation_custom(grobTree(textGrob("moderate", x=0.7,  y=0.05, hjust=0,gp=gpar(col="red", fontsize=8))))+
  geom_hline(yintercept = 18.7, colour = "brown4")+
  annotation_custom(grobTree(textGrob("strong", x=0.94,  y=0.05, hjust=0,gp=gpar(col="brown4", fontsize=8))))


#final figure

epdc_prop_legend = get_legend(ggplot(EPDnew_df_cont, aes(INFO, prop_exp.group, fill = EPDnew))+ 
                                geom_col(position = "dodge",colour="#333333") + labs(x = "", y = "EPDnew overlap (proportion)")+
                                theme(legend.position="right") +
                                scale_fill_manual(name = "Control",values=c("#00BFC4","darkslategray","#9999CC"), na.value="#999999"))

epdd_prop_legend = get_legend(ggplot(EPDnew_df_dis, aes(INFO, prop_exp.group,fill = EPDnew))+ 
                                geom_col(position = "dodge",colour="#333333") + labs(x = "", y = "EPDnew overlap (proportion)")+
                                theme(legend.position="right") +
                                scale_fill_manual(name = "Disease",values=c("#F8766D","brown4","#9999CC"), na.value= "#999999"))

plot_grid(plot_grid(epdc_prop,epdd_prop,(plot_grid
                                         (epdc_prop_legend,epdd_prop_legend, ncol = 1, nrow = 2,rel_widths = c(2,2,1))),
                    labels = c("A"),ncol = 3, nrow = 1),
          plot_grid(epdplot, caddplot, remmplot, epdcaddremmplot,
          labels = c("B","C","D","E"),
          align = "v", rel_widths = c(1, 1, 1,1,1),
          rel_heights = c(4,4,4,6),
          ncol = 1, nrow = 4),
          ncol = 1, nrow = 2, rel_heights = c(2,5))

ggsave(path = "D:/cisreg_manuscript/figures", filename = "Fig_allLR_scores.png", width = 8, height = 12, device='png', dpi=1200)

