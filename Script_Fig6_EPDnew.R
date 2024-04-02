##cis reg variants
##fig 6 - EPDnew analysis
##input datasets in data folder
##visualisation and analysis
## Created by Rehan Villani - 21.03.24

library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(grid)
library(cowplot)

#-setup-------------------------------------------
#setwd 

## Load datasets
#load disease-associated variant set
cisregDraw = read.table("Data/Supplemental_DisVar.txt", header= TRUE, sep = "\t")
#load control variant set
cisregCraw = read.table("Data/Supplemental_ContVar.txt", header= TRUE, sep = "\t")

#load complete annotation file, cisreg_annot for continued  analysis
cisregANNOT_c = read.table("Data/cisreg_testset_ANNOTcategories.txt", stringsAsFactors = TRUE , sep = "\t", header = TRUE)
#thresholds - CADD defined based on the independent calibration via wiggins scores
CADD_upper_threshold = 10
CADD_lower_threshold = 8
#thresholds - defined based on the independent calibration via wiggins scores
REMM_upper_threshold = 0.86
REMM_lower_threshold = 0.80

####################basic reference set characteristics
head(cisregANNOT_c)
cisregANNOT_c$Group = as.factor(ifelse(cisregANNOT_c$INFO == "Cont", "Control", "Disease"))

sum(is.na(cisregANNOT_c$CADD_PHRED))
sum(is.na(cisregANNOT_c$CADD_scorebin))

sum(is.na(cisregANNOT_c$REMM_score))
sum(is.na(cisregANNOT_c$REMM_scorebin))

#count the total number of variants in the 
n_distinct(cisregANNOT_c$ID)

##################setting up LR calculations
# create a summary table for your diagnostic  evaluations
colnames(cisregANNOT_c)
results = cisregANNOT_c %>% select(.,"ID","INFO","Group",
                                   "CADD_RawScore","CADD_PHRED","CADD_scorebin","CADD_testbin",
                                   "REMM_grc38","REMM_score","REMM_scorebin","REMM_testbin")

#basic values for stats calculations
#Control n in test set (all variants)
sum(results$INFO=="Cont")
sum(cisregANNOT_c$Group=="Control")
#disease n in test set (all variants)
sum(results$INFO=="Dis")
sum(cisregANNOT_c$Group=="Disease")

#How many variants in test set scored
#Control n
cont_n = sum(results$Group=="Control")

#disease n
dis_n = sum(results$Group=="Disease")

#test n
all_n = nrow(results)

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
cisregANNOT_c  %>% count(Group, EPDnew)

sum(cisregANNOT_c$Group=="Disease" & cisregANNOT_c$EPDnew == "EPDnew overlap")
sum(cisregANNOT_c$Group=="Control" & cisregANNOT_c$EPDnew == "EPDnew overlap")

###-EPDnew-diseaseVcontrol-summary-------------------------------------

#EPDnew
ggplot(cisregANNOT_c, aes(Group, fill = EPDnew))+ 
  geom_bar(position = "fill") + labs(y = "# in EPDnew regions")

ggplot(cisregANNOT_c, aes(Group, fill = EPDnew))+ 
  geom_bar(position = "dodge") + 
  labs(y = "# in EPDnew regions", x = "")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","#999999","#9999CC"), na.value= "#999999")

#calculate and graph proportion
EPDnew_df = cisregANNOT_c %>% group_by(Group) %>% count(EPDnew)
EPDnew_df$prop_exp.group = ifelse(EPDnew_df$Group == "Control", EPDnew_df$n/(sum(cisregANNOT_c$Group=="Control")), EPDnew_df$n/(sum(cisregANNOT_c$Group=="Disease")))
print(EPDnew_df)

#graph both groups in one plot
ggplot(EPDnew_df, aes(Group, prop_exp.group, fill = EPDnew))+ 
  geom_col(position = "dodge") + labs(x = "", y = "EPDnew overlap (proportion)") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#666666","#999999","#9999CC"), na.value= "#999999")

#graphing the reg feature overlap with control and disease separately
EPDnew_df_dis = EPDnew_df[EPDnew_df$Group == "Disease", ] 
EPDnew_df_cont = EPDnew_df[EPDnew_df$Group == "Control", ]


#EPD summary

EPDsumm = cisregANNOT_c %>% group_by(Group, EPDnew) %>%
  summarise (count = n()) %>% mutate (perc= count/sum(count))

ggplot(EPDsumm, aes(x = factor(Group), y = perc, fill = factor(EPDnew))) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) + 
  labs(x = "", y = "EPDnew region overlap (proportion)")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("brown4","#999999","#9999CC"), na.value= "#999999")

ggplot(EPDsumm, aes(x = factor(Group), y = perc, fill = factor(EPDnew))) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7, color = "black") + 
  labs(x = "", y = "EPDnew region overlap (proportion)")+
  theme(legend.title = element_blank(), text = element_text(size = 14),axis.text.x = element_text(size = 16)) + 
  scale_fill_manual(values = c("brown4","lightgray"))
  
#chi square test to determine if variants associated with EPDnew
EPDnew_df = data.frame(cisregANNOT_c$Group,cisregANNOT_c$EPDnew)
EPDnew_df = table(cisregANNOT_c$Group,cisregANNOT_c$EPDnew)
print(EPDnew_df)

print(chisq.test(EPDnew_df))


epdc_prop = ggplot(EPDnew_df_cont, aes(Group, prop_exp.group, fill = EPDnew))+ 
geom_col(position = "dodge",colour="#333333") + labs(x = "", y = "EPDnew overlap (proportion)")+
  theme(legend.title = element_blank(),legend.position="none") +
  scale_fill_manual(values=c("#00BFC4","darkslategray","#9999CC"), na.value="#999999")

epdd_prop = ggplot(EPDnew_df_dis, aes(Group, prop_exp.group,fill = EPDnew))+ 
geom_col(position = "dodge",colour="#333333") + labs(x = "", y = "EPDnew overlap (proportion)")+
  theme(legend.title = element_blank(),legend.position="none") +
  scale_fill_manual(values=c("#F8766D","brown4","#9999CC"), na.value= "#999999")

##########################

# create a summary table for your diagnostic  evaluations
colnames(cisregANNOT_c)
results = cisregANNOT_c %>% select(.,"ID","Group","EPDnew",
                                   "CADD_RawScore","CADD_PHRED","CADD_scorebin","CADD_testbin",
                                   "REMM_grc38","REMM_score","REMM_scorebin","REMM_testbin")

#basic values for stats calculations
#Control n in test set (all variants)
sum(results$Group=="Control")
#disease n in test set (all variants)
sum(results$Group=="Disease")

#How many variants in test set scored
#Control n
cont_n = sum(results$Group=="Control")

#disease n
dis_n = sum(results$Group=="Disease")

#test n
all_n = nrow(results)

####TOOL#PERFORMANCE#SUMMARY###CALCULATIONS#####
#calculate LR for binary EPDnew colocalisation
test = "EPDnew"

#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test1 = results %>% select(., "Group","EPDnew")
table(df_test1)
summary(df_test1)

#Calculate n True Pos (Event pos)
test1_TP = df_test1 %>% filter(Group=="Disease" & EPDnew == "EPDnew overlap") %>% count()
print(test1_TP)

#calculate n True neg (no event neg)
test1_TN = df_test1 %>% filter(Group=="Control" & EPDnew == "no overlap") %>% count()
print(test1_TN)

#calculate false positive (no event pos)
test1_FP = df_test1 %>% filter(Group=="Control" & EPDnew == "EPDnew overlap") %>% count()
print(test1_FP)

#calculate false negative (event neg)
test1_FN = df_test1 %>% filter(Group=="Disease" & EPDnew == "no overlap") %>% count()
print(test1_FN)

# calculate sensitivity
test1_sens = test1_TP/(test1_TP + test1_FN)
print(test1_sens)

# cal specificity
test1_spec = test1_TN/(test1_TN + test1_FP)
print(test1_spec)

df_test1 %>% group_by(Group) %>% summary()

#LIKELIHOOD_RATIO_CALCULATIONS#####
##test1-LR-calculation-EPDnew------------------------------------

# calculate the LR for pos +
test1_LRpos = (test1_TP/dis_n)/(test1_FP/cont_n)
print(test1_LRpos)

#std var LR+
test1_seLRpos = sqrt(((1/test1_TP)-(1/dis_n))+((1/test1_FP)-(1/cont_n)))

#calc 95% CI for the LR+
#var(logLR) for the LR +
test1_varLRpos = (1/test1_TP-1/cont_n)+(1/test1_FP-1/dis_n)
print(test1_varLRpos)

#calc low CI for the LR+
exp(log(test1_LRpos) - 1.96 * sqrt(test1_varLRpos))

#calc high CI for the LR+
exp(log(test1_LRpos) + 1.96 * sqrt(test1_varLRpos))

epdnew_overlap = c(test = test,
                   group = "Epos",
                   LR = as.numeric(test1_LRpos), 
                   Var = as.numeric(test1_varLRpos), 
                   CI_low = as.numeric(exp(log(test1_LRpos) - 1.96 * sqrt(test1_varLRpos))), 
                   CI_high = as.numeric(exp(log(test1_LRpos) + 1.96 * sqrt(test1_varLRpos))))

epdnew_overlap

# calculate the LR for negative -
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

epdnew_noover = c(test = test,
                   group = "Eneg",
                   LR = as.numeric(test1_LRneg), 
                   Var = as.numeric(test1_varLRneg), 
                   CI_low = as.numeric(exp(log(test1_LRneg) - 1.96 * sqrt(test1_varLRneg))), 
                   CI_high = as.numeric(exp(log(test1_LRneg) + 1.96 * sqrt(test1_varLRneg))))

epdnew_noover

##test2-EPDnew-negative-CADDcalibration------------------------
#establishing the groups
results_test2 = results[results$EPDnew == "no overlap", ]
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

results_test2 %>% select("Group","test2_scorebin") %>% table()
results_test2 %>% select("Group","EPDnew") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test2 = results_test2[!is.na(results_test2$test2_scorebin), ] %>% select(., "Group","test2_scorebin")

table(df_test2)
summary(df_test2)
#Control n in test set (all variants)
sum(df_test2$Group=="Cont")
#disease n in test set (all variants)
sum(df_test2$Group=="Disease")

#How many variants in test set scored
#Control n
test2_cont_n = sum(df_test2$Group=="Control")

#Diseaseease n
test2_dis_n = sum(df_test2$Group=="Disease")

#test n
test2_all_n = nrow(df_test2)

#test1-evaluation numbers

##evaluating the LRs for CADD in EPDnew neg variants

#designate test
test = "EPDnew_neg/CADD"

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(Group=="Disease" & test2_scorebin == "test2_CADDhighEPneg") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(Group=="Control" & test2_scorebin == "test2_CADDlowEPneg") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(Group=="Control" & test2_scorebin == "test2_CADDhighEPneg") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(Group=="Disease" & test2_scorebin == "test2_CADDlowEPneg") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(Group=="Disease" & test2_scorebin == "test2_CADDintEPneg") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(Group=="Control" & test2_scorebin == "test2_CADDintEPneg") %>% count()

# calculate sensitivity
test2_sens = test2_TP/(test2_TP + test2_FN)
print(test2_sens)

# cal specificity
test2_spec = test2_TN/(test2_TN + test2_FP)
print(test2_spec)

# calculate the LR for pos +
# LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
test2_LRpos = (test2_TP/test2_dis_n)/(test2_FP/test2_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR +
test2_varLRpos = (1/test2_TP-1/test2_cont_n)+(1/test2_FP-1/test2_dis_n)
print(test2_varLRpos)

#calc low CI for the LR+
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR+
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

#generate summary row for test
EPDnewnegCADDpos = c(test = test,
               group = "CADDpos",
               LR = as.numeric(test2_LRpos), 
               Var = as.numeric(test2_varLRpos), 
               CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
               CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

EPDnewnegCADDpos


# calculate the LR for negative -
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

#generate row information
EPDnewnegCADDneg = c(test = test,
                     group = "CADDneg",
                     LR = as.numeric(test2_LRneg), 
                     Var = as.numeric(test2_varLRneg), 
                     CI_low = as.numeric(exp(log(test2_LRneg) - 1.96 * sqrt(test2_varLRneg))), 
                     CI_high = as.numeric(exp(log(test2_LRneg) + 1.96 * sqrt(test2_varLRneg))))

EPDnewnegCADDneg

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

#generate row information
EPDnewnegCADDint = c(test = test,
                     group = "CADDint",
                     LR = as.numeric(test2_LRnocall), 
                     Var = as.numeric(test2_varLRnocall), 
                     CI_low = as.numeric(exp(log(test2_LRnocall) - 1.96 * sqrt(test2_varLRnocall))), 
                     CI_high = as.numeric(exp(log(test2_LRnocall) + 1.96 * sqrt(test2_varLRnocall))))

EPDnewnegCADDint

##test3-EPDnew-negative-REMMcalibration############################
#TEST 3
#establishing the groups

results_test3 = results[results$EPDnew == "no overlap", ]
summary(as.factor(results_test3$EPDnew))

#designate test
test = "EPDnew_neg/REMM"

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

results_test3 %>% select("Group","test3_scorebin") %>% table()
results_test3 %>% select("Group","EPDnew") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test3 = results_test3[!is.na(results_test3$test3_scorebin), ] %>% select(., "Group","test3_scorebin")

table(df_test3)
summary(df_test3)
#Control n in test set (all variants)
sum(df_test3$Group=="Control")
#disease n in test set (all variants)
sum(df_test3$Group=="Disease")

#How many variants in test set scored
#Control n
test3_cont_n = sum(df_test3$Group=="Control")

#disease n
test3_dis_n = sum(df_test3$Group=="Disease")

#test n
test3_all_n = nrow(df_test3)

##test3-LR-calculation-REMM-EPDnew-negative calculation----------------------------------------------
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic

#Calculate n True Pos (Event pos)
test3_TP = df_test3 %>% filter(Group=="Disease" & test3_scorebin == "test3_REMMhighEPneg") %>% count()
print(test3_TP)

#calculate n True neg (no event neg)
test3_TN = df_test3 %>% filter(Group=="Control" & test3_scorebin == "test3_REMMlowEPneg") %>% count()
print(test3_TN)

#calculate false positive (no event pos)
test3_FP = df_test3 %>% filter(Group=="Control" & test3_scorebin == "test3_REMMhighEPneg") %>% count()
print(test3_FP)

#calculate false negative (event neg)
test3_FN = df_test3 %>% filter(Group=="Disease" & test3_scorebin == "test3_REMMlowEPneg") %>% count()
print(test3_FN)

#Calculate Uncertain (event no call)
test3_UNP = df_test3 %>% filter(Group=="Disease" & test3_scorebin == "test3_REMMintEPneg") %>% count()

#Calculate Uncertain (no event event no call )
test3_UNN = df_test3 %>% filter(Group=="Control" & test3_scorebin == "test3_REMMintEPneg") %>% count()

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
#var(logLR) for the LR +
test3_varLRpos = (1/test3_TP-1/test3_cont_n)+(1/test3_FP-1/test3_dis_n)
print(test3_varLRpos)

#calc low CI for the LR+
exp(log(test3_LRpos) - 1.96 * sqrt(test3_varLRpos))

#calc high CI for the LR-
exp(log(test3_LRpos) + 1.96 * sqrt(test3_varLRpos))

#generate summary row for test
EPDnewnegREMMpos = c(test = test,
                     group = "REMMpos",
                     LR = as.numeric(test2_LRpos), 
                     Var = as.numeric(test2_varLRpos), 
                     CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
                     CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

EPDnewnegREMMpos

# calculate the LR for negative -
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

#generate row information
EPDnewnegREMMneg = c(test = test,
                     group = "REMMneg",
                     LR = as.numeric(test3_LRneg), 
                     Var = as.numeric(test3_varLRneg), 
                     CI_low = as.numeric(exp(log(test3_LRneg) - 1.96 * sqrt(test3_varLRneg))), 
                     CI_high = as.numeric(exp(log(test3_LRneg) + 1.96 * sqrt(test3_varLRneg))))

EPDnewnegREMMneg

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

EPDnewnegREMMint = c(test = test,
                     group = "REMMint",
                     LR = as.numeric(test3_LRnocall), 
                     Var = as.numeric(test3_varLRnocall), 
                     CI_low = as.numeric(exp(log(test3_LRnocall) - 1.96 * sqrt(test3_varLRnocall))), 
                     CI_high = as.numeric(exp(log(test3_LRnocall) + 1.96 * sqrt(test3_varLRnocall))))

EPDnewnegREMMint

######test4-LR-calculation-combinedEPDnew-REMM-CADD--------------------------------
#CADD/REMM categories combined in figure 5, repeat the category allocation to combine with  EPDnew

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

########repeat calculations but with EPDnew negative variants ----------

#need to convert this to a table
summary(as.factor(results$test_comb1))
summary(as.factor(results$CADD_scorebin))

results %>% filter(EPDnew == "no overlap") %>% select("Group","test_comb1") %>% table()

#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_comb1_2 = results[!is.na(results$test_comb1), ] %>% filter(EPDnew == "no overlap") %>% select(., "Group","test_comb1")
summary(as.factor(df_comb1_2$test_comb1))

table(df_comb1_2)
summary(df_comb1_2)

#check that the new EPDnew overlapping removed variant table has correct number of variants
summary(as.factor(results$EPDnew))
nrow(df_comb1_2)

test = "EPDnew_neg/CADD/REMM"

#Control n in test set (all variants)
sum(df_comb1_2$Group=="Control")
#disease n in test set (all variants)
sum(df_comb1_2$Group=="Disease")

#How many variants in test set scored
#Control n
comb1_cont_n = sum(df_comb1_2$Group=="Control")

#disease n
comb1_dis_n = sum(df_comb1_2$Group=="Disease")

#test n
comb1_all_n = nrow(df_comb1_2)

#test1-evaluation 

#######LR-for-2high--------
#Calculate n True Pos (Event pos)
comb1_TP = df_comb1_2 %>% filter(Group=="Disease" & test_comb1 == "2high") %>% count()
print(comb1_TP)

#calculate n True neg (no event neg)
comb1_TN = df_comb1_2 %>% filter(Group=="Control" & test_comb1 == "2low") %>% count()
print(comb1_TN)

#calculate false positive (no event pos)
comb1_FP = df_comb1_2 %>% filter(Group=="Control" & test_comb1 == "2high") %>% count()
print(comb1_FP)

#calculate false negative (event neg)
comb1_FN = df_comb1_2 %>% filter(Group=="Disease" & test_comb1 == "2low") %>% count()
print(comb1_FN)

# calculate the LR for positive +
comb1_LRpos = (comb1_TP/comb1_dis_n)/(comb1_FP/comb1_cont_n)
print(comb1_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR +
comb1_varLRpos = (1/comb1_TP-1/comb1_cont_n)+(1/comb1_FP-1/comb1_dis_n)
print(comb1_varLRpos)

#calc low CI for the LR+
exp(log(comb1_LRpos) - 1.96 * sqrt(comb1_varLRpos))

#calc high CI for the LR+
exp(log(comb1_LRpos) + 1.96 * sqrt(comb1_varLRpos))

EPDnegCADDREMM_2high = c(test = test,
                group = "Both_pos",
                LR = as.numeric(comb1_LRpos), 
                Var = as.numeric(comb1_varLRpos), 
                CI_low = as.numeric(exp(log(comb1_LRpos) - 1.96 * sqrt(comb1_varLRpos))), 
                CI_high = as.numeric(exp(log(comb1_LRpos) + 1.96 * sqrt(comb1_varLRpos))))

EPDnegCADDREMM_2high

#calculation for '1high' group---------------------
group = "One_pos"

#Calculate n True Pos 1 (Event pos - both scores high)
test2_TP1 = df_comb1_2 %>% filter(Group=="Disease" & test_comb1 == "1high") %>% count()
print(test2_TP1)

#calculate false positive (no event pos)
test2_FP1 = df_comb1_2 %>% filter(Group=="Control" & test_comb1 == "1high") %>% count()
print(test2_FP1)

# calculate the LR for pos +
comb1_LRpos = (test2_TP1/comb1_dis_n)/(test2_FP1/comb1_cont_n)
print(comb1_LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
comb1_varLRpos = (1/test2_TP1-1/comb1_cont_n)+(1/test2_FP1-1/comb1_dis_n)
print(comb1_varLRpos)

#calc low CI for the LR-
exp(log(comb1_LRpos) - 1.96 * sqrt(comb1_varLRpos))

#calc high CI for the LR-
exp(log(comb1_LRpos) + 1.96 * sqrt(comb1_varLRpos))

EPDnegCADDREMM_1high = c(test = test,
                group = "One_pos",
                LR = as.numeric(comb1_LRpos), 
                Var = as.numeric(comb1_varLRpos), 
                CI_low = as.numeric(exp(log(comb1_LRpos) - 1.96 * sqrt(comb1_varLRpos))), 
                CI_high = as.numeric(exp(log(comb1_LRpos) + 1.96 * sqrt(comb1_varLRpos))))

EPDnegCADDREMM_1high

#calculation for 'conflict' group---------------------
group = "Conflict"

#Calculate n True Pos 1 (Event pos - both scores high)
comb1__TP1 = df_comb1_2 %>% filter(Group=="Disease" & test_comb1 == "conflict") %>% count()
print(comb1__TP1)

#calculate false positive (no event pos)
comb1__FP1 = df_comb1_2 %>% filter(Group=="Control" & test_comb1 == "conflict") %>% count()
print(comb1__FP1)

# calculate the LR for pos +
comb1__LRpos = (comb1__TP1/comb1_dis_n)/(comb1__FP1/comb1_cont_n)
print(comb1__LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
comb1__varLRpos = (1/comb1__TP1-1/comb1_cont_n)+(1/comb1__FP1-1/comb1_dis_n)
print(comb1__varLRpos)

#calc low CI for the LR-
exp(log(comb1__LRpos) - 1.96 * sqrt(comb1__varLRpos))

#calc high CI for the LR-
exp(log(comb1__LRpos) + 1.96 * sqrt(comb1__varLRpos))

EPDnegCADDREMM_conflict = c(test = test,
                   group = "Conflict",
                   LR = as.numeric(comb1__LRpos), 
                   Var = as.numeric(comb1__varLRpos), 
                   CI_low = as.numeric(exp(log(comb1__LRpos) - 1.96 * sqrt(comb1__varLRpos))), 
                   CI_high = as.numeric(exp(log(comb1__LRpos) + 1.96 * sqrt(comb1__varLRpos))))

EPDnegCADDREMM_conflict

#calculation for 'int' group---------------------
group = "Both_uninf"

#Calculate n True Pos 1 (Event pos - both scores high)
comb1__TP1 = df_comb1_2 %>% filter(Group=="Disease" & test_comb1 == "int") %>% count()
print(comb1__TP1)

#calculate false positive (no event pos)
comb1__FP1 = df_comb1_2 %>% filter(Group=="Control" & test_comb1 == "int") %>% count()
print(comb1__FP1)

# calculate the LR for pos +
comb1__LRpos = (comb1__TP1/comb1_dis_n)/(comb1__FP1/comb1_cont_n)
print(comb1__LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
comb1__varLRpos = (1/comb1__TP1-1/comb1_cont_n)+(1/comb1__FP1-1/comb1_dis_n)
print(comb1__varLRpos)

#calc low CI for the LR-
exp(log(comb1__LRpos) - 1.96 * sqrt(comb1__varLRpos))

#calc high CI for the LR-
exp(log(comb1__LRpos) + 1.96 * sqrt(comb1__varLRpos))

EPDnegCADDREMM_int = c(test = test,
              group = "Both_uninf",
              LR = as.numeric(comb1__LRpos), 
              Var = as.numeric(comb1__varLRpos), 
              CI_low = as.numeric(exp(log(comb1__LRpos) - 1.96 * sqrt(comb1__varLRpos))), 
              CI_high = as.numeric(exp(log(comb1__LRpos) + 1.96 * sqrt(comb1__varLRpos))))

EPDnegCADDREMM_int

#calculation for '1low' group---------------------
group = "One_neg"

#Calculate n True Pos 1 
comb1__TP1 = df_comb1_2 %>% filter(Group=="Disease" & test_comb1 == "1low") %>% count()
print(comb1__TP1)

#calculate false positive (no event pos)
comb1__FP1 = df_comb1_2 %>% filter(Group=="Control" & test_comb1 == "1low") %>% count()
print(comb1__FP1)

# calculate the LR for pos +
comb1__LRpos = (comb1__TP1/comb1_dis_n)/(comb1__FP1/comb1_cont_n)
print(comb1__LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
comb1__varLRpos = (1/comb1__TP1-1/comb1_cont_n)+(1/comb1__FP1-1/comb1_dis_n)
print(comb1__varLRpos)

#calc low CI for the LR-
exp(log(comb1__LRpos) - 1.96 * sqrt(comb1__varLRpos))

#calc high CI for the LR-
exp(log(comb1__LRpos) + 1.96 * sqrt(comb1__varLRpos))

EPDnegCADDREMM_1low = c(test = test,
               group = "One_neg",
               LR = as.numeric(comb1__LRpos), 
               Var = as.numeric(comb1__varLRpos), 
               CI_low = as.numeric(exp(log(comb1__LRpos) - 1.96 * sqrt(comb1__varLRpos))), 
               CI_high = as.numeric(exp(log(comb1__LRpos) + 1.96 * sqrt(comb1__varLRpos))))

EPDnegCADDREMM_1low


#calculation for '2low' group---------------------
group = "Both_neg"

#Calculate n True Pos 1 (Event pos - both scores high)
comb1__TP1 = df_comb1_2 %>% filter(Group=="Disease" & test_comb1 == "2low") %>% count()
print(comb1__TP1)

#calculate false positive (no event pos)
comb1__FP1 = df_comb1_2 %>% filter(Group=="Control" & test_comb1 == "2low") %>% count()
print(comb1__FP1)

# calculate the LR for pos +
comb1__LRpos = (comb1__TP1/comb1_dis_n)/(comb1__FP1/comb1_cont_n)
print(comb1__LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
comb1__varLRpos = (1/comb1__TP1-1/comb1_cont_n)+(1/comb1__FP1-1/comb1_dis_n)
print(comb1__varLRpos)

#calc low CI for the LR-
exp(log(comb1__LRpos) - 1.96 * sqrt(comb1__varLRpos))

#calc high CI for the LR-
exp(log(comb1__LRpos) + 1.96 * sqrt(comb1__varLRpos))

EPDnegCADDREMM_2low = c(test = test,
               group = "Both_neg",
               LR = as.numeric(comb1__LRpos), 
               Var = as.numeric(comb1__varLRpos), 
               CI_low = as.numeric(exp(log(comb1__LRpos) - 1.96 * sqrt(comb1__varLRpos))), 
               CI_high = as.numeric(exp(log(comb1__LRpos) + 1.96 * sqrt(comb1__varLRpos))))

EPDnegCADDREMM_2low

###################################

#join to create a combined table of all tool results
cisregLRs = type.convert(as.data.frame(rbind(
  epdnew_overlap, epdnew_noover,
  EPDnewnegCADDpos,EPDnewnegCADDint,EPDnewnegCADDneg,
  EPDnewnegREMMpos,EPDnewnegREMMint,EPDnewnegREMMneg,
  EPDnegCADDREMM_2high,EPDnegCADDREMM_1high, EPDnegCADDREMM_conflict, EPDnegCADDREMM_int, EPDnegCADDREMM_1low, EPDnegCADDREMM_2low
  ), row.names = TRUE), as.is = TRUE)

cisregLRs
colnames(cisregLRs) = c("annotation","LR_category","LR","Var","CI_low","CI_high")
#nb in SUPP12 score combination (test name) = annotation and combined score category name = LR category  

#save data table, nb information included in Supp table 12 
write.table(cisregLRs,"Data/240321_reg_variant_tool_evaluation_LRs_allcombined.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

######EVALUATION OF COMBINED SCORE GROUPINGS###################
colnames(cisregANNOT_c)
cisregANNOT_c$epd_bin = as.factor(case_when(
  cisregANNOT_c$EPDnew == "EPDnew overlap" ~ "pred_path",
  cisregANNOT_c$EPDnew == "no overlap" ~ "no_call"))
summary(cisregANNOT_c$epd_bin)
sum(is.na(cisregANNOT_c$epd_bin))

#create 6 bins for evidence categories
cisregANNOT_c$epd_testbin = as.factor(paste(cisregANNOT_c$Group,cisregANNOT_c$epd_bin ,sep = "_"))
summary(cisregANNOT_c$epd_testbin)
sum(is.na(cisregANNOT_c$epd_testbin))

epdnew_summary = c(tool = "EPDnew",
                     all_n = nrow(cisregANNOT_c),
                     total_n_cont = sum(cisregANNOT_c$Group=="Control"),
                     total_n_dis = sum(cisregANNOT_c$Group=="Disease"),
                     n_scored = sum(!is.na(cisregANNOT_c$epd_bin)),
                     n_scored_cont = sum(cisregANNOT_c$Group=="Control" & !is.na(cisregANNOT_c$epd_testbin)),
                     n_scored_dis = sum(cisregANNOT_c$Group=="Disease" & !is.na(cisregANNOT_c$epd_testbin)),
                     n_unscored = sum(is.na(cisregANNOT_c$epd_bin)),
                     n_uninformative = sum(cisregANNOT_c$epd_bin == "no_call"),
                     neg_cont_TN = sum(cisregANNOT_c$epd_testbin == "Control_pred_ben"),
                     int_cont_uninf = sum(cisregANNOT_c$epd_testbin == "Control_no_call"),
                     pos_cont_FP = sum(cisregANNOT_c$epd_testbin == "Control_pred_path"),
                     neg_dis_FN = sum(cisregANNOT_c$epd_testbin == "Disease_pred_ben"),
                     int_dis_uninf = sum(cisregANNOT_c$epd_testbin == "Disease_no_call"),
                     pos_dis_TP = sum(cisregANNOT_c$epd_testbin == "Disease_pred_path"),
                     n_correct = sum(cisregANNOT_c$epd_testbin == "Control_pred_ben")+sum(cisregANNOT_c$epd_testbin == "Disease_pred_path"),
                     n_incorrect = sum(cisregANNOT_c$epd_testbin == "Control_pred_path") + sum(cisregANNOT_c$epd_testbin == "Disease_pred_ben"),
                     n_undetermined = sum(cisregANNOT_c$epd_testbin == "Control_no_call" | 
                                            cisregANNOT_c$epd_testbin == "Disease_no_call" |
                                            cisregANNOT_c$epd_testbin == "Cont_NA" |
                                            cisregANNOT_c$epd_testbin == "Dis_NA"),
                     percent_sens = (sum(cisregANNOT_c$epd_testbin == "Disease_pred_path"))/(sum(cisregANNOT_c$Group=="Disease"))*100,
                     percent_spec = sum(cisregANNOT_c$epd_testbin == "Control_pred_ben")/(sum(cisregANNOT_c$Group=="Control"))*100,
                     percent_accuracy = (((sum(cisregANNOT_c$epd_testbin == "Control_pred_ben")+
                                            sum(cisregANNOT_c$epd_testbin == "Disease_pred_path")))/sum(!is.na(cisregANNOT_c$epd_testbin))*100),
                     percent_unscored = sum(is.na(cisregANNOT_c$test_comb1))/nrow(cisregANNOT_c)*100,
                     percent_uninformative = (sum(cisregANNOT_c$epd_testbin == "Control_no_call" | cisregANNOT_c$epd_testbin == "Disease_no_call"))/nrow(cisregANNOT_c)*100,
                     percent_undetermined = (sum(cisregANNOT_c$epd_testbin == "Control_no_call" | 
                                                   cisregANNOT_c$epd_testbin == "Disease_no_call" |
                                                   cisregANNOT_c$epd_testbin == "Cont_NA" |
                                                   cisregANNOT_c$epd_testbin == "Dis_NA"))/nrow(cisregANNOT_c)*100,
                     CLIN_percent_sens = (sum(cisregANNOT_c$epd_testbin == "Disease_pred_path"))/(sum(cisregANNOT_c$Group=="Disease"))*100,
                     CLIN_percent_spec = sum(cisregANNOT_c$epd_testbin == "Control_pred_ben")/(sum(cisregANNOT_c$Group=="Control"))*100,
                     CLIN_percent_accuracy = ((sum(cisregANNOT_c$epd_testbin == "Control_pred_ben")+
                                                 sum(cisregANNOT_c$epd_testbin == "Disease_pred_path"))/(nrow(cisregANNOT_c)))*100,
                     CLIN_percent_correct = (sum(cisregANNOT_c$epd_testbin == "Control_pred_ben")+sum(cisregANNOT_c$epd_testbin == "Disease_pred_path"))/nrow(cisregANNOT_c)*100,
                     CLIN_percent_incorrect = (sum(cisregANNOT_c$epd_testbin == "Control_pred_path") + sum(cisregANNOT_c$epd_testbin == "Disease_pred_ben"))/nrow(cisregANNOT_c)*100,
                     CLIN_percent_undetermined = (sum(cisregANNOT_c$epd_testbin == "Control_no_call" | 
                                                        cisregANNOT_c$epd_testbin == "Disease_no_call" |
                                                        cisregANNOT_c$epd_testbin == "Cont_NA" |
                                                        cisregANNOT_c$epd_testbin == "Dis_NA")/nrow(cisregANNOT_c)*100))
epdnew_summary

cadd_summary = c(tool = "CADD",
                 all_n = nrow(cisregANNOT_c),
                 total_n_cont = sum(cisregANNOT_c$Group=="Control"),
                 total_n_dis = sum(cisregANNOT_c$Group=="Disease"),
                 n_scored = sum(!is.na(cisregANNOT_c$CADD_PHRED)),
                 n_scored_cont = sum(cisregANNOT_c$Group=="Control" & !is.na(cisregANNOT_c$CADD_PHRED)),
                 n_scored_dis = sum(cisregANNOT_c$Group=="Disease" & !is.na(cisregANNOT_c$CADD_PHRED)),
                 n_unscored = sum(is.na(cisregANNOT_c$CADD_PHRED)),
                 n_uninformative = sum(cisregANNOT_c$CADD_testbin == "Cont_score_no_call" | cisregANNOT_c$CADD_testbin == "Dis_score_no_call"),
                 neg_cont_TN = sum(cisregANNOT_c$CADD_testbin == "Cont_score_low"),
                 int_cont_uninf = sum(cisregANNOT_c$CADD_testbin == "Cont_score_no_call"),
                 pos_cont_FP = sum(cisregANNOT_c$CADD_testbin == "Cont_score_high"),
                 neg_dis_FN = sum(cisregANNOT_c$CADD_testbin == "Dis_score_low"),
                 int_dis_uninf = sum(cisregANNOT_c$CADD_testbin == "Dis_score_no_call"),
                 pos_dis_TP = sum(cisregANNOT_c$CADD_testbin == "Dis_score_high"),
                 n_correct = sum(cisregANNOT_c$CADD_testbin == "Cont_score_low")+sum(cisregANNOT_c$CADD_testbin == "Dis_score_high"),
                 n_incorrect = sum(cisregANNOT_c$CADD_testbin == "Cont_score_high") + sum(cisregANNOT_c$CADD_testbin == "Dis_score_low"),
                 n_undetermined = sum(cisregANNOT_c$CADD_testbin == "Cont_score_no_call" | 
                                        cisregANNOT_c$CADD_testbin == "Dis_score_no_call" |
                                        cisregANNOT_c$CADD_testbin == "Cont_NA" |
                                        cisregANNOT_c$CADD_testbin == "Dis_NA"),
                 percent_sens = (sum(cisregANNOT_c$CADD_testbin == "Dis_score_high"))/((sum(cisregANNOT_c$CADD_testbin == "Dis_score_high") + sum(cisregANNOT_c$CADD_testbin == "Dis_score_low")))*100,
                 percent_spec = sum(cisregANNOT_c$CADD_testbin == "Cont_score_low")/(sum(cisregANNOT_c$CADD_testbin == "Cont_score_low") + sum(cisregANNOT_c$CADD_testbin == "Cont_score_high"))*100,
                 percent_acuracy = (((sum(cisregANNOT_c$CADD_testbin == "Cont_score_low")+
                                        sum(cisregANNOT_c$CADD_testbin == "Dis_score_high")))/sum(!is.na(cisregANNOT_c$CADD_PHRED))*100),
                 percent_unscored = sum(is.na(cisregANNOT_c$CADD_PHRED))/nrow(cisregANNOT_c)*100,
                 percent_uninformative = (sum(cisregANNOT_c$CADD_testbin == "Cont_score_no_call" | cisregANNOT_c$CADD_testbin == "Dis_score_no_call"))/nrow(cisregANNOT_c)*100,
                 percent_undetermined = (sum(cisregANNOT_c$CADD_testbin == "Cont_score_no_call" | 
                                               cisregANNOT_c$CADD_testbin == "Dis_score_no_call" |
                                               cisregANNOT_c$CADD_testbin == "Cont_NA" |
                                               cisregANNOT_c$CADD_testbin == "Dis_NA"))/nrow(cisregANNOT_c)*100,
                 CLIN_percent_sens = (sum(cisregANNOT_c$CADD_testbin == "Dis_score_high"))/(sum(cisregANNOT_c$Group=="Disease"))*100,
                 CLIN_percent_spec = sum(cisregANNOT_c$CADD_testbin == "Cont_score_low")/(sum(cisregANNOT_c$Group=="Control"))*100,
                 CLIN_percent_accuracy = (((sum(cisregANNOT_c$CADD_testbin == "Cont_score_low")+
                                              sum(cisregANNOT_c$CADD_testbin == "Dis_score_high")))/all_n)*100,
                 CLIN_percent_correct = (sum(cisregANNOT_c$CADD_testbin == "Cont_score_low")+sum(cisregANNOT_c$CADD_testbin == "Dis_score_high"))/nrow(cisregANNOT_c)*100,
                 CLIN_percent_incorrect = (sum(cisregANNOT_c$CADD_testbin == "Cont_score_high") + sum(cisregANNOT_c$CADD_testbin == "Dis_score_low"))/nrow(cisregANNOT_c)*100,
                 CLIN_percent_undetermined = (sum(cisregANNOT_c$CADD_testbin == "Cont_score_no_call" | 
                                                    cisregANNOT_c$CADD_testbin == "Dis_score_no_call" |
                                                    cisregANNOT_c$CADD_testbin == "Cont_NA" |
                                                    cisregANNOT_c$CADD_testbin == "Dis_NA")/nrow(cisregANNOT_c)*100))


cadd_summary

remm_summary = c(tool = "REMM",
                 all_n = nrow(cisregANNOT_c),
                 total_n_cont = sum(cisregANNOT_c$Group=="Control"),
                 total_n_dis = sum(cisregANNOT_c$Group=="Disease"),
                 n_scored = sum(!is.na(cisregANNOT_c$REMM_score)),
                 n_scored_cont = sum(cisregANNOT_c$Group=="Control" & !is.na(cisregANNOT_c$REMM_score)),
                 n_scored_dis = sum(cisregANNOT_c$Group=="Disease" & !is.na(cisregANNOT_c$REMM_score)),                 
                 n_unscored = sum(is.na(cisregANNOT_c$REMM_score)),
                 n_uninformative = sum(cisregANNOT_c$REMM_testbin == "Cont_score_no_call" | cisregANNOT_c$REMM_testbin == "Dis_score_no_call"),
                 neg_cont_TN = sum(cisregANNOT_c$REMM_testbin == "Cont_score_low"),
                 int_cont_uninf = sum(cisregANNOT_c$REMM_testbin == "Cont_score_no_call"),
                 pos_cont_FP = sum(cisregANNOT_c$REMM_testbin == "Cont_score_high"),
                 neg_dis_FN = sum(cisregANNOT_c$REMM_testbin == "Dis_score_low"),
                 int_dis_uninf = sum(cisregANNOT_c$REMM_testbin == "Dis_score_no_call"),
                 pos_dis_TP = sum(cisregANNOT_c$REMM_testbin == "Dis_score_high"),
                 n_correct = sum(cisregANNOT_c$REMM_testbin == "Cont_score_low")+sum(cisregANNOT_c$REMM_testbin == "Dis_score_high"),
                 n_incorrect = sum(cisregANNOT_c$REMM_testbin == "Cont_score_high") + sum(cisregANNOT_c$REMM_testbin == "Dis_score_low"),
                 n_undetermined = sum(cisregANNOT_c$REMM_testbin == "Cont_score_no_call" | 
                                        cisregANNOT_c$REMM_testbin == "Dis_score_no_call" |
                                        cisregANNOT_c$REMM_testbin == "Cont_NA" |
                                        cisregANNOT_c$REMM_testbin == "Dis_NA"),
                 percent_sens = (sum(cisregANNOT_c$REMM_testbin == "Dis_score_high"))/((sum(cisregANNOT_c$REMM_testbin == "Dis_score_high") + sum(cisregANNOT_c$REMM_testbin == "Dis_score_low")))*100,
                 percent_spec = sum(cisregANNOT_c$REMM_testbin == "Cont_score_low")/(sum(cisregANNOT_c$REMM_testbin == "Cont_score_low") + sum(cisregANNOT_c$REMM_testbin == "Cont_score_high"))*100,
                 percent_acuracy = (((sum(cisregANNOT_c$REMM_testbin == "Cont_score_low")+
                                        sum(cisregANNOT_c$REMM_testbin == "Dis_score_high")))/sum(!is.na(cisregANNOT_c$REMM_score))*100),
                 percent_unscored = sum(is.na(cisregANNOT_c$REMM_score))/nrow(cisregANNOT_c)*100,
                 percent_uninformative = (sum(cisregANNOT_c$REMM_testbin == "Cont_score_no_call" | cisregANNOT_c$REMM_testbin == "Dis_score_no_call"))/nrow(cisregANNOT_c)*100,
                 percent_undetermined = (sum(cisregANNOT_c$REMM_testbin == "Cont_score_no_call" | 
                                               cisregANNOT_c$REMM_testbin == "Dis_score_no_call" |
                                               cisregANNOT_c$REMM_testbin == "Cont_NA" |
                                               cisregANNOT_c$REMM_testbin == "Dis_NA"))/nrow(cisregANNOT_c)*100,
                 CLIN_percent_sens = (sum(cisregANNOT_c$REMM_testbin == "Dis_score_high"))/(sum(cisregANNOT_c$Group=="Disease"))*100,
                 CLIN_percent_spec = sum(cisregANNOT_c$REMM_testbin == "Cont_score_low")/(sum(cisregANNOT_c$Group=="Control"))*100,
                 CLIN_percent_accuracy = (((sum(cisregANNOT_c$REMM_testbin == "Cont_score_low")+
                                              sum(cisregANNOT_c$REMM_testbin == "Dis_score_high")))/(nrow(cisregANNOT_c)))*100,
                 CLIN_percent_correct = (sum(cisregANNOT_c$REMM_testbin == "Cont_score_low")+sum(cisregANNOT_c$REMM_testbin == "Dis_score_high"))/nrow(cisregANNOT_c)*100,
                 CLIN_percent_incorrect = (sum(cisregANNOT_c$REMM_testbin == "Cont_score_high") + sum(cisregANNOT_c$REMM_testbin == "Dis_score_low"))/nrow(cisregANNOT_c)*100,
                 CLIN_percent_undetermined = (sum(cisregANNOT_c$REMM_testbin == "Cont_score_no_call" | 
                                                    cisregANNOT_c$REMM_testbin == "Dis_score_no_call" |
                                                    cisregANNOT_c$REMM_testbin == "Cont_NA" |
                                                    cisregANNOT_c$REMM_testbin == "Dis_NA")/nrow(cisregANNOT_c)*100))


remm_summary

#combining CADD and REMM
cisregANNOT_c$test_comb1 = case_when(
  cisregANNOT_c$CADD_PHRED >= 10 & cisregANNOT_c$REMM_score >= 0.86 ~ "2high",
  cisregANNOT_c$CADD_PHRED >= 10 & between(cisregANNOT_c$REMM_score,0.80,0.86) ~ "1high",
  cisregANNOT_c$CADD_PHRED >= 10 & cisregANNOT_c$REMM_score <= 0.80 ~ "conflict",
  between(cisregANNOT_c$CADD_PHRED,8,10) & cisregANNOT_c$REMM_score >= 0.86 ~ "1high",
  between(cisregANNOT_c$CADD_PHRED,8,10) & between(cisregANNOT_c$REMM_score,0.80,0.86) ~ "int",
  between(cisregANNOT_c$CADD_PHRED,8,10) & cisregANNOT_c$REMM_score <= 0.80 ~ "1low",
  cisregANNOT_c$CADD_PHRED <= 8 & cisregANNOT_c$REMM_score >= 0.86 ~ "conflict",
  cisregANNOT_c$CADD_PHRED <= 8 & between(cisregANNOT_c$REMM_score,0.80,0.86) ~ "1low",
  cisregANNOT_c$CADD_PHRED <= 8 & cisregANNOT_c$REMM_score <= 0.80 ~ "2low",
)

summary(as.factor(cisregANNOT_c$test_comb1))

cisregANNOT_c$test_comb1_bin = as.factor(case_when(
  cisregANNOT_c$test_comb1 == "2high" | cisregANNOT_c$test_comb1 == "1high" ~ "pred_path",
  cisregANNOT_c$test_comb1 == "1low" | cisregANNOT_c$test_comb1 == "int" | cisregANNOT_c$test_comb1 == "conflict"  ~ "no_call",
  cisregANNOT_c$test_comb1 == "2low" ~ "pred_ben"))
summary(cisregANNOT_c$test_comb1_bin)
sum(is.na(cisregANNOT_c$test_comb1_bin))

#create 6 bins for evidence categories
cisregANNOT_c$test_comb1_testbin = as.factor(paste(cisregANNOT_c$Group,cisregANNOT_c$test_comb1_bin ,sep = "_"))
summary(cisregANNOT_c$test_comb1_testbin)
sum(is.na(cisregANNOT_c$test_comb1_testbin))

CADDREMM_summary = c(tool = "CADD_REMM",
                     all_n = nrow(cisregANNOT_c),
                     total_n_cont = sum(cisregANNOT_c$Group=="Control"),
                     total_n_dis = sum(cisregANNOT_c$Group=="Disease"),
                     n_scored = sum(!is.na(cisregANNOT_c$test_comb1)),
                     n_scored_cont = sum(cisregANNOT_c$Group=="Control" & !is.na(cisregANNOT_c$test_comb1)),
                     n_scored_dis = sum(cisregANNOT_c$Group=="Disease" & !is.na(cisregANNOT_c$test_comb1)),
                     n_unscored = sum(is.na(cisregANNOT_c$test_comb1)),
                     n_uninformative = sum(cisregANNOT_c$test_comb1_bin == "no_call"),
                     neg_cont_TN = sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_ben"),
                     int_cont_uninf = sum(cisregANNOT_c$test_comb1_testbin == "Control_no_call"),
                     pos_cont_FP = sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_path"),
                     neg_dis_FN = sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_ben"),
                     int_dis_uninf = sum(cisregANNOT_c$test_comb1_testbin == "Disease_no_call"),
                     pos_dis_TP = sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_path"),
                     n_correct = sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_ben")+sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_path"),
                     n_incorrect = sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_path") + sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_ben"),
                     n_undetermined = sum(cisregANNOT_c$test_comb1_testbin == "Control_no_call" | 
                                            cisregANNOT_c$test_comb1_testbin == "Disease_no_call" |
                                            cisregANNOT_c$test_comb1_testbin == "Cont_NA" |
                                            cisregANNOT_c$test_comb1_testbin == "Dis_NA"),
                     percent_sens = (sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_path"))/((sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_path") + sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_ben")))*100,
                     percent_spec = sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_ben")/(sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_ben") + sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_path"))*100,
                     percent_acuracy = (((sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_ben")+
                                            sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_path")))/sum(!is.na(cisregANNOT_c$test_comb1))*100),
                     percent_unscored = sum(is.na(cisregANNOT_c$test_comb1))/nrow(cisregANNOT_c)*100,
                     percent_uninformative = (sum(cisregANNOT_c$test_comb1_testbin == "Control_no_call" | cisregANNOT_c$test_comb1_testbin == "Disease_no_call"))/nrow(cisregANNOT_c)*100,
                     percent_undetermined = (sum(cisregANNOT_c$test_comb1_testbin == "Control_no_call" | 
                                                   cisregANNOT_c$test_comb1_testbin == "Disease_no_call" |
                                                   cisregANNOT_c$test_comb1_testbin == "Cont_NA" |
                                                   cisregANNOT_c$test_comb1_testbin == "Dis_NA"))/nrow(cisregANNOT_c)*100,
                     CLIN_percent_sens = (sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_path"))/(sum(cisregANNOT_c$Group=="Disease"))*100,
                     CLIN_percent_spec = sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_ben")/(sum(cisregANNOT_c$Group=="Control"))*100,
                     CLIN_percent_accuracy = ((sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_ben")+
                                                 sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_path"))/(nrow(cisregANNOT_c)))*100,
                     CLIN_percent_correct = (sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_ben")+sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_path"))/nrow(cisregANNOT_c)*100,
                     CLIN_percent_incorrect = (sum(cisregANNOT_c$test_comb1_testbin == "Control_pred_path") + sum(cisregANNOT_c$test_comb1_testbin == "Disease_pred_ben"))/nrow(cisregANNOT_c)*100,
                     CLIN_percent_undetermined = (sum(cisregANNOT_c$test_comb1_testbin == "Control_no_call" | 
                                                        cisregANNOT_c$test_comb1_testbin == "Disease_no_call" |
                                                        cisregANNOT_c$test_comb1_testbin == "Cont_NA" |
                                                        cisregANNOT_c$test_comb1_testbin == "Dis_NA")/nrow(cisregANNOT_c)*100))
CADDREMM_summary

#repeating evaluation on EPDnew negative variants

epdneg = cisregANNOT_c[cisregANNOT_c$EPDnew == "no overlap", ]
summary(as.factor(epdneg$EPDnew))

#check remaining bins
summary(epdneg$test_comb1_bin)
sum(is.na(epdneg$test_comb1_bin))
summary(epdneg$test_comb1_testbin)
sum(is.na(epdneg$test_comb1_testbin))

EPDnegcadd_summary = c(tool = "EPDneg_CADD",
                 all_n = nrow(epdneg),
                 total_n_cont = sum(epdneg$Group=="Control"),
                 total_n_dis = sum(epdneg$Group=="Disease"),
                 n_scored = sum(!is.na(epdneg$CADD_PHRED)),
                 n_scored_cont = sum(epdneg$Group=="Control" & !is.na(epdneg$CADD_PHRED)),
                 n_scored_dis = sum(epdneg$Group=="Disease" & !is.na(epdneg$CADD_PHRED)),
                 n_unscored = sum(is.na(epdneg$CADD_PHRED)),
                 n_uninformative = sum(epdneg$CADD_testbin == "Cont_score_no_call" | epdneg$CADD_testbin == "Dis_score_no_call"),
                 neg_cont_TN = sum(epdneg$CADD_testbin == "Cont_score_low"),
                 int_cont_uninf = sum(epdneg$CADD_testbin == "Cont_score_no_call"),
                 pos_cont_FP = sum(epdneg$CADD_testbin == "Cont_score_high"),
                 neg_dis_FN = sum(epdneg$CADD_testbin == "Dis_score_low"),
                 int_dis_uninf = sum(epdneg$CADD_testbin == "Dis_score_no_call"),
                 pos_dis_TP = sum(epdneg$CADD_testbin == "Dis_score_high"),
                 n_correct = sum(epdneg$CADD_testbin == "Cont_score_low")+sum(epdneg$CADD_testbin == "Dis_score_high"),
                 n_incorrect = sum(epdneg$CADD_testbin == "Cont_score_high") + sum(epdneg$CADD_testbin == "Dis_score_low"),
                 n_undetermined = sum(epdneg$CADD_testbin == "Cont_score_no_call" | 
                                        epdneg$CADD_testbin == "Dis_score_no_call" |
                                        epdneg$CADD_testbin == "Cont_NA" |
                                        epdneg$CADD_testbin == "Dis_NA"),
                 percent_sens = (sum(epdneg$CADD_testbin == "Dis_score_high"))/((sum(epdneg$CADD_testbin == "Dis_score_high") + sum(epdneg$CADD_testbin == "Dis_score_low")))*100,
                 percent_spec = sum(epdneg$CADD_testbin == "Cont_score_low")/(sum(epdneg$CADD_testbin == "Cont_score_low") + sum(epdneg$CADD_testbin == "Cont_score_high"))*100,
                 percent_acuracy = (((sum(epdneg$CADD_testbin == "Cont_score_low")+
                                        sum(epdneg$CADD_testbin == "Dis_score_high")))/sum(!is.na(epdneg$CADD_PHRED))*100),
                 percent_unscored = sum(is.na(epdneg$CADD_PHRED))/nrow(epdneg)*100,
                 percent_uninformative = (sum(epdneg$CADD_testbin == "Cont_score_no_call" | epdneg$CADD_testbin == "Dis_score_no_call"))/nrow(epdneg)*100,
                 percent_undetermined = (sum(epdneg$CADD_testbin == "Cont_score_no_call" | 
                                               epdneg$CADD_testbin == "Dis_score_no_call" |
                                               epdneg$CADD_testbin == "Cont_NA" |
                                               epdneg$CADD_testbin == "Dis_NA"))/nrow(epdneg)*100,
                 CLIN_percent_sens = (sum(epdneg$CADD_testbin == "Dis_score_high"))/(sum(epdneg$Group=="Disease"))*100,
                 CLIN_percent_spec = sum(epdneg$CADD_testbin == "Cont_score_low")/(sum(epdneg$Group=="Control"))*100,
                 CLIN_percent_accuracy = (((sum(epdneg$CADD_testbin == "Cont_score_low")+
                                              sum(epdneg$CADD_testbin == "Dis_score_high")))/all_n)*100,
                 CLIN_percent_correct = (sum(epdneg$CADD_testbin == "Cont_score_low")+sum(epdneg$CADD_testbin == "Dis_score_high"))/nrow(epdneg)*100,
                 CLIN_percent_incorrect = (sum(epdneg$CADD_testbin == "Cont_score_high") + sum(epdneg$CADD_testbin == "Dis_score_low"))/nrow(epdneg)*100,
                 CLIN_percent_undetermined = (sum(epdneg$CADD_testbin == "Cont_score_no_call" | 
                                                    epdneg$CADD_testbin == "Dis_score_no_call" |
                                                    epdneg$CADD_testbin == "Cont_NA" |
                                                    epdneg$CADD_testbin == "Dis_NA")/nrow(epdneg)*100))


EPDnegcadd_summary

EPDnegremm_summary = c(tool = "EPDneg_REMM",
                 all_n = nrow(epdneg),
                 total_n_cont = sum(epdneg$Group=="Control"),
                 total_n_dis = sum(epdneg$Group=="Disease"),
                 n_scored = sum(!is.na(epdneg$REMM_score)),
                 n_scored_cont = sum(epdneg$Group=="Control" & !is.na(epdneg$REMM_score)),
                 n_scored_dis = sum(epdneg$Group=="Disease" & !is.na(epdneg$REMM_score)),                 
                 n_unscored = sum(is.na(epdneg$REMM_score)),
                 n_uninformative = sum(epdneg$REMM_testbin == "Cont_score_no_call" | epdneg$REMM_testbin == "Dis_score_no_call"),
                 neg_cont_TN = sum(epdneg$REMM_testbin == "Cont_score_low"),
                 int_cont_uninf = sum(epdneg$REMM_testbin == "Cont_score_no_call"),
                 pos_cont_FP = sum(epdneg$REMM_testbin == "Cont_score_high"),
                 neg_dis_FN = sum(epdneg$REMM_testbin == "Dis_score_low"),
                 int_dis_uninf = sum(epdneg$REMM_testbin == "Dis_score_no_call"),
                 pos_dis_TP = sum(epdneg$REMM_testbin == "Dis_score_high"),
                 n_correct = sum(epdneg$REMM_testbin == "Cont_score_low")+sum(epdneg$REMM_testbin == "Dis_score_high"),
                 n_incorrect = sum(epdneg$REMM_testbin == "Cont_score_high") + sum(epdneg$REMM_testbin == "Dis_score_low"),
                 n_undetermined = sum(epdneg$REMM_testbin == "Cont_score_no_call" | 
                                        epdneg$REMM_testbin == "Dis_score_no_call" |
                                        epdneg$REMM_testbin == "Cont_NA" |
                                        epdneg$REMM_testbin == "Dis_NA"),
                 percent_sens = (sum(epdneg$REMM_testbin == "Dis_score_high"))/((sum(epdneg$REMM_testbin == "Dis_score_high") + sum(epdneg$REMM_testbin == "Dis_score_low")))*100,
                 percent_spec = sum(epdneg$REMM_testbin == "Cont_score_low")/(sum(epdneg$REMM_testbin == "Cont_score_low") + sum(epdneg$REMM_testbin == "Cont_score_high"))*100,
                 percent_acuracy = (((sum(epdneg$REMM_testbin == "Cont_score_low")+
                                        sum(epdneg$REMM_testbin == "Dis_score_high")))/sum(!is.na(epdneg$REMM_score))*100),
                 percent_unscored = sum(is.na(epdneg$REMM_score))/nrow(epdneg)*100,
                 percent_uninformative = (sum(epdneg$REMM_testbin == "Cont_score_no_call" | epdneg$REMM_testbin == "Dis_score_no_call"))/nrow(epdneg)*100,
                 percent_undetermined = (sum(epdneg$REMM_testbin == "Cont_score_no_call" | 
                                               epdneg$REMM_testbin == "Dis_score_no_call" |
                                               epdneg$REMM_testbin == "Cont_NA" |
                                               epdneg$REMM_testbin == "Dis_NA"))/nrow(epdneg)*100,
                 CLIN_percent_sens = (sum(epdneg$REMM_testbin == "Dis_score_high"))/(sum(epdneg$Group=="Disease"))*100,
                 CLIN_percent_spec = sum(epdneg$REMM_testbin == "Cont_score_low")/(sum(epdneg$Group=="Control"))*100,
                 CLIN_percent_accuracy = (((sum(epdneg$REMM_testbin == "Cont_score_low")+
                                              sum(epdneg$REMM_testbin == "Dis_score_high")))/(nrow(epdneg)))*100,
                 CLIN_percent_correct = (sum(epdneg$REMM_testbin == "Cont_score_low")+sum(epdneg$REMM_testbin == "Dis_score_high"))/nrow(epdneg)*100,
                 CLIN_percent_incorrect = (sum(epdneg$REMM_testbin == "Cont_score_high") + sum(epdneg$REMM_testbin == "Dis_score_low"))/nrow(epdneg)*100,
                 CLIN_percent_undetermined = (sum(epdneg$REMM_testbin == "Cont_score_no_call" | 
                                                    epdneg$REMM_testbin == "Dis_score_no_call" |
                                                    epdneg$REMM_testbin == "Cont_NA" |
                                                    epdneg$REMM_testbin == "Dis_NA")/nrow(epdneg)*100))


EPDnegremm_summary

EPDnegCADDREMM_summary = c(tool = "EPDneg_CADD_REMM",
                     all_n = nrow(epdneg),
                     total_n_cont = sum(epdneg$Group=="Control"),
                     total_n_dis = sum(epdneg$Group=="Disease"),
                     n_scored = sum(!is.na(epdneg$test_comb1)),
                     n_scored_cont = sum(epdneg$Group=="Control" & !is.na(epdneg$test_comb1)),
                     n_scored_dis = sum(epdneg$Group=="Disease" & !is.na(epdneg$test_comb1)),
                     n_unscored = sum(is.na(epdneg$test_comb1)),
                     n_uninformative = sum(epdneg$test_comb1_bin == "no_call"),
                     neg_cont_TN = sum(epdneg$test_comb1_testbin == "Control_pred_ben"),
                     int_cont_uninf = sum(epdneg$test_comb1_testbin == "Control_no_call"),
                     pos_cont_FP = sum(epdneg$test_comb1_testbin == "Control_pred_path"),
                     neg_dis_FN = sum(epdneg$test_comb1_testbin == "Disease_pred_ben"),
                     int_dis_uninf = sum(epdneg$test_comb1_testbin == "Disease_no_call"),
                     pos_dis_TP = sum(epdneg$test_comb1_testbin == "Disease_pred_path"),
                     n_correct = sum(epdneg$test_comb1_testbin == "Control_pred_ben")+sum(epdneg$test_comb1_testbin == "Disease_pred_path"),
                     n_incorrect = sum(epdneg$test_comb1_testbin == "Control_pred_path") + sum(epdneg$test_comb1_testbin == "Disease_pred_ben"),
                     n_undetermined = sum(epdneg$test_comb1_testbin == "Control_no_call" | 
                                            epdneg$test_comb1_testbin == "Disease_no_call" |
                                            epdneg$test_comb1_testbin == "Cont_NA" |
                                            epdneg$test_comb1_testbin == "Dis_NA"),
                     percent_sens = (sum(epdneg$test_comb1_testbin == "Disease_pred_path"))/((sum(epdneg$test_comb1_testbin == "Disease_pred_path") + sum(epdneg$test_comb1_testbin == "Disease_pred_ben")))*100,
                     percent_spec = sum(epdneg$test_comb1_testbin == "Control_pred_ben")/(sum(epdneg$test_comb1_testbin == "Control_pred_ben") + sum(epdneg$test_comb1_testbin == "Control_pred_path"))*100,
                     percent_acuracy = (((sum(epdneg$test_comb1_testbin == "Control_pred_ben")+
                                            sum(epdneg$test_comb1_testbin == "Disease_pred_path")))/sum(!is.na(epdneg$test_comb1))*100),
                     percent_unscored = sum(is.na(epdneg$test_comb1))/nrow(epdneg)*100,
                     percent_uninformative = (sum(epdneg$test_comb1_testbin == "Control_no_call" | epdneg$test_comb1_testbin == "Disease_no_call"))/nrow(epdneg)*100,
                     percent_undetermined = (sum(epdneg$test_comb1_testbin == "Control_no_call" | 
                                                   epdneg$test_comb1_testbin == "Disease_no_call" |
                                                   epdneg$test_comb1_testbin == "Cont_NA" |
                                                   epdneg$test_comb1_testbin == "Dis_NA"))/nrow(epdneg)*100,
                     CLIN_percent_sens = (sum(epdneg$test_comb1_testbin == "Disease_pred_path"))/(sum(epdneg$Group=="Disease"))*100,
                     CLIN_percent_spec = sum(epdneg$test_comb1_testbin == "Control_pred_ben")/(sum(epdneg$Group=="Control"))*100,
                     CLIN_percent_accuracy = ((sum(epdneg$test_comb1_testbin == "Control_pred_ben")+
                                                 sum(epdneg$test_comb1_testbin == "Disease_pred_path"))/(nrow(epdneg)))*100,
                     CLIN_percent_correct = (sum(epdneg$test_comb1_testbin == "Control_pred_ben")+sum(epdneg$test_comb1_testbin == "Disease_pred_path"))/nrow(epdneg)*100,
                     CLIN_percent_incorrect = (sum(epdneg$test_comb1_testbin == "Control_pred_path") + sum(epdneg$test_comb1_testbin == "Disease_pred_ben"))/nrow(epdneg)*100,
                     CLIN_percent_undetermined = (sum(epdneg$test_comb1_testbin == "Control_no_call" | 
                                                        epdneg$test_comb1_testbin == "Disease_no_call" |
                                                        epdneg$test_comb1_testbin == "Cont_NA" |
                                                        epdneg$test_comb1_testbin == "Dis_NA")/nrow(epdneg)*100))
EPDnegCADDREMM_summary

summary_all_comb = type.convert(as.data.frame(rbind(
  epdnew_summary, cadd_summary, remm_summary, CADDREMM_summary, 
  EPDnegcadd_summary,EPDnegremm_summary, EPDnegCADDREMM_summary), row.names = TRUE), as.is = TRUE)

summary_all_comb_rounded = summary_all_comb %>% dplyr::mutate(across(where(is.numeric), round, 2))
write.table(summary_all_comb_rounded,"Data/tool_comb_eval.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


######################ALL LR GRAPH#############################
#2 data tables as saved in combined score evaluation
#test = read.table("Data/EpdCaddRemm3.txt", header= TRUE, sep = "\t")

#graphtest
ggplot(data = cisregLRs, aes(x=LR_category, y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip()

#graph
#convert the table to separate annotation groups
colnames(test)
test_epd = cisregLRs[cisregLRs$annotation == "EPDnew", ]
print(test_epd)
test_CADD = cisregLRs[cisregLRs$annotation == "EPDnew_neg/CADD", ]
print(test_CADD)
test_REMM = cisregLRs[cisregLRs$annotation == "EPDnew_neg/REMM", ]
print(test_REMM)
test_EpdCADDREMM = cisregLRs[cisregLRs$annotation == "EPDnew_neg/CADD/REMM", ]
print(test_EpdCADDREMM)


epdplot = ggplot(data = test_epd, aes(x=reorder(LR_category,LR), y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "EPDnew", y = "Likelihood ratio log10 (95% CI)")+
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

caddplot = ggplot(data = test_CADD, aes(x=reorder(LR_category,LR), y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "EPDnew/CADD", y = "Likelihood ratio log10 (95% CI)")+
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

remmplot = ggplot(data = test_REMM, aes(x=reorder(LR_category,LR), y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "EPDnew/REMM", y = "Likelihood ratio log10 (95% CI)")+
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

#this graph is now allocated to Figure 4
#caddremmplot = ggplot(data = test_CADDREMM, aes(x=reorder(LR_category,LR), y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
#  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
#  labs(x = "CADD/REMM", y = "Likelihood ratio log10 (95% CI)")+
#  geom_hline(yintercept = 0.23,  colour = "darkgreen") +
#  annotation_custom(grobTree(textGrob("moderate", x=0.08,  y=0.05, hjust=0,gp=gpar(col="darkgreen", fontsize=8))))+
#  geom_hline(yintercept = 0.48, colour = "green") +
#  annotation_custom(grobTree(textGrob("supporting", x=0.20,  y=0.05, hjust=0,gp=gpar(col="chartreuse3", fontsize=8))))+
#  geom_hline(yintercept = 1, colour = "gold") +
#  geom_hline(yintercept = 2.1, colour = "orange") +
#  annotation_custom(grobTree(textGrob("supporting", x=0.57,  y=0.05, hjust=0,gp=gpar(col="orange", fontsize=8))))+
#  geom_hline(yintercept = 4.3, colour = "red")+
#  annotation_custom(grobTree(textGrob("moderate", x=0.7,  y=0.05, hjust=0,gp=gpar(col="red", fontsize=8))))+
#  geom_hline(yintercept = 18.7, colour = "brown4")+
#  annotation_custom(grobTree(textGrob("strong", x=0.94,  y=0.05, hjust=0,gp=gpar(col="brown4", fontsize=8))))

epdcaddremmplot = ggplot(data = test_EpdCADDREMM, aes(x= reorder(LR_category,LR), y =LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange(stat = "identity")+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "EPDnew/CADD/REMM", y = "Likelihood ratio log10 (95% CI)")+
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

ggplot(data = test_EpdCADDREMM, aes(x= reorder(LR_category,LR), y =LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange(stat = "identity")+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "EPDnew/CADD/REMM", y = "Likelihood ratio log10 (95% CI)")+
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

epdc_prop_legend = get_legend(ggplot(EPDnew_df_cont, aes(Group, prop_exp.group, fill = EPDnew))+ 
                                geom_col(position = "dodge",colour="#333333") + labs(x = "", y = "EPDnew overlap (proportion)")+
                                theme(legend.position="right") +
                                scale_fill_manual(name = "Control",values=c("#00BFC4","darkslategray","#9999CC"), na.value="#999999"))

epdd_prop_legend = get_legend(ggplot(EPDnew_df_dis, aes(Group, prop_exp.group,fill = EPDnew))+ 
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

ggsave(path = "Data/", filename = "Fig6.png", width = 8, height = 14, device='png', dpi=1200)
