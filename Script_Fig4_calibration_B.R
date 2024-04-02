##cis reg bioinformatic tool calibration - second half, LR calculation
##web based input files for annotations
##annotations and categories annotated in Fig4_calibration
##calculate LRs

library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(cowplot)
library(grid)
#library(liftOver)

#set up files--------------------------------------------
#setwd
date = strftime(Sys.Date(),"%y%m%d")

###ASSOCIATED DATA SETS
##raw datasets-Cisreg variants, control and disease to combine dataset and create into VCF file
# disease-associated variant set
#cisregDraw = read.table("Data/Supplemental_DisVar.txt", header= TRUE, sep = "\t")
# control variant set
#cisregCraw = read.table("Data/Supplemental_ContVar.txt", header= TRUE, sep = "\t")
#reference set - from FILTER script
#refset = read.table("temp/cisreg_refset_grc38.txt", header= TRUE, sep = "\t",stringsAsFactors = TRUE)

#load complete annotation file, cisreg_annot for continued  analysis
cisregANNOT_c = read.table("temp/cisreg_testset_ANNOTcategories.txt", stringsAsFactors = TRUE , sep = "\t", header = TRUE)

#thresholds - CADD defined based on the independent calibration via wiggins scores
CADD_upper_threshold = 10
CADD_lower_threshold = 8
#thresholds - defined based on the independent calibration via wiggins scores
linsight_upper_threshold = 0.24
linsight_lower_threshold = 0.16
#thresholds - defined based on the independent calibration via wiggins scores
fathmmlk_upper_threshold = 0.59
fathmmlk_lower_threshold = 0.39
#thresholds - defined based on the independent calibration via wiggins scores
fathmxf_upper_threshold = 0.14
fathmxf_lower_threshold = 0.12
#eigen - score designated via eigen_Eigen.raw
#thresholds - defined based on the independent calibration via wiggins scores
eigen_upper_threshold = 0.594
eigen_lower_threshold = 0.394
#thresholds - defined based on the independent calibration via wiggins scores
REMM_upper_threshold = 0.86
REMM_lower_threshold = 0.80

####################basic reference set characteristics############
head(cisregANNOT_c)
cisregANNOT_c$Group = as.factor(ifelse(cisregANNOT_c$INFO == "Cont", "Control", "Disease"))

sum(is.na(cisregANNOT_c$CADD_PHRED))
sum(is.na(cisregANNOT_c$CADD_scorebin))

sum(is.na(cisregANNOT_c$REMM_score))
sum(is.na(cisregANNOT_c$REMM_scorebin))

#How many variants in test set scored
#Control n
cont_n = sum(cisregANNOT_c$INFO=="Cont")

#disease n
dis_n = sum(cisregANNOT_c$INFO=="Dis")

#test n
all_n = nrow(cisregANNOT_c)

###CADD analysis#############################################################
cisregANNOT_c$CADD_PHRED = as.numeric(cisregANNOT_c$CADD_PHRED)
cisregANNOT_c %>% group_by(CADD_scorebin) %>% summary()

#contingency table of exp.group vs the score bin 
table(cisregANNOT_c$INFO,cisregANNOT_c$CADD_scorebin)

#create a summary stats
#summary stats of the CADD Phred by group
cisregANNOT_c %>% group_by(INFO) %>% 
  summarize(min = min(CADD_PHRED),
            q1 = quantile(CADD_PHRED, 0.25),
            median = median(CADD_PHRED),
            mean = mean(CADD_PHRED),
            q3 = quantile(CADD_PHRED, 0.75),
            max = max(CADD_PHRED))

#graphing CADD score control versus disease PHRED - WITH PEJAVER MISSENSE THRESHOLDS
cadd = ggplot(cisregANNOT_c, aes(INFO, CADD_PHRED))
cadd + geom_boxplot() + geom_hline(yintercept = 25.3) + geom_hline(yintercept = 22.7)

#graphing test score control versus disease  - WITH THRESHOLD 8 and 10 indciated
ggplot(cisregANNOT_c, aes(CADD_PHRED, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  geom_vline(xintercept = CADD_lower_threshold) + geom_vline(xintercept = CADD_upper_threshold)+ 
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

#REMM_score
summary(cisregANNOT_c$CADD_score)
cisregANNOT_c$CADD_score = as.numeric(cisregANNOT_c$CADD_score)
cisregANNOT_c %>% group_by(CADD_scorebin) %>% summary()

#contingency table of exp.group vs the score bin 
table(cisregANNOT_c$INFO,cisregANNOT_c$CADD_scorebin)

######CADD LRs###########################
#extract a table of REMM_score LR, LR CI low and LR CIhigh
#creating new groups
colnames(cisregANNOT_c)

#need to convert this to a table
summary(as.factor(cisregANNOT_c$CADD_score))
cisregANNOT_c %>% select("Group","CADD_score") %>% table()

#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
#test
test = "CADD"
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$CADD_scorebin), ] %>% select(., "Group","CADD_scorebin")

table(df_test2)
summary(df_test2)
#Control n in test set (all variants)
sum(df_test2$Group=="Control")
#disease n in test set (all variants)
sum(df_test2$Group=="Disease")

#How many variants in test set scored
#Controlrol n
test2_cont_n = sum(df_test2$Group=="Control")

#disease n
test2_dis_n = sum(df_test2$Group=="Disease")

#test n
test2_all_n = nrow(df_test2)

#test1-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(Group=="Disease" & CADD_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(Group=="Control" & CADD_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(Group=="Control" & CADD_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(Group=="Disease" & CADD_scorebin == "score_low") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(Group=="Disease" & CADD_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(Group=="Control" & CADD_scorebin == "score_no_call") %>% count()

# calculate sensitivity
test2_sens = test2_TP/(test2_TP + test2_FN)
print(test2_sens)

# cal specificity
test2_spec = test2_TN/(test2_TN + test2_FP)
print(test2_spec)

# calculate the LR for pos
test2_LRpos = (test2_TP/test2_dis_n)/(test2_FP/test2_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR
test2_varLRpos = (1/test2_TP-1/test2_cont_n)+(1/test2_FP-1/test2_dis_n)
print(test2_varLRpos)

#calc low CI for the LR-
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR-
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

# calculate the LR for negative
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

#create 1st value with test type
#create 2nd value with test group ()
#next three columns as LR, var, low CI and high CI
CADD_pos = c(test = test,
             group = "LRpos",
             LR = as.numeric(test2_LRpos), 
             Var = as.numeric(test2_varLRpos), 
             CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
             CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

CADD_int = c(test = test,
             group = "LRuninf",
             LR = as.numeric(test2_LRnocall), 
             Var = as.numeric(test2_varLRnocall), 
             CI_low = as.numeric(exp(log(test2_LRnocall) - 1.96 * sqrt(test2_varLRnocall))), 
             CI_high = as.numeric(exp(log(test2_LRnocall) + 1.96 * sqrt(test2_varLRnocall))))

CADD_neg = c(test = test,
            group = "LRneg",
            LR = as.numeric(test2_LRneg), 
            Var = as.numeric(test2_varLRneg), 
            CI_low = as.numeric(exp(log(test2_LRneg) - 1.96 * sqrt(test2_varLRneg))), 
            CI_high = as.numeric(exp(log(test2_LRneg) + 1.96 * sqrt(test2_varLRneg))))

CADD_results = type.convert(as.data.frame(rbind(CADD_pos,CADD_int,CADD_neg), row.names = FALSE), as.is = TRUE)
CADD_results

###REMM analysis#############################################################
colnames(cisregANNOT_c)
#REMM_score
summary(cisregANNOT_c$REMM_score)
cisregANNOT_c$REMM_score = as.numeric(cisregANNOT_c$REMM_score)
cisregANNOT_c %>% group_by(REMM_scorebin) %>% summary()

#contingency table of exp.group vs the score bin 
table(cisregANNOT_c$INFO,cisregANNOT_c$REMM_scorebin)

#create a summary stats
#summary stats of the CADD Phred by group
cisregANNOT_c %>% group_by(INFO) %>% 
  summarize(min = min(REMM_score),
            q1 = quantile(REMM_score, 0.25),
            median = median(REMM_score),
            mean = mean(REMM_score),
            q3 = quantile(REMM_score, 0.75),
            max = max(REMM_score))

#graphing CADD score control versus disease PHRED
ggplot(cisregANNOT_c, aes(INFO, REMM_score))+
  geom_boxplot() + geom_hline(yintercept = 0.86) + geom_hline(yintercept = 0.8)

#graphing test score control versus disease  - WITH THRESHOLD 8 and 10 indciated
ggplot(cisregANNOT_c, aes(REMM_score, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  geom_vline(xintercept = REMM_lower_threshold) + geom_vline(xintercept = REMM_upper_threshold)+ 
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

######REMM LRs###############################################
#extract a table of REMM_score LR, LR CI low and LR CIhigh
#creating new groups
colnames(cisregANNOT_c)

#need to convert this to a table
summary(as.factor(cisregANNOT_c$REMM_score))
cisregANNOT_c %>% select("Group","REMM_score") %>% table()

#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
#test
test = "REMM"
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$REMM_scorebin), ] %>% select(., "Group","REMM_scorebin")

table(df_test2)
summary(df_test2)
#Control n in test set (all variants)
sum(df_test2$Group=="Control")
#disease n in test set (all variants)
sum(df_test2$Group=="Disease")

#How many variants in test set scored
#Controlrol n
test2_cont_n = sum(df_test2$Group=="Control")

#disease n
test2_dis_n = sum(df_test2$Group=="Disease")

#test n
test2_all_n = nrow(df_test2)

#test1-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(Group=="Disease" & REMM_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(Group=="Control" & REMM_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(Group=="Control" & REMM_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(Group=="Disease" & REMM_scorebin == "score_low") %>% count()
print(test2_FN)


#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(Group=="Disease" & REMM_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(Group=="Control" & REMM_scorebin == "score_no_call") %>% count()

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
# LR method used in brnich - prop Control variants pred path/prop disease variants predicted pathogenic
test2_LRpos = (test2_TP/test2_dis_n)/(test2_FP/test2_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR -
test2_varLRpos = (1/test2_TP-1/test2_cont_n)+(1/test2_FP-1/test2_dis_n)
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

#create 1st value with test type
#create 2nd value with test group ()
#next three columns as LR, var, low CI and high CI

REMM_pos = c(test = test,
             group = "LRpos",
             LR = as.numeric(test2_LRpos), 
             Var = as.numeric(test2_varLRpos), 
             CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
             CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

REMM_int = c(test = test,
             group = "LRuninf",
             LR = as.numeric(test2_LRnocall), 
             Var = as.numeric(test2_varLRnocall), 
             CI_low = as.numeric(exp(log(test2_LRnocall) - 1.96 * sqrt(test2_varLRnocall))), 
             CI_high = as.numeric(exp(log(test2_LRnocall) + 1.96 * sqrt(test2_varLRnocall))))

REMM_neg = c(test = test,
            group = "LRneg",
            LR = as.numeric(test2_LRneg), 
            Var = as.numeric(test2_varLRneg), 
            CI_low = as.numeric(exp(log(test2_LRneg) - 1.96 * sqrt(test2_varLRneg))), 
            CI_high = as.numeric(exp(log(test2_LRneg) + 1.96 * sqrt(test2_varLRneg))))

REMM_results = type.convert(as.data.frame(rbind(REMM_pos,REMM_int,REMM_neg), row.names = FALSE), as.is = TRUE)
REMM_results

###FATHMM-MKL analysis#############################################################
colnames(cisregANNOT_c)
#fathmmmkl
summary(cisregANNOT_c$fthmmkl_Non.Coding.Score_grc37)
cisregANNOT_c$fthmmkl_Non.Coding.Score_grc37 = as.numeric(cisregANNOT_c$fthmmkl_Non.Coding.Score_grc37)
cisregANNOT_c %>% group_by(fathmmlk_scorebin) %>% summary()

#contingency table of exp.group vs the score bin 
table(cisregANNOT_c$INFO,cisregANNOT_c$fathmmlk_scorebin)

#create a summary stats
#summary stats of the score by group
cisregANNOT_c %>% group_by(INFO) %>%
  summarize(min = min(fthmmkl_Non.Coding.Score_grc37, na.rm= TRUE),
            q1 = quantile(fthmmkl_Non.Coding.Score_grc37, 0.25, na.rm= TRUE),
            median = median(fthmmkl_Non.Coding.Score_grc37, na.rm= TRUE),
            mean = mean(fthmmkl_Non.Coding.Score_grc37, na.rm= TRUE),
            q3 = quantile(fthmmkl_Non.Coding.Score_grc37, 0.75, na.rm= TRUE),
            max = max(fthmmkl_Non.Coding.Score_grc37, na.rm= TRUE))

#graphing score density plot - control versus disease 
ggplot(cisregANNOT_c, aes(INFO, fthmmkl_Non.Coding.Score_grc37))+
  geom_boxplot() + geom_hline(yintercept = 0.86) + geom_hline(yintercept = 0.8)

#graphing test score control versus disease  - WITH THRESHOLD  indciated
ggplot(cisregANNOT_c, aes(fthmmkl_Non.Coding.Score_grc37, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  geom_vline(xintercept = fathmmlk_lower_threshold) + geom_vline(xintercept = fathmmlk_upper_threshold)+ 
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

######FATHMM-MKL LRs#############################################################
##test
test = "FATHMMMKL"
sum(is.na(cisregANNOT_c$fthmmkl_Non.Coding.Score_grc37))
sum(is.na(cisregANNOT_c$fathmmlk_scorebin))
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$fathmmlk_scorebin), ] %>% select(., "Group","fathmmlk_scorebin")

table(df_test2)
summary(df_test2)
#Control n in test set (all variants)
sum(df_test2$Group=="Control")
#disease n in test set (all variants)
sum(df_test2$Group=="Disease")

#How many variants in test set scored
#Controlrol n
test2_cont_n = sum(df_test2$Group=="Control")

#disease n
test2_dis_n = sum(df_test2$Group=="Disease")

#test n
test2_all_n = nrow(df_test2)

#test1-evaluation numbers
#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(Group=="Disease" & fathmmlk_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(Group=="Control" & fathmmlk_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(Group=="Control" & fathmmlk_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(Group=="Disease" & fathmmlk_scorebin == "score_low") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(Group=="Disease" & fathmmlk_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(Group=="Control" & fathmmlk_scorebin == "score_no_call") %>% count()

# calculate sensitivity
test2_sens = test2_TP/(test2_TP + test2_FN)
print(test2_sens)

# cal specificity
test2_spec = test2_TN/(test2_TN + test2_FP)
print(test2_spec)

# calculate the LR for pos +
test2_LRpos = (test2_TP/test2_dis_n)/(test2_FP/test2_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR -
test2_varLRpos = (1/test2_TP-1/test2_cont_n)+(1/test2_FP-1/test2_dis_n)
print(test2_varLRpos)

#calc low CI for the LR-
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR-
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

# calculate the LR for negative -
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

#create 1st value with test type
#create 2nd value with test group ()
#next three columns as LR, var, low CI and high CI

FATHMMMKL_pos = c(test = test,
             group = "LRpos",
             LR = as.numeric(test2_LRpos), 
             Var = as.numeric(test2_varLRpos), 
             CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
             CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

FATHMMMKL_int = c(test = test,
             group = "LRuninf",
             LR = as.numeric(test2_LRnocall), 
             Var = as.numeric(test2_varLRnocall), 
             CI_low = as.numeric(exp(log(test2_LRnocall) - 1.96 * sqrt(test2_varLRnocall))), 
             CI_high = as.numeric(exp(log(test2_LRnocall) + 1.96 * sqrt(test2_varLRnocall))))

FATHMMMKL_neg = c(test = test,
            group = "LRneg",
            LR = as.numeric(test2_LRneg), 
            Var = as.numeric(test2_varLRneg), 
            CI_low = as.numeric(exp(log(test2_LRneg) - 1.96 * sqrt(test2_varLRneg))), 
            CI_high = as.numeric(exp(log(test2_LRneg) + 1.96 * sqrt(test2_varLRneg))))

FATHMMMKL_results = type.convert(as.data.frame(rbind(FATHMMMKL_pos,FATHMMMKL_int,FATHMMMKL_neg), row.names = FALSE), as.is = TRUE)
FATHMMMKL_results

###FATHMM XF analysis#############################################################
colnames(cisregANNOT_c)
#fthmxf_Non.Coding.Score_grc38
#fathmxf_scorebin

#fathmmmkl
summary(cisregANNOT_c$fthmxf_Non.Coding.Score_grc38)
cisregANNOT_c$fthmxf_Non.Coding.Score_grc38 = as.numeric(cisregANNOT_c$fthmxf_Non.Coding.Score_grc38)
cisregANNOT_c %>% group_by(fathmxf_scorebin) %>% summary()

#contingency table of exp.group vs the score bin 
table(cisregANNOT_c$INFO,cisregANNOT_c$fathmxf_scorebin)

#create a summary stats
#summary stats of the score by group
cisregANNOT_c %>% group_by(INFO) %>% 
  summarize(min = min(fthmxf_Non.Coding.Score_grc38, na.rm= TRUE),
            q1 = quantile(fthmxf_Non.Coding.Score_grc38, 0.25, na.rm= TRUE),
            median = median(fthmxf_Non.Coding.Score_grc38, na.rm= TRUE),
            mean = mean(fthmxf_Non.Coding.Score_grc38, na.rm= TRUE),
            q3 = quantile(fthmxf_Non.Coding.Score_grc38, 0.75, na.rm= TRUE),
            max = max(fthmxf_Non.Coding.Score_grc38, na.rm= TRUE))

#graphing score density plot - control versus disease 
ggplot(cisregANNOT_c, aes(INFO, fthmxf_Non.Coding.Score_grc38))+
  geom_boxplot() + geom_hline(yintercept = 0.86) + geom_hline(yintercept = 0.8)

#graphing test score control versus disease  - WITH THRESHOLD  indciated
ggplot(cisregANNOT_c, aes(fthmxf_Non.Coding.Score_grc38, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  geom_vline(xintercept = fathmmlk_lower_threshold) + geom_vline(xintercept = fathmmlk_upper_threshold)+ 
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

#####FATHMM XF LRs#############################################################
test = "FATHMMXF"
sum(is.na(cisregANNOT_c$fthmxf_Non.Coding.Score_grc38))
sum(is.na(cisregANNOT_c$fathmxf_scorebin))
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$fathmxf_scorebin), ] %>% select(., "Group","fathmxf_scorebin")

table(df_test2)
summary(df_test2)
#Control n in test set (all variants)
sum(df_test2$Group=="Control")
#disease n in test set (all variants)
sum(df_test2$Group=="Disease")

#How many variants in test set scored
#Controlrol n
test2_cont_n = sum(df_test2$Group=="Control")

#disease n
test2_dis_n = sum(df_test2$Group=="Disease")

#test n
test2_all_n = nrow(df_test2)

#test1-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(Group=="Disease" & fathmxf_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(Group=="Control" & fathmxf_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(Group=="Control" & fathmxf_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(Group=="Disease" & fathmxf_scorebin == "score_low") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(Group=="Disease" & fathmxf_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(Group=="Control" & fathmxf_scorebin == "score_no_call") %>% count()

# calculate sensitivity
test2_sens = test2_TP/(test2_TP + test2_FN)
print(test2_sens)

# cal specificity
test2_spec = test2_TN/(test2_TN + test2_FP)
print(test2_spec)

# calculate the LR for pos +
test2_LRpos = (test2_TP/test2_dis_n)/(test2_FP/test2_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR+
#var(logLR) for the LR 
test2_varLRpos = (1/test2_TP-1/test2_cont_n)+(1/test2_FP-1/test2_dis_n)
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

#create 1st value with test type
#create 2nd value with test group ()
#next three columns as LR, var, low CI and high CI

FATHMMXF_pos = c(test = test,
                  group = "LRpos",
                  LR = as.numeric(test2_LRpos), 
                  Var = as.numeric(test2_varLRpos), 
                  CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
                  CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

FATHMMXF_int = c(test = test,
                  group = "LRuninf",
                  LR = as.numeric(test2_LRnocall), 
                  Var = as.numeric(test2_varLRnocall), 
                  CI_low = as.numeric(exp(log(test2_LRnocall) - 1.96 * sqrt(test2_varLRnocall))), 
                  CI_high = as.numeric(exp(log(test2_LRnocall) + 1.96 * sqrt(test2_varLRnocall))))

FATHMMXF_neg = c(test = test,
                  group = "LRneg",
                  LR = as.numeric(test2_LRneg), 
                  Var = as.numeric(test2_varLRneg), 
                  CI_low = as.numeric(exp(log(test2_LRneg) - 1.96 * sqrt(test2_varLRneg))), 
                  CI_high = as.numeric(exp(log(test2_LRneg) + 1.96 * sqrt(test2_varLRneg))))

FATHMMXF_results = type.convert(as.data.frame(rbind(FATHMMXF_pos,FATHMMXF_int,FATHMMXF_neg), row.names = FALSE), as.is = TRUE)
FATHMMXF_results

###Eigen analysis#############################################################
colnames(cisregANNOT_c)

#eigen
summary(cisregANNOT_c$eigen_Eigen.PC.raw)
cisregANNOT_c$eigen_Eigen.PC.raw = as.numeric(cisregANNOT_c$eigen_Eigen.PC.raw)
cisregANNOT_c %>% group_by(eigen_scorebin) %>% summary()

#contingency table of exp.group vs the score bin 
table(cisregANNOT_c$INFO,cisregANNOT_c$eigen_scorebin)

#create a summary stats
#summary stats of the score by group
cisregANNOT_c %>% group_by(INFO) %>% 
  summarize(min = min(eigen_Eigen.PC.raw, na.rm = TRUE),
            q1 = quantile(eigen_Eigen.PC.raw, 0.25, na.rm = TRUE),
            median = median(eigen_Eigen.PC.raw, na.rm = TRUE),
            mean = mean(eigen_Eigen.PC.raw, na.rm = TRUE),
            q3 = quantile(eigen_Eigen.PC.raw, 0.75, na.rm = TRUE),
            max = max(eigen_Eigen.PC.raw, na.rm = TRUE))

#graphing score density plot - control versus disease 
ggplot(cisregANNOT_c, aes(INFO, eigen_Eigen.PC.raw))+
  geom_boxplot() + geom_hline(yintercept = eigen_lower_threshold) + geom_hline(yintercept = eigen_upper_threshold)

#graphing test score control versus disease  - WITH THRESHOLD  indciated
ggplot(cisregANNOT_c, aes(eigen_Eigen.PC.raw, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  geom_vline(xintercept = eigen_lower_threshold) + geom_vline(xintercept = eigen_upper_threshold)+ 
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

###Eigen LRs#############################################################
#determine LR for each category, plus CI low and LR CIhigh
test = "Eigen"
sum(is.na(cisregANNOT_c$eigen_Eigen.PC.raw))
sum(is.na(cisregANNOT_c$eigen_scorebin))
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$eigen_scorebin), ] %>% select(., "Group","eigen_scorebin")

table(df_test2)
summary(df_test2)
#Control n in test set (all variants)
sum(df_test2$Group=="Control")
#disease n in test set (all variants)
sum(df_test2$Group=="Disease")

#How many variants in test set scored
#Controlrol n
test2_cont_n = sum(df_test2$Group=="Control")

#disease n
test2_dis_n = sum(df_test2$Group=="Disease")

#test n
test2_all_n = nrow(df_test2)

#test1-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(Group=="Disease" & eigen_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(Group=="Control" & eigen_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(Group=="Control" & eigen_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(Group=="Disease" & eigen_scorebin == "score_low") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(Group=="Disease" & eigen_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(Group=="Control" & eigen_scorebin == "score_no_call") %>% count()

# calculate sensitivity
test2_sens = test2_TP/(test2_TP + test2_FN)
print(test2_sens)

# cal specificity
test2_spec = test2_TN/(test2_TN + test2_FP)
print(test2_spec)

# calculate the LR for pos +
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

# calculate the LR for negative -
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
test2_LRneg = (test2_FN/test2_dis_n)/(test2_TN/test2_cont_n)
print(test2_LRneg)

#std error of LR -
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

#create 1st value with test type
#create 2nd value with test group ()
#next three columns as LR, var, low CI and high CI

Eigen_pos = c(test = test,
                 group = "LRpos",
                 LR = as.numeric(test2_LRpos), 
                 Var = as.numeric(test2_varLRpos), 
                 CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
                 CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

Eigen_int = c(test = test,
                 group = "LRuninf",
                 LR = as.numeric(test2_LRnocall), 
                 Var = as.numeric(test2_varLRnocall), 
                 CI_low = as.numeric(exp(log(test2_LRnocall) - 1.96 * sqrt(test2_varLRnocall))), 
                 CI_high = as.numeric(exp(log(test2_LRnocall) + 1.96 * sqrt(test2_varLRnocall))))

Eigen_neg = c(test = test,
                 group = "LRneg",
                 LR = as.numeric(test2_LRneg), 
                 Var = as.numeric(test2_varLRneg), 
                 CI_low = as.numeric(exp(log(test2_LRneg) - 1.96 * sqrt(test2_varLRneg))), 
                 CI_high = as.numeric(exp(log(test2_LRneg) + 1.96 * sqrt(test2_varLRneg))))

Eigen_results = type.convert(as.data.frame(rbind(Eigen_pos,Eigen_int,Eigen_neg), row.names = FALSE), as.is = TRUE)
Eigen_results

##Linsight analysis#############################################################
colnames(cisregANNOT_c)
#fathmmmkl
summary(cisregANNOT_c$linsight_Lscore_hg19)
cisregANNOT_c$linsight_Lscore_hg19 = as.numeric(cisregANNOT_c$linsight_Lscore_hg19)
cisregANNOT_c %>% group_by(linsight_scorebin) %>% summary()

#contingency table of exp.group vs the score bin 
table(cisregANNOT_c$INFO,cisregANNOT_c$linsight_scorebin)

#create a summary stats
#summary stats of the score by group
cisregANNOT_c %>% group_by(INFO) %>% 
  summarize(min = min(linsight_Lscore_hg19, na.rm = TRUE),
            q1 = quantile(linsight_Lscore_hg19, 0.25, na.rm = TRUE),
            median = median(linsight_Lscore_hg19, na.rm = TRUE),
            mean = mean(linsight_Lscore_hg19, na.rm = TRUE),
            q3 = quantile(linsight_Lscore_hg19, 0.75, na.rm = TRUE),
            max = max(linsight_Lscore_hg19),na.rm = TRUE)

#graphing score density plot - control versus disease 
ggplot(cisregANNOT_c, aes(INFO, linsight_Lscore_hg19))+
  geom_boxplot() + geom_hline(yintercept = linsight_lower_threshold) + geom_hline(yintercept = linsight_upper_threshold)

#graphing test score control versus disease  - WITH THRESHOLD  indicated
ggplot(cisregANNOT_c, aes(linsight_Lscore_hg19, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  geom_vline(xintercept = linsight_lower_threshold) + geom_vline(xintercept = linsight_upper_threshold)+ 
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

##Linsight LRs#############################################################
#determine LR for each category, plus CI low and LR CIhigh
test = "LINSIGHT"
sum(is.na(cisregANNOT_c$linsight_Lscore_hg19))
sum(is.na(cisregANNOT_c$linsight_scorebin))
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$linsight_scorebin), ] %>% select(., "Group","linsight_scorebin")

table(df_test2)
summary(df_test2)
#Control n in test set (all variants)
sum(df_test2$Group=="Control")
#disease n in test set (all variants)
sum(df_test2$Group=="Disease")

#How many variants in test set scored
#Controlrol n
test2_cont_n = sum(df_test2$Group=="Control")

#disease n
test2_dis_n = sum(df_test2$Group=="Disease")

#test n
test2_all_n = nrow(df_test2)

#test1-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(Group=="Disease" & linsight_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(Group=="Control" & linsight_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(Group=="Control" & linsight_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(Group=="Disease" & linsight_scorebin == "score_low") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(Group=="Disease" & linsight_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(Group=="Control" & linsight_scorebin == "score_no_call") %>% count()

# calculate sensitivity
test2_sens = test2_TP/(test2_TP + test2_FN)
print(test2_sens)

# cal specificity
test2_spec = test2_TN/(test2_TN + test2_FP)
print(test2_spec)

# calculate the LR for pos +
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

# calculate the LR for negative -
## LR method used in brnich - prop cont variants pred path/prop disease variants predicted pathogenic
test2_LRneg = (test2_FN/test2_dis_n)/(test2_TN/test2_cont_n)
print(test2_LRneg)

#std error of LR-
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

#create 1st value with test type
#create 2nd value with test group ()
#next three columns as LR, var, low CI and high CI

LINSIGHT_pos = c(test = test,
              group = "LRpos",
              LR = as.numeric(test2_LRpos), 
              Var = as.numeric(test2_varLRpos), 
              CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
              CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

LINSIGHT_int = c(test = test,
              group = "LRuninf",
              LR = as.numeric(test2_LRnocall), 
              Var = as.numeric(test2_varLRnocall), 
              CI_low = as.numeric(exp(log(test2_LRnocall) - 1.96 * sqrt(test2_varLRnocall))), 
              CI_high = as.numeric(exp(log(test2_LRnocall) + 1.96 * sqrt(test2_varLRnocall))))

LINSIGHT_neg = c(test = test,
              group = "LRneg",
              LR = as.numeric(test2_LRneg), 
              Var = as.numeric(test2_varLRneg), 
              CI_low = as.numeric(exp(log(test2_LRneg) - 1.96 * sqrt(test2_varLRneg))), 
              CI_high = as.numeric(exp(log(test2_LRneg) + 1.96 * sqrt(test2_varLRneg))))

LINSIGHT_results = type.convert(as.data.frame(rbind(LINSIGHT_pos,LINSIGHT_int,LINSIGHT_neg), row.names = FALSE), as.is = TRUE)
LINSIGHT_results


###############################################################
#reg_variant_tool_evaluation - LRs
#thresholds - defined based on the independent calibration 
#graph the LRs of each group

#data sourced from Supplemental table 10_Combined_score_evaluation
#cisreg_eval = read.table("Data/reg_variant_tool_evaluation - summary.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
#colnames(cisreg_eval)

#cisregLRs = read.table("Data/reg_variant_tool_evaluation - LRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
#colnames(cisregLRs)

#join to create a combined table of all tool results
cisregLRs = rbind(CADD_results,REMM_results, FATHMMMKL_results, FATHMMXF_results, Eigen_results, LINSIGHT_results)
colnames(cisregLRs) = c("annotation","LR_category","LR","Var","CI_low","CI_high")

write.table(cisregLRs,(paste("Data/",date,"_reg_variant_tool_evaluation_LRs.txt")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Figure 4 - Bioinformatic tool calibration results

#CADD

cadd_density = ggplot(cisregANNOT_c, aes(CADD_PHRED, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  labs(x = "CADD", y = "density")+
  geom_vline(xintercept = CADD_lower_threshold) + geom_vline(xintercept = CADD_upper_threshold)+ 
  theme(legend.title = element_blank(), legend.position = "none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

cadd_LRplot = ggplot(data = (cisregLRs[cisregLRs$annotation == "CADD", ]), aes(x=reorder(LR_category,LR), y = LR, ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "", y = "Likelihood ratio log10 (95% CI)")+
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

#remm

remm_density = ggplot(cisregANNOT_c, aes(REMM_score, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  labs(x = "REMM", y = "density")+
  geom_vline(xintercept = REMM_lower_threshold) + geom_vline(xintercept = REMM_upper_threshold)+ 
  theme(legend.title = element_blank(),legend.position = "none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

remm_LRplot = ggplot(data = (cisregLRs[cisregLRs$annotation == "REMM", ]), aes(x=reorder(LR_category,LR), y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "", y = "Likelihood ratio log10 (95% CI)")+
  geom_hline(yintercept = 0.23,  colour = "darkgreen") +
  annotation_custom(grobTree(textGrob("moderate", x=0.08,  y=0.05, hjust=0,gp=gpar(col="darkgreen", fontsize=8))))+
  geom_hline(yintercept = 0.48, colour = "green") +
  annotation_custom(grobTree(textGrob("supporting", x=0.20,  y=0.05, hjust=0,gp=gpar(col="chartreuse3", fontsize=8))))+
  geom_hline(yintercept = 1, colour = "gold") +
  geom_hline(yintercept = 2.1, colour = "orange") +
  annotation_custom(grobTree(textGrob("supporting", x=0.57,  y=0.05, hjust=0,gp=gpar(col="orange", fontsize=8))))+
  geom_hline(yintercept = 4.3, colour = "red")+
  annotation_custom(grobTree(textGrob("moderate", x=0.7,  y=0.05, hjust=0,gp=gpar(col="red", fontsize=8))))+  geom_hline(yintercept = 18.7, colour = "brown4")+
  annotation_custom(grobTree(textGrob("strong", x=0.94,  y=0.05, hjust=0,gp=gpar(col="brown4", fontsize=8))))

#fathmmkl

fathmmkl_density = ggplot(cisregANNOT_c, aes(fthmmkl_Non.Coding.Score_grc37, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  labs(x = "FATHMM-MKL", y = "density")+
  geom_vline(xintercept = fathmmlk_lower_threshold) + geom_vline(xintercept = fathmmlk_upper_threshold)+ 
  theme(legend.title = element_blank(),legend.position = "none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

fathmmkl_LRplot = ggplot(data = (cisregLRs[cisregLRs$annotation == "FATHMMMKL", ]), aes(x=reorder(LR_category,LR), y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "", y = "Likelihood ratio log10 (95% CI)")+
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

#fathmmxf
fathmmxf_density = ggplot(cisregANNOT_c, aes(fthmxf_Non.Coding.Score_grc38, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  labs(x = "FATHMM-XF", y = "density")+
  geom_vline(xintercept = fathmxf_lower_threshold) + geom_vline(xintercept = fathmxf_upper_threshold)+ 
  theme(legend.title = element_blank(),legend.position = "none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

fathmmxf_LRplot = ggplot(data = (cisregLRs[cisregLRs$annotation == "FATHMMXF", ]), aes(x=reorder(LR_category,LR), y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "", y = "Likelihood ratio log10 (95% CI)")+
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

#eigen
eigen_density = ggplot(cisregANNOT_c, aes(eigen_Eigen.PC.raw, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  labs(x = "Eigen", y = "density")+
  geom_vline(xintercept = eigen_lower_threshold) + geom_vline(xintercept = eigen_upper_threshold)+ 
  theme(legend.title = element_blank(),legend.position = "none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

eigen_LRplot = ggplot(data = (cisregLRs[cisregLRs$annotation == "Eigen", ]), aes(x=reorder(LR_category,LR), y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "", y = "Likelihood ratio log10 (95% CI)")+
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

#linsight

linsight_density = ggplot(cisregANNOT_c, aes(linsight_Lscore_hg19, fill = INFO))+
  geom_density(position = "identity", alpha = 0.2)+
  labs(x = "LINSIGHT", y = "density")+
  geom_vline(xintercept = linsight_lower_threshold) + geom_vline(xintercept = linsight_upper_threshold)+ 
  theme(legend.title = element_blank(),legend.position = "none")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999")

linsight_LRplot = ggplot(data = (cisregLRs[cisregLRs$annotation == "LINSIGHT", ]), aes(x=reorder(LR_category,LR), y = LR,ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
  coord_flip() + scale_y_continuous(trans = "log10", limits = c(0.1,23))+
  labs(x = "", y = "Likelihood ratio log10 (95% CI)")+
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

#Figure 4 - summary of all tool scores and calibrated LRs

legend = get_legend(ggplot(cisregANNOT_c, aes(CADD_PHRED, fill = Group))+
                      geom_density(position = "identity", alpha = 0.2)+
                      labs(x = "CADD", y = "density")+
                      theme(legend.position = "bottom", legend.title = element_blank())+
                      geom_vline(xintercept = CADD_lower_threshold) + geom_vline(xintercept = CADD_upper_threshold)+
                      scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999"))

ggdraw(plot_grid((plot_grid(cadd_density, cadd_LRplot,
                 remm_density, remm_LRplot,
                 fathmmkl_density, fathmmkl_LRplot,
                 fathmmxf_density, fathmmxf_LRplot,
                 eigen_density, eigen_LRplot,
                 linsight_density, linsight_LRplot,
                 align = "hv", axis = "b", rel_widths = c(0.75, 1.5),
                 labels = c("A","B", "C", "D","E", "F", "G", "H", "I", "J", "K", "L" ),ncol = 2, nrow = 6)),
                 (plot_grid(legend, align = "v", axis = "b", rel_widths = c(1,1,1), ncol = 3, nrow = 1)),
                 rel_heights = c(6,1),
                 ncol = 1, 
                 nrow = 2))

#ggsave(path = "Data/", filename = "Fig4wlegend.png", width = 10, height = 12, device='png', dpi=1200)
