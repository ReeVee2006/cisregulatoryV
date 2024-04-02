##annotations and categories annotated in Fig5_calibration
##########-checking-methods-to-combine------------------
##calculate LRs

library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(cowplot)
library(grid)

#set up files--------------------------------------------
#setwd 
#combine scores mathematically, see if calibration is better

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

##################setting up LR calculations
# create a summary table for your diagnostic  evaluations
colnames(cisregANNOT_c)
results = cisregANNOT_c %>% select(.,"ID","INFO","Group",
                                   "CADD_RawScore","CADD_PHRED","CADD_scorebin","CADD_testbin",
                                   "REMM_grc38","REMM_score","REMM_scorebin","REMM_testbin")

#basic values for stats calculations
#Control n in test set (all variants)
sum(results$INFO=="Cont")
#disease n in test set (all variants)
sum(results$INFO=="Dis")

#How many variants in test set scored
#Control n
cont_n = sum(results$Group=="Control")

#disease n
dis_n = sum(results$Group=="Disease")

#test n
all_n = nrow(results)


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

#check and summary
summary(as.factor(results$test_comb1))
summary(as.factor(results$CADD_scorebin))
results %>% select("INFO","test_comb1") %>% table()
sum(is.na(results$test_comb1))

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
comb1_cont_n

#disease n
comb1_dis_n = sum(df_comb1$INFO=="Dis")

#test n
comb1_all_n = nrow(df_comb1)

#test1-evaluation 
#determine LR for each category, plus CI low and LR CIhigh
test = "CADD/REMM"
sum(is.na(results$test_comb1))

#calculation for '2high' group---------------------
group = "2high"

#Calculate n True Pos 1 (Event pos - both scores high)
test2_TP1 = results %>% filter(Group=="Disease" & test_comb1 == "2high") %>% count()
print(test2_TP1)

#calculate false positive (no event pos)
test2_FP1 = results %>% filter(Group=="Control" & test_comb1 == "2high") %>% count()
print(test2_FP1)

# calculate the LR for pos +
test2_LRpos = (test2_TP1/comb1_dis_n)/(test2_FP1/comb1_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
test2_varLRpos = (1/test2_TP1-1/comb1_cont_n)+(1/test2_FP1-1/comb1_dis_n)
print(test2_varLRpos)

#calc low CI for the LR-
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR-
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

comb1_2high = c(test = test,
                group = "2high",
                LR = as.numeric(test2_LRpos), 
                Var = as.numeric(test2_varLRpos), 
                CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
                CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

comb1_2high

#calculation for '1high' group---------------------
group = "1high"

#Calculate n True Pos 1 (Event pos - both scores high)
test2_TP1 = results %>% filter(Group=="Disease" & test_comb1 == "1high") %>% count()
print(test2_TP1)

#calculate false positive (no event pos)
test2_FP1 = results %>% filter(Group=="Control" & test_comb1 == "1high") %>% count()
print(test2_FP1)

# calculate the LR for pos +
test2_LRpos = (test2_TP1/comb1_dis_n)/(test2_FP1/comb1_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
test2_varLRpos = (1/test2_TP1-1/comb1_cont_n)+(1/test2_FP1-1/comb1_dis_n)
print(test2_varLRpos)

#calc low CI for the LR-
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR-
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

comb1_1high = c(test = test,
                group = "1high",
                LR = as.numeric(test2_LRpos), 
                Var = as.numeric(test2_varLRpos), 
                CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
                CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

comb1_1high

#calculation for 'conflict' group---------------------
group = "conflict"

#Calculate n True Pos 1 (Event pos - both scores high)
test2_TP1 = results %>% filter(Group=="Disease" & test_comb1 == "conflict") %>% count()
print(test2_TP1)

#calculate false positive (no event pos)
test2_FP1 = results %>% filter(Group=="Control" & test_comb1 == "conflict") %>% count()
print(test2_FP1)

# calculate the LR for pos +
test2_LRpos = (test2_TP1/comb1_dis_n)/(test2_FP1/comb1_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
test2_varLRpos = (1/test2_TP1-1/comb1_cont_n)+(1/test2_FP1-1/comb1_dis_n)
print(test2_varLRpos)

#calc low CI for the LR-
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR-
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

comb1_conflict = c(test = test,
                group = "conflict",
                LR = as.numeric(test2_LRpos), 
                Var = as.numeric(test2_varLRpos), 
                CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
                CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

comb1_conflict

#calculation for 'int' group---------------------
group = "int"

#Calculate n True Pos 1 (Event pos - both scores high)
test2_TP1 = results %>% filter(Group=="Disease" & test_comb1 == "int") %>% count()
print(test2_TP1)

#calculate false positive (no event pos)
test2_FP1 = results %>% filter(Group=="Control" & test_comb1 == "int") %>% count()
print(test2_FP1)

# calculate the LR for pos +
test2_LRpos = (test2_TP1/comb1_dis_n)/(test2_FP1/comb1_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
test2_varLRpos = (1/test2_TP1-1/comb1_cont_n)+(1/test2_FP1-1/comb1_dis_n)
print(test2_varLRpos)

#calc low CI for the LR-
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR-
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

comb1_int = c(test = test,
                   group = "int",
                   LR = as.numeric(test2_LRpos), 
                   Var = as.numeric(test2_varLRpos), 
                   CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
                   CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

comb1_int

#calculation for '1low' group---------------------
group = "1low"

#Calculate n True Pos 1 (Event pos - both scores high)
test2_TP1 = results %>% filter(Group=="Disease" & test_comb1 == "1low") %>% count()
print(test2_TP1)

#calculate false positive (no event pos)
test2_FP1 = results %>% filter(Group=="Control" & test_comb1 == "1low") %>% count()
print(test2_FP1)

# calculate the LR for pos +
test2_LRpos = (test2_TP1/comb1_dis_n)/(test2_FP1/comb1_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
test2_varLRpos = (1/test2_TP1-1/comb1_cont_n)+(1/test2_FP1-1/comb1_dis_n)
print(test2_varLRpos)

#calc low CI for the LR-
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR-
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

comb1_1low = c(test = test,
              group = "1low",
              LR = as.numeric(test2_LRpos), 
              Var = as.numeric(test2_varLRpos), 
              CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
              CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

comb1_1low


#calculation for '2low' group---------------------
group = "2low"

#Calculate n True Pos 1 (Event pos - both scores high)
test2_TP1 = results %>% filter(Group=="Disease" & test_comb1 == "2low") %>% count()
print(test2_TP1)

#calculate false positive (no event pos)
test2_FP1 = results %>% filter(Group=="Control" & test_comb1 == "2low") %>% count()
print(test2_FP1)

# calculate the LR for pos +
test2_LRpos = (test2_TP1/comb1_dis_n)/(test2_FP1/comb1_cont_n)
print(test2_LRpos)

#calc 95% CI for the LR
#var(logLR) for the LR 
test2_varLRpos = (1/test2_TP1-1/comb1_cont_n)+(1/test2_FP1-1/comb1_dis_n)
print(test2_varLRpos)

#calc low CI for the LR-
exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))

#calc high CI for the LR-
exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))

comb1_2low = c(test = test,
               group = "2low",
               LR = as.numeric(test2_LRpos), 
               Var = as.numeric(test2_varLRpos), 
               CI_low = as.numeric(exp(log(test2_LRpos) - 1.96 * sqrt(test2_varLRpos))), 
               CI_high = as.numeric(exp(log(test2_LRpos) + 1.96 * sqrt(test2_varLRpos))))

comb1_2low

######SAVING FULL SET

#join to create a combined table of all tool results
cisregLRs = type.convert(as.data.frame(rbind(comb1_2high,comb1_1high, comb1_conflict, comb1_int, comb1_1low, comb1_2low), row.names = FALSE), as.is = TRUE)

cisregLRs
colnames(cisregLRs) = c("annotation","LR_category","LR","Var","CI_low","CI_high")

write.table(cisregLRs,"Data/tool_evaluation_LRs_CADDREMMcomb.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


##########################

#Figure 5 - LRs of combined CADD REMM
caddremm_LRplot = ggplot(data = (cisregLRs[cisregLRs$annotation == "CADD/REMM", ]), aes(x=reorder(LR_category,LR), y = LR, ymin=CI_low, ymax=CI_high)) + geom_pointrange()+
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

ggsave(path = "D:/cisreg_manuscript/figures", filename = "Fig_5_CADDREMM_LRs.png", width = 10, height = 6, device='png', dpi=1200)

############################

colnames(results)
summary(as.factor(results$test_comb1))

results$test_comb1_bin = as.factor(case_when(
  results$test_comb1 == "2high" | results$test_comb1 == "1high" ~ "pred_path",
  results$test_comb1 == "1low" | results$test_comb1 == "int" | results$test_comb1 == "conflict"  ~ "no_call",
  results$test_comb1 == "2low" ~ "pred_ben"))
summary(results$test_comb1_bin)
sum(is.na(results$test_comb1_bin))

#create 6 bins for evidence categories
results$test_comb1_testbin = as.factor(paste(results$Group,results$test_comb1_bin ,sep = "_"))
summary(results$test_comb1_testbin)
sum(is.na(results$test_comb1_testbin))

CADDREMM_summary = c(tool = "CADD/REMM",
                     all_n = nrow(results),
                     total_n_cont = sum(results$Group=="Control"),
                     total_n_dis = sum(results$Group=="Disease"),
                     n_scored = sum(!is.na(results$test_comb1)),
                     n_scored_cont = sum(results$Group=="Control" & !is.na(results$test_comb1)),
                     n_scored_dis = sum(results$Group=="Disease" & !is.na(results$test_comb1)),
                     n_unscored = sum(is.na(results$test_comb1)),
                     n_uninformative = sum(results$test_comb1_bin == "no_call"),
                     neg_cont_TN = sum(results$test_comb1_testbin == "Control_pred_ben"),
                     int_cont_uninf = sum(results$test_comb1_testbin == "Control_no_call"),
                     pos_cont_FP = sum(results$test_comb1_testbin == "Control_pred_path"),
                     neg_dis_FN = sum(results$test_comb1_testbin == "Disease_pred_ben"),
                     int_dis_uninf = sum(results$test_comb1_testbin == "Disease_no_call"),
                     pos_dis_TP = sum(results$test_comb1_testbin == "Disease_pred_path"),
                     n_correct = sum(results$test_comb1_testbin == "Control_pred_ben")+sum(results$test_comb1_testbin == "Disease_pred_path"),
                     n_incorrect = sum(results$test_comb1_testbin == "Control_pred_path") + sum(results$test_comb1_testbin == "Disease_pred_ben"),
                     n_undetermined = sum(results$test_comb1_testbin == "Control_no_call" | 
                                            results$test_comb1_testbin == "Disease_no_call" |
                                            results$test_comb1_testbin == "Cont_NA" |
                                            results$test_comb1_testbin == "Dis_NA"),
                     percent_sens = (sum(results$test_comb1_testbin == "Disease_pred_path"))/((sum(results$test_comb1_testbin == "Disease_pred_path") + sum(results$test_comb1_testbin == "Disease_pred_ben")))*100,
                     percent_spec = sum(results$test_comb1_testbin == "Control_pred_ben")/(sum(results$test_comb1_testbin == "Control_pred_ben") + sum(results$test_comb1_testbin == "Control_pred_path"))*100,
                     percent_acuracy = (((sum(results$test_comb1_testbin == "Control_pred_ben")+
                                            sum(results$test_comb1_testbin == "Disease_pred_path")))/sum(!is.na(results$test_comb1))*100),
                     percent_unscored = sum(is.na(results$test_comb1))/nrow(results)*100,
                     percent_uninformative = (sum(results$test_comb1_testbin == "Control_no_call" | results$test_comb1_testbin == "Disease_no_call"))/nrow(results)*100,
                     percent_undetermined = (sum(results$test_comb1_testbin == "Control_no_call" | 
                                                   results$test_comb1_testbin == "Disease_no_call" |
                                                   results$test_comb1_testbin == "Cont_NA" |
                                                   results$test_comb1_testbin == "Dis_NA"))/nrow(results)*100,
                     CLIN_percent_sens = (sum(results$test_comb1_testbin == "Disease_pred_path"))/(sum(results$Group=="Disease"))*100,
                     CLIN_percent_spec = sum(results$test_comb1_testbin == "Control_pred_ben")/(sum(results$Group=="Control"))*100,
                     CLIN_percent_accuracy = ((sum(results$test_comb1_testbin == "Control_pred_ben")+
                                                 sum(results$test_comb1_testbin == "Disease_pred_path"))/all_n)*100,
                     CLIN_percent_correct = (sum(results$test_comb1_testbin == "Control_pred_ben")+sum(results$test_comb1_testbin == "Disease_pred_path"))/nrow(results)*100,
                     CLIN_percent_incorrect = (sum(results$test_comb1_testbin == "Control_pred_path") + sum(results$test_comb1_testbin == "Disease_pred_ben"))/nrow(results)*100,
                     CLIN_percent_undetermined = (sum(results$test_comb1_testbin == "Control_no_call" | 
                                                        results$test_comb1_testbin == "Disease_no_call" |
                                                        results$test_comb1_testbin == "Cont_NA" |
                                                        results$test_comb1_testbin == "Dis_NA")/nrow(results)*100))
CADDREMM_summary

summary_CADDREMM = type.convert(as.data.frame(rbind(CADDREMM_summary), row.names = TRUE), as.is = TRUE)

CADDREMM_summary_rounded = CADDREMM_summary %>% dplyr::mutate(across(where(is.numeric), round, 2))
write.table(CADDREMM_summary_rounded,"Data/240322_caddremm_tool_sum_r.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
