##cis reg bioinformatic tool calibration - B - calibration step
##web based input files for annotations
##annotation of variant files
##annotate and translate to categories, calculate LRs
##last executed 26.10.23 RV

library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(cowplot)
library(grid)
#library(liftOver)

#set up files--------------------------------------------
#setwd 
setwd("D:/cisreg_manuscript/SUPP")

## Load datasets-Cisreg variants, control and disease to combine dataset and create into VCF file
#load disease-associated variant set
cisregDraw = read.table("Data/Supplemental_DisVar.txt", header= TRUE, sep = "\t")
#load control variant set
cisregCraw = read.table("Data/Supplemental_ContVar.txt", header= TRUE, sep = "\t")

##Load datasets for annotation###################
##  
#Cisreg variants, control and disease-associated combined dataset and create into VCF file
refset = read.table("temp/cisreg_refset_grc38.txt", header= TRUE, sep = "\t",stringsAsFactors = TRUE)
summary(refset)
colnames(refset)

cisregANNOT = refset[ ,c(1:22)]

#assess the dataset characteristics
#number of SNVs in full dataset
sum(cisregANNOT$allelesize==0)

#number of indels in full dataset
sum(cisregANNOT$allelesize!=0)

#count variant type in each exp.group(INFO)
cisregANNOT %>% count(INFO,allelesize ==0)

summary(as.factor(cisregANNOT$INFO))
refset$INFO %>% as.factor() %>% summary()
refset$ROIgene %>% as.factor() %>% n_distinct()
refset$variant38 %>% as.factor() %>% n_distinct()
refset$ID %>% as.factor() %>% n_distinct()

#
#adding GRCh37/hg19 to the cisregANNOT

#getting the grch37 locations - create a location bed file
refset_bed = refset %>% select(CHROM,POS,REF,ALT)
refset_bed$end = refset_bed$POS+1
refset_bed$bedstart = paste(refset$CHROM,refset$POS,sep = ":")
refset_bed$bedstart = paste("chr",refset_bed$bedstart, sep = "")
refset_bed$bed38 = paste(refset_bed$bedstart, refset_bed$end, sep = "-")
bed38_input = as.data.frame(refset_bed$bed38)

#save grch38 bed file - liftover via https://genome.ucsc.edu/cgi-bin/hgLiftOver
write.table(bed38_input, "temp/bed38_input.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#input grch38 location file above into liftove and collect liftover locations
bed37_input= read.table("Data/hglft_genome_3ad85_b177f0_hg19.bed", sep = "\t", header = FALSE)
colnames(bed37_input) = paste("bed37")

#create a liftover file
refset_liftover = cbind(refset_bed,bed37_input)

#add grch37 locations to the ref set
refset_liftover$location38 = paste(refset_liftover$CHROM, refset_liftover$POS, sep = ":")
refset_liftover$location38ref = paste(refset_liftover$location38, refset_liftover$REF, sep = "")
refset_liftover$variant38 = paste(refset_liftover$location38ref, refset_liftover$ALT, sep = ">")
refset_liftover$bed37b = refset_liftover$bed37
refset_liftover = refset_liftover %>% separate(bed37b, c("chr37","start37","end37"))
refset_liftover$CHROM37 = str_replace(refset_liftover$chr37, "chr","")
refset_liftover$location37 = paste(refset_liftover$CHROM37, refset_liftover$start37,sep = ":")
grc37_cols = refset_liftover %>% select("variant38", "location37","chr37","CHROM37","start37")

cisregANNOT = dplyr::left_join(cisregANNOT, grc37_cols, by = c("variant38" = "variant38"))
cisregANNOT$variant37ref = paste(cisregANNOT$location37, cisregANNOT$REF, sep = "")
cisregANNOT$variant37 = paste(cisregANNOT$variant37ref, cisregANNOT$ALT, sep = ">")

#create GRanges object from cisreg  variant list for input to annotations
#converting dataframe into GRanges object - from vcf like file - grc38 locations
grcisregANNOT =  makeGRangesFromDataFrame(cisregANNOT,
                                          keep.extra.columns=FALSE,
                                          ignore.strand=TRUE,
                                          seqinfo=NULL,
                                          seqnames.field=c("grchr_38"),
                                          start.field=c("POS"),
                                          end.field=c("end.field"),
                                          starts.in.df.are.0based=FALSE)

##VEP annotation#################################################---------------------------------------
##ANNOTATIONS - 2
# VEP annotation - input above file from web GUI, 
##VEP annotation - created a vcf file in FILTER script as input for online GUI VEP see INPUT script, 
##VEP annotation input - cisreg_all38 -  full set of ~13K

##selecting a single annotation transcript database Ensembl
#remove # at start of VEP annotation to read
#VEP annotation performed 11.102023
VEPfull = read.delim("Data/OaVt9eEsWgLnckRd.txt", sep = "\t", header= TRUE)
#formatting VEP to support downstream process
VEPfull$Uploaded_variation %>% as.factor() %>% n_distinct()
colnames(VEPfull)

#count the number of unique variants based on project identified "CR#"
VEPfull$Uploaded_variation %>% as.factor() %>% n_distinct()
VEPfull$MANE_SELECT %>% as.factor() %>% n_distinct()
VEPfull$Gene %>% as.factor() %>% n_distinct()
VEPfull$testsetID = VEPfull$Uploaded_variation

#select the VEP annotations for only the reference set variants
VEP = subset(VEPfull, (Uploaded_variation %in% c(refset$ID)))
VEP$Uploaded_variation %>% as.factor() %>% n_distinct()
VEP$MANE_SELECT %>% as.factor() %>% n_distinct()
VEP$Gene %>% as.factor() %>% n_distinct()
VEP$testsetID = VEP$Uploaded_variation

###selecting the VEP annotation according to the clinically identified (selected) gene---------------------

#ensure that "Selected_gene_name" column in present on the ANNOT file - currently sourced from manual annotation
colnames(cisregANNOT)

#create a table with the selected gene symbol and gene information for each variant
colnames(refset)
var_geneinfo = refset[ ,c("ID","ROIgene","Selected_gene_name")]

# annotate the VEP file with the ROI gene name
VEP = left_join(x = VEP, y= var_geneinfo, by = c("Uploaded_variation" = "ID"))

#select the VEP annotation rows that match the selected gene 
VEP_select = VEP %>% filter(VEP$SYMBOL == VEP$ROIgene)
VEP$SYMBOL %>% as.factor() %>% n_distinct()
VEP_select$SYMBOL %>% as.factor() %>% n_distinct()
VEP_select$Uploaded_variation %>% as.factor() %>% n_distinct()

###selecting single annotation based on MANE or MANE_CLINICAL----------
#select MANEselect transcripts only
VEP_select_m = VEP_select[VEP_select$MANE_SELECT != "-", ]
VEP_select_m_u = distinct(VEP_select_m, Uploaded_variation, .keep_all = TRUE)
VEP_select_m$Uploaded_variation %>% n_distinct()
VEP_select_m$SYMBOL %>% as.factor() %>% n_distinct()

#join with MANE plus clinical (if not already a MANE transcript present) (566 annotation selected)
VEP_select_c = VEP_select[VEP_select$MANE_PLUS_CLINICAL != "-", ]
VEP_select_c1 = subset(VEP_select_c, !(Uploaded_variation %in% VEP_select_m$Uploaded_variation))
VEP_select_mc = rbind(VEP_select_m,VEP_select_c1)
VEP_select_mc = VEP_select_mc %>% distinct(Uploaded_variation, .keep_all = TRUE)
VEP_select_mc$Uploaded_variation %>% n_distinct()

#selecting only the refseq variants potentially not necessary if only ensembl transcript annotation - selection in web gui upload
#VEP_select_mc = VEP_select_mc[VEP_select_mc$SYMBOL_SOURCE == "HGNC", ]

###selecting single annotation for remainder of variants - those without MANE or MANE_CLINICAL----------
#identify the variants without a VEP annotation line and select their VEP annotations
VEPselectvar = VEP_select_mc$Uploaded_variation
testsetID = cisregANNOT$ID
nonmc_var = setdiff(cisregANNOT$ID,VEP_select_mc$Uploaded_variation)

#step by step inclusion by logical priory
#subset out the nonincluded variants (16 reamining)
VEP_rem= subset(VEP, (Uploaded_variation %in% nonmc_var))
as.factor(VEP_rem$Uploaded_variation) %>% n_distinct()

#selecting only the refseq variants, however potentially not necessary if only ensembl transcript annotation - selection in web gui upload
#VEP_rem = VEP_rem[VEP_rem$SYMBOL_SOURCE == "HGNC", ]
#as.factor(VEP_rem$Uploaded_variation) %>% n_distinct()

#selecting the remaining variants by 
VEP_rem_select = VEP_rem %>% filter(VEP_rem$SYMBOL == VEP_rem$ROIgene)
VEP_rem_select$Uploaded_variation %>% n_distinct()

#selected manually the clinically relevant transcript, ie standard transcript in clinvar
#APC = ENST00000507379.6 or "NM_001127511.3" as this is an additional promoter/alt transcript but relevant
VEP_select_rem_APC = filter(VEP_rem_select, Feature == "ENST00000507379.6")
VEP_select_mcc = rbind(VEP_select_mc,VEP_select_rem_APC)
VEP_select_mcc$Uploaded_variation %>% n_distinct()

#finding the transcripts for the 'non-MANE' and annotations
#remove from the VEP_select_rem list (2 to go)
VEP_rem2= subset(VEP, !(Uploaded_variation %in% VEP_select_mcc$Uploaded_variation))
#only required if only ensembl transcript annotation
VEP_rem2 = VEP_rem2[VEP_rem2$SYMBOL_SOURCE == "HGNC", ]
VEP_rem2b = VEP_rem2[VEP_rem2$MANE_SELECT != "-", ]
VEP_rem2c = VEP_rem2b %>% filter(VEP_rem2b$SYMBOL == VEP_rem2b$Selected_gene_name)
as.factor(VEP_rem2c$Uploaded_variation) %>% n_distinct()
VEP_rem2c$Uploaded_variation %>% n_distinct()

#add the variants to the ANNOT list (9542 + 7 = 9556 variants annotated)
VEP_select_mccr = rbind(VEP_select_mcc,VEP_rem2c)

VEP_rem3= subset(VEP, !(Uploaded_variation %in% VEP_select_mccr$Uploaded_variation))
#remaining - CR variants CR10 and CR512
#select the relevant annotation transcript for CR10

#CR10 no mane transcript, MANE transcript name different to ROI gene MYOC = MYOCOS
VEP_CR10 = VEP_rem3 %>% filter(Uploaded_variation =="CR10" & MANE_SELECT != "-" & SYMBOL_SOURCE == "HGNC") 
#CR512 - no mane transcript, select the gene matching the clinical report
VEP_CR512 = VEP_rem3 %>% filter(Uploaded_variation =="CR512" & SYMBOL == Selected_gene_name & SYMBOL_SOURCE == "HGNC") 

VEP_select_mcrcf = rbind(VEP_select_mccr,VEP_CR10, VEP_CR512)
VEP_select_mcrcf$Uploaded_variation %>% n_distinct()

#final check
count(subset(VEP, !(Uploaded_variation %in% VEP_select_mcrcf$Uploaded_variation)))

#saving a copy of the VEP annotations selected by the testset selected gene
write.table(VEP_select_mcrcf, "temp/cisreg_VEP_selectedannot.txt", sep = "\t", col.names = TRUE, row.names = FALSE ,quote = FALSE)

#select vep columns to add to annotation base file and add a header distinguishing
colnames(VEP_select_mcrcf)
cols2a = c("Uploaded_variation","Consequence","IMPACT","SYMBOL","Gene",
           "STRAND","Feature","EXON","INTRON","DISTANCE","X5UTR_annotation","X5UTR_consequence",
           "MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS")
VEP_select_mcrcf_cols = VEP_select_mcrcf[ ,c(cols2a)]
VEP_select_mcrcf_cols_u = distinct(VEP_select_mcrcf_cols, Uploaded_variation, .keep_all = TRUE)
summary(VEP_select_mcrcf_cols)

#rename colnames so that clear sourced from closest gene
colnames(VEP_select_mcrcf_cols) = paste("VEP", colnames(VEP_select_mcrcf_cols),sep="_")

#add columns to annotation
cisregANNOT1 = dplyr::left_join(cisregANNOT, VEP_select_mcrcf_cols, by = c("ID" = "VEP_Uploaded_variation"), keep = TRUE)

write.table(cisregANNOT1, "temp/ANNOT1_cisregall38.txt", sep = "\t", col.names = TRUE, quote = FALSE)

### ANNOTATIONS################################################
#-CADD---------------------------------------

### ANNOTATION - 2 -  

#CADD annotation addition - using grc38 v1.6 collected via INPUT.R
CADDfull= read.table("Data/GRCh38-v1.6_anno_ed4d5aceb9f7cd988ac26527738120ea.tsv", sep = "\t", header = TRUE)
colnames(CADDfull)

#make a variant key (to match to ANNOT set as not variant ID in CADD annotation)
CADDfull$location38 = paste(CADDfull$Chrom, CADDfull$Pos, sep = ":")
CADDfull$variant38ref = paste(CADDfull$location38, CADDfull$Ref, sep = "")
CADDfull$variant38 = paste(CADDfull$variant38ref, CADDfull$Alt, sep = ">")
CADDfull$variant38 %>% n_distinct()

CADDrefset = subset(CADDfull, (variant38 %in% c(refset$variant38)))
n_distinct(CADDrefset$variant38 %>% unique())

#first select the columns to add
#cisregDraw_tsin = select(cisregDraw)
CADD = unique(CADDrefset[ ,c(11:16,34:35,80:81,104:105,128,131:137)])

#adding datasource to dataframe header
colnames(CADD) = paste("CADD", colnames(CADD),sep="_")
cisregANNOT2 = left_join(cisregANNOT1, CADD, by = c("variant38" = "CADD_variant38"))

cisregANNOT2$ID %>% n_distinct()
cisregANNOT2$variant38 %>% n_distinct()

#count number of NAs
sum(is.na(cisregANNOT2$CADD_PHRED))

#-Linsight - hg19------------------------------------------
#annotation by rehan via VEP cmd line on hg19 custom annotation, extracted by Maddison from resulting VCF in R,  28.06.23
#load raw annotation
linsightfull= read.table("Data/linsightscoresFINAL_df.txt", sep = "\t", header = TRUE)
summary(linsightfull$Lscore_hg19)

#select columns to add to annotation base file and add a header distinguishing
colnames(linsightfull)
colnames(linsightfull) = paste(colnames(linsightfull),"hg19",sep="_")
ls_cols = linsightfull[ ,c("ID_hg19","Lscore_hg19")]
ls_cols$ID %>% n_distinct()
summary(ls_cols)

#rename colnames so that clear sourced from closest gene
colnames(ls_cols) = paste("linsight", colnames(ls_cols),sep="_")

#add columns to annotation
cisregANNOT3 = dplyr::left_join(cisregANNOT2, ls_cols, by = c("ID" = "linsight_ID_hg19"), keep = TRUE)
cisregANNOT3$ID %>% n_distinct()

#check by comparing original result
nrow(linsightfull)
summary(linsightfull$Lscore_hg19)
summary(cisregANNOT3)

#-FATHMM-MKL--------------------------------------------------

#create input format(e,g,11,247320,G,T on grc37)
refset_liftover$FATHMMinput37 = paste(refset_liftover$CHROM37,refset_liftover$start37,refset_liftover$REF,refset_liftover$ALT, sep = ",")
FATHMMinput37 = as.data.frame(refset_liftover$FATHMMinput37)

#save FATHMM-MKL input file
write.table(FATHMMinput37, "temp/FATHMMinput37.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#annotation via web gui by Rehan grc37 input - https://fathmm.biocompute.org.uk/fathmm-xf/
#load raw annotation
FATHMMMKLfull= read.table("Data/FATHMM-MKLresults_grc37_271023.txt", sep = "\t", header = TRUE)

#exclude NA's/empty scores
FATHMMMKL = FATHMMMKLfull[!is.na(FATHMMMKLfull$Non.Coding.Score), ]

#make a variant key (to match to ANNOT set as not variant ID in FATHMMKL annotation)
colnames(FATHMMMKL)
colnames(FATHMMMKL) = paste(colnames(FATHMMMKL),"grc37",sep="_")
FATHMMMKL$location_grc37 = paste(FATHMMMKL$Chromosome_grc37, FATHMMMKL$Position_grc37, sep = ":")
FATHMMMKL$variant_grc37ref = paste(FATHMMMKL$location_grc37, FATHMMMKL$Ref..Base_grc37, sep = "")
FATHMMMKL$variant_grc37 = paste(FATHMMMKL$variant_grc37ref, FATHMMMKL$Mutant.Base_grc37, sep = ">")

#select columns to add to annotation base file and add a header distinguishing
fthm_cols = FATHMMMKL[ ,c("variant_grc37","Non.Coding.Score_grc37","Non.Coding.Groups_grc37","Coding.Score_grc37","Coding.Groups_grc37")]
fthm_cols$variant_grc38 %>% n_distinct()
summary(fthm_cols)

#rename colnames so that clear sourced from closest gene
colnames(fthm_cols) = paste("fthmmkl", colnames(fthm_cols),sep="_")

#add columns to annotation
cisregANNOT4 = dplyr::left_join(cisregANNOT3, fthm_cols, by = c("variant37" = "fthmmkl_variant_grc37"), keep = TRUE)
cisregANNOT4$ID %>% n_distinct()
summary(cisregANNOT4)

#-FATHMM-XF---------------------------------------------------------------------------------------------

#annotation by Rehan 26.10.23 - grc38
fathmxf= read.table("Data/FATHMM-XFresults_261023.txt", sep = "\t", header = TRUE)

#exclude the rows with NA (in fathmxf this is designated as a '--'), these are all indels
colnames(fathmxf)
fathmxf = fathmxf[fathmxf$Non.Coding.Score != "--", ]

#make a variant key (to match to ANNOT set as not variant ID in FATHMMKL annotation)
colnames(fathmxf) = paste(colnames(fathmxf),"grc38",sep="_")
colnames(fathmxf)
fathmxf$location_grc38 = paste(fathmxf$Chromosome_grc38, fathmxf$Position_grc38, sep = ":")
fathmxf$variant_grc38ref = paste(fathmxf$location_grc38, fathmxf$Ref..Base_grc38, sep = "")
fathmxf$variant_grc38 = paste(fathmxf$variant_grc38ref, fathmxf$Mutant.Base_grc38, sep = ">")

#select columns to add to annotation base file and add a header distinguishing
fthx_cols = fathmxf[ ,c("variant_grc38","Non.Coding.Score_grc38","Warning_grc38")]
fthx_cols$Non.Coding.Score_grc38 = as.numeric(fthx_cols$Non.Coding.Score_grc38)
fthx_cols$variant_grc38 %>% n_distinct()
summary(fthx_cols)

#rename colnames so that clear sourced from closest gene
colnames(fthx_cols) = paste("fthmxf", colnames(fthx_cols),sep="_")

#add columns to annotation
cisregANNOT4 = dplyr::left_join(cisregANNOT4, fthx_cols, by = c("variant38" = "fthmxf_variant_grc38"), keep = TRUE)
summary(cisregANNOT4)

#checks
summary(as.numeric(fathmxf$Non.Coding.Score_grc38))
summary(as.numeric(cisregANNOT4$fthmxf_Non.Coding.Score_grc38))
nrow(cisregANNOT) - nrow(fathmxf)
cisregANNOT4$ID %>% n_distinct()

#-Eigen ---------------------------------------------------------------------------------------------
#annotation by Rehan, cmd line vep custom annotation on hg38, 5.7.23 plus chr 8 repeated 10.07.23
#extracted by Maddison, cmd line results combined 
eigenfull= read.table("Data/Eigenscores_results_full13k_df_RV.txt", sep = "\t", header = TRUE)
eigen_refset = subset(eigenfull, (ID %in% c(refset$ID)))

#exclude the rows with NA (in eigen this is designated as a '--'), these are all indels
colnames(eigen_refset)
summary(eigen_refset)

#select columns to add to annotation base file and add a header distinguishing
eigen_cols = eigen_refset[ ,c("ID","Eigen.raw","Eigen.phred","Eigen.PC.raw","Eigen.PC.phred")]
eigen_cols$Eigen.raw = as.numeric(eigen_cols$Eigen.raw)
eigen_cols$Eigen.phred = as.numeric(eigen_cols$Eigen.phred)
eigen_cols$Eigen.PC.raw = as.numeric(eigen_cols$Eigen.PC.raw)
eigen_cols$Eigen.PC.phred = as.numeric(eigen_cols$Eigen.PC.phred)
summary(eigen_cols)

#rename colnames so that clear sourced from closest gene
colnames(eigen_cols) = paste("eigen", colnames(eigen_cols),sep="_")

#add columns to annotation
colnames(eigen_cols)
cisregANNOT4 = dplyr::left_join(cisregANNOT4, eigen_cols, by = c("ID" = "eigen_ID"))
cisregANNOT4$ID %>% n_distinct()
summary(cisregANNOT4)

#-REMM---------------------------------------------------------------------------------------------
#annotation by maddison, 6.6.23 - hg38 + manual correction of 4 indels using sngl lookup
REMM= read.table("Data/Full REMM score 13k variants_df2.txt", sep = "\t", header = TRUE)

# check no rows with NA
summary(REMM$score)

#dvariant key (to match to ANNOT)
REMM$location38 = paste(REMM$chrom, REMM$location, sep = ":")
REMM$variant38ref = paste(REMM$location38, REMM$ref, sep = "")
REMM$variant38 = paste(REMM$variant38ref, REMM$alt, sep = ">")

#select columns to add to annotation base file and add a header distinguishing
REMM_cols = REMM[ ,c("grc38","score")]
REMM_cols$grc38 %>% n_distinct()
summary(REMM_cols)

#rename colnames so that clear sourced from closest gene
colnames(REMM_cols) = paste("REMM", colnames(REMM_cols),sep="_")

#add columns to annotation
cisregANNOT4 = dplyr::left_join(cisregANNOT4, REMM_cols, by = c("variant38" = "REMM_grc38"), keep = TRUE)

#checks
cisregANNOT4$ID %>% n_distinct()
summary(cisregANNOT4$REMM_score)
cisregANNOT4[is.na(cisregANNOT4$REMM_score),]

write.table(cisregANNOT4, "temp/cisreg_ANNOT4.txt", sep = "\t", col.names = TRUE, quote = FALSE)

###############################################################
#SAVE ALL ANNOTATION SCORES 
#################################################################

cisregANNOT_c = cisregANNOT4

#save complete annotation file, cisreg_annot, add dataset, add date to indicate most recent version, save in "Current_cisreg_annot" folder
write.table(cisregANNOT_c, "temp/cisreg_testset_ANNOTallscores.txt", sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

#list annotations/columns for creating data dictionary
annotations = as.data.frame(colnames(cisregANNOT_c))
write.table(annotations, "temp/cisreg_testset_ANNOTallscore_INFO.txt", sep = "\t", col.names = TRUE, quote = FALSE)

#save datatable for wiggins analysis
colnames(cisregANNOT4)
cisreg_wiggins = cisregANNOT4[ , c("ID","INFO","CADD_PHRED", "linsight_Lscore_hg19",
                                   "fthmmkl_Non.Coding.Score_grc37","fthmxf_Non.Coding.Score_grc38",
                                   "eigen_Eigen.raw","REMM_score")]

#converting INFO/exp column into binary o-1 for wiggins tool
cisreg_wiggins$observed.event = if_else(cisreg_wiggins$INFO == "Cont", 0,1)

#save files - csv for wiggins tool
write.csv(cisreg_wiggins, "temp/cisreg_refset_wigginsinput.csv", row.names = FALSE)

#save tool specific data for wiggins input
cisreg_wiggins_eigen = cisreg_wiggins[ , c("ID","INFO","observed.event","eigen_Eigen.raw")]
cisreg_wiggins_eigen <- cisreg_wiggins_eigen[rowSums(is.na(cisreg_wiggins_eigen)) == 0, ]
write.csv(cisreg_wiggins, "temp/wigginsinput_eigen.csv", row.names = FALSE)


###################################################################
##Creating score categories
###################################################################
#--CADD categories-----------------------------------------------
#translating scores to categories
colnames(cisregANNOT_c)

#CADD - score designated via CADD PHRED
#thresholds - defined based on the independent calibration via wiggins scores
CADD_upper_threshold = 10
CADD_lower_threshold = 8

#designate categories(via bins) for each testset
summary(cisregANNOT_c$CADD_PHRED)
CADD_binvalues = c((min(na.omit(cisregANNOT_c$CADD_PHRED))), CADD_lower_threshold, CADD_upper_threshold, (max(na.omit(cisregANNOT_c$CADD_PHRED))))
binlabels = c("score_low","score_no_call","score_high")
cisregANNOT_c$CADD_scorebin  = cut(cisregANNOT_c$CADD_PHRED, breaks = CADD_binvalues, labels = binlabels)
cisregANNOT_c %>% group_by(CADD_scorebin) %>% summary()

#create 6 bins for CADD score
cisregANNOT_c$CADD_testbin  = as.factor(paste(cisregANNOT_c$INFO,cisregANNOT_c$CADD_scorebin ,sep = "_"))
cisregANNOT_c$CADD_testbin %>% summary()

#--linsight categories--------------------------------------------------------
#linsight - score designated via linsight_Lscore_hg19
#thresholds - defined based on the independent calibration via wiggins scores
linsight_upper_threshold = 0.24
linsight_lower_threshold = 0.16

#designate categories(via bins) for each testset
summary(cisregANNOT_c$linsight_Lscore_hg19)
linsight_binvalues = c((min(na.omit(cisregANNOT_c$linsight_Lscore_hg19))), linsight_lower_threshold, linsight_upper_threshold, (max(na.omit(cisregANNOT_c$linsight_Lscore_hg19))))
#binlabels = c("score_ben","score_int","score_path")
cisregANNOT_c$linsight_scorebin  = cut(cisregANNOT_c$linsight_Lscore_hg19, breaks = linsight_binvalues, labels = binlabels)
cisregANNOT_c %>% group_by(linsight_scorebin) %>% summary()

#create 6 bins for linsight score
cisregANNOT_c$linsight_testbin  = as.factor(paste(cisregANNOT_c$INFO,cisregANNOT_c$linsight_scorebin ,sep = "_"))
cisregANNOT_c$linsight_testbin %>% summary()

#--fathmmkl categories--------------------------------------------------------
#fathmmlk - score designated via fthmmkl_Non.Coding.Score_grc37
#thresholds - defined based on the independent calibration via wiggins scores
fathmmlk_upper_threshold = 0.59
fathmmlk_lower_threshold = 0.39

#designate categories(via bins) for each testset
summary(cisregANNOT_c$fthmmkl_Non.Coding.Score_grc37)
fthmmkl_binvalues = c((min(na.omit(cisregANNOT_c$fthmmkl_Non.Coding.Score_grc37))), fathmmlk_lower_threshold, fathmmlk_upper_threshold, (max(na.omit(cisregANNOT_c$fthmmkl_Non.Coding.Score_grc37))))
#binlabels = c("score_ben","score_int","score_path")
cisregANNOT_c$fathmmlk_scorebin  = cut(cisregANNOT_c$fthmmkl_Non.Coding.Score_grc37, breaks = fthmmkl_binvalues, labels = binlabels)
cisregANNOT_c %>% group_by(fathmmlk_scorebin) %>% summary()

#create 6 bins for fathmmlk score
cisregANNOT_c$fathmmlk_testbin  = as.factor(paste(cisregANNOT_c$INFO,cisregANNOT_c$fathmmlk_scorebin ,sep = "_"))
cisregANNOT_c$fathmmlk_testbin %>% summary()

#--fathmxf categories--------------------------------------------------------
#fathmxf - score designated via fthmxf_Non.Coding.Score_grc38
#thresholds - defined based on the independent calibration via wiggins scores
fathmxf_upper_threshold = 0.14
fathmxf_lower_threshold = 0.12

#designate categories(via bins) for each testset
summary(cisregANNOT_c$fthmxf_Non.Coding.Score_grc38)
fathmxf_binvalues = c((min(na.omit(cisregANNOT_c$fthmxf_Non.Coding.Score_grc38))), fathmxf_lower_threshold, fathmxf_upper_threshold, (max(na.omit(cisregANNOT_c$fthmxf_Non.Coding.Score_grc38))))
#binlabels = c("score_ben","score_int","score_path")
cisregANNOT_c$fathmxf_scorebin  = cut(cisregANNOT_c$fthmxf_Non.Coding.Score_grc38, breaks = fathmxf_binvalues, labels = binlabels)
cisregANNOT_c %>% group_by(fathmxf_scorebin) %>% summary()

#create 6 bins for fathmxf score
cisregANNOT_c$fathmxf_testbin  = as.factor(paste(cisregANNOT_c$INFO,cisregANNOT_c$fathmxf_scorebin ,sep = "_"))
cisregANNOT_c$fathmxf_testbin %>% summary()

#--eigen categories---------------------------------------------------------
#eigen - score designated via eigen_Eigen.raw
#thresholds - defined based on the independent calibration via wiggins scores
eigen_upper_threshold = 0.594
eigen_lower_threshold = 0.394

#designate categories(via bins) for each testset
summary(cisregANNOT_c$eigen_Eigen.raw)
eigen_binvalues = c((min(na.omit(cisregANNOT_c$eigen_Eigen.raw))), eigen_lower_threshold, eigen_upper_threshold, (max(na.omit(cisregANNOT_c$eigen_Eigen.raw))))
#binlabels = c("score_ben","score_int","score_path")
cisregANNOT_c$eigen_scorebin  = cut(cisregANNOT_c$eigen_Eigen.raw, breaks = eigen_binvalues, labels = binlabels)
cisregANNOT_c %>% group_by(eigen_scorebin) %>% summary()

#create 6 bins for fathmxf score
cisregANNOT_c$eigen_testbin  = as.factor(paste(cisregANNOT_c$INFO,cisregANNOT_c$eigen_scorebin ,sep = "_"))
cisregANNOT_c$eigen_testbin %>% summary()

#--ReMM categories--------------------------------------------------------
#REMM - score designated via REMM_score
#thresholds - defined based on the independent calibration via wiggins scores
REMM_upper_threshold = 0.86
REMM_lower_threshold = 0.80

#designate categories(via bins) for each testset
summary(cisregANNOT_c$REMM_score)
REMM_binvalues = c((min(na.omit(cisregANNOT_c$REMM_score))), REMM_lower_threshold, REMM_upper_threshold, (max(na.omit(cisregANNOT_c$REMM_score))))
#binlabels = c("score_ben","score_int","score_path")
cisregANNOT_c$REMM_scorebin  = cut(cisregANNOT_c$REMM_score, breaks = REMM_binvalues, labels = binlabels, include.lowest = TRUE)
cisregANNOT_c %>% group_by(REMM_scorebin) %>% summary()

#create 6 bins for fathmxf score
cisregANNOT_c$REMM_testbin  = as.factor(paste(cisregANNOT_c$INFO,cisregANNOT_c$REMM_scorebin ,sep = "_"))
cisregANNOT_c$REMM_testbin %>% summary()

#summary-of-results------------------------------------

summary_eval = cisregANNOT_c[ , c("ID","INFO","CADD_scorebin","CADD_testbin", "linsight_scorebin" ,"linsight_testbin",
                                  "fathmmlk_scorebin", "fathmmlk_testbin","fathmxf_scorebin","fathmxf_testbin",
                                  "eigen_scorebin","eigen_testbin","REMM_scorebin", "REMM_testbin")]
summary(summary_eval)

#summary tables - don't include NA's currently
table(summary_eval$INFO, summary_eval$CADD_scorebin)
table(summary_eval$INFO, summary_eval$linsight_scorebin)
table(summary_eval$INFO, summary_eval$fathmmlk_scorebin)
table(summary_eval$INFO, summary_eval$fathmxf_scorebin)
table(summary_eval$INFO, summary_eval$eigen_scorebin)
table(summary_eval$INFO, summary_eval$REMM_scorebin)

#---------------------------------------------------
###############################################################
#SAVE  ANNOTATION CATEORIES 
#################################################################

#save complete annotation file, cisreg_annot for continued  analysis
write.table(cisregANNOT_c, "temp/cisreg_testset_ANNOTcategories.txt", sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

#save complete annotation file, for annotation file
write.table(cisregANNOT_c, "Data/cisreg_testset_ANNOTcategories.txt", sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

#list annotations/columns for creating data dictionary
annotations = as.data.frame(colnames(cisregANNOT_c))
write.table(annotations, "Data/cisreg_testset_ANNOTallscore_INFO.txt", sep = "\t", col.names = TRUE, quote = FALSE)


####################################################################

####################basic reference set characteristics
#How many variants in test set scored
#Control n
cont_n = sum(cisregANNOT_c$INFO=="Cont")

#disease n
dis_n = sum(cisregANNOT_c$INFO=="Dis")

#test n
all_n = nrow(cisregANNOT_c)

###CADD analysis ########################
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

#extract a table of CADD_scorebin LR, LR CI low and LR CIhigh
#creating new groups
colnames(cisregANNOT_c)

#need to convert this to a table
summary(as.factor(cisregANNOT_c$CADD_scorebin))

cisregANNOT_c %>% select("INFO","CADD_scorebin") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$CADD_scorebin), ] %>% select(., "INFO","CADD_scorebin")

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
test2_TP = df_test2 %>% filter(INFO=="Dis" & CADD_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(INFO=="Cont" & CADD_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(INFO=="Cont" & CADD_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(INFO=="Dis" & CADD_scorebin == "score_low") %>% count()
print(test2_FN)


#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(INFO=="Dis" & CADD_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(INFO=="Cont" & CADD_scorebin == "score_no_call") %>% count()

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

###REMM analysis#############################################################
colnames(cisregANNOT_c)
#REMM_score
summary(cisregANNOT_c$REMM_score)
summary(REMM$score)
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

#extract a table of REMM_score LR, LR CI low and LR CIhigh
#creating new groups
colnames(cisregANNOT_c)

#need to convert this to a table
summary(as.factor(cisregANNOT_c$REMM_score))

cisregANNOT_c %>% select("INFO","REMM_score") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$REMM_scorebin), ] %>% select(., "INFO","REMM_scorebin")

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
test2_TP = df_test2 %>% filter(INFO=="Dis" & REMM_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(INFO=="Cont" & REMM_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(INFO=="Cont" & REMM_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(INFO=="Dis" & REMM_scorebin == "score_low") %>% count()
print(test2_FN)


#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(INFO=="Dis" & REMM_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(INFO=="Cont" & REMM_scorebin == "score_no_call") %>% count()

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

#determine LR for each category, plus CI low and LR CIhigh
#creating new groups
colnames(cisregANNOT_c)

#need to convert this to a table
summary(as.factor(cisregANNOT_c$fathmmlk_scorebin))

cisregANNOT_c %>% select("INFO","fathmmlk_scorebin") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$fathmmlk_scorebin), ] %>% select(., "INFO","fathmmlk_scorebin")

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

#test-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(INFO=="Dis" & fathmmlk_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(INFO=="Cont" & fathmmlk_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(INFO=="Cont" & fathmmlk_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(INFO=="Dis" & fathmmlk_scorebin == "score_low") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(INFO=="Dis" & fathmmlk_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(INFO=="Cont" & fathmmlk_scorebin == "score_no_call") %>% count()

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
test2_varLRpos = (1/test2_TN-1/test2_cont_n)+(1/test2_FN-1/test2_dis_n)
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

#determine LR for each category, plus CI low and LR CIhigh
#creating new groups
colnames(cisregANNOT_c)

#need to convert this to a table
summary(as.factor(cisregANNOT_c$fathmxf_scorebin))

cisregANNOT_c %>% select("INFO","fathmxf_scorebin") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$fathmxf_scorebin), ] %>% select(., "INFO","fathmxf_scorebin")

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

#test-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(INFO=="Dis" & fathmxf_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(INFO=="Cont" & fathmxf_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(INFO=="Cont" & fathmxf_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(INFO=="Dis" & fathmxf_scorebin == "score_low") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(INFO=="Dis" & fathmxf_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(INFO=="Cont" & fathmxf_scorebin == "score_no_call") %>% count()

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
test2_varLRpos = (1/test2_TN-1/test2_cont_n)+(1/test2_FN-1/test2_dis_n)
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

#determine LR for each category, plus CI low and LR CIhigh
#creating new groups
colnames(cisregANNOT_c)

#need to convert this to a table
summary(as.factor(cisregANNOT_c$eigen_scorebin))

cisregANNOT_c %>% select("INFO","eigen_scorebin") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$eigen_scorebin), ] %>% select(., "INFO","eigen_scorebin")

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

#test-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(INFO=="Dis" & eigen_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(INFO=="Cont" & eigen_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(INFO=="Cont" & eigen_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(INFO=="Dis" & eigen_scorebin == "score_low") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(INFO=="Dis" & eigen_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(INFO=="Cont" & eigen_scorebin == "score_no_call") %>% count()

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
test2_varLRpos = (1/test2_TN-1/test2_cont_n)+(1/test2_FN-1/test2_dis_n)
print(test2_varLRpos)

#calc low CI for the LR
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

#determine LR for each category, plus CI low and LR CIhigh
#creating new groups
colnames(cisregANNOT_c)

#need to convert this to a table
summary(as.factor(cisregANNOT_c$linsight_scorebin))

cisregANNOT_c %>% select("INFO","linsight_scorebin") %>% table()
#selecting the relevant results for the test evaluation from scores - categories designated in previous script based on scores/process
df_test2 = cisregANNOT_c[!is.na(cisregANNOT_c$linsight_scorebin), ] %>% select(., "INFO","linsight_scorebin")

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

#test-evaluation numbers

#Calculate n True Pos (Event pos)
test2_TP = df_test2 %>% filter(INFO=="Dis" & linsight_scorebin == "score_high") %>% count()
print(test2_TP)

#calculate n True neg (no event neg)
test2_TN = df_test2 %>% filter(INFO=="Cont" & linsight_scorebin == "score_low") %>% count()
print(test2_TN)

#calculate false positive (no event pos)
test2_FP = df_test2 %>% filter(INFO=="Cont" & linsight_scorebin == "score_high") %>% count()
print(test2_FP)

#calculate false negative (event neg)
test2_FN = df_test2 %>% filter(INFO=="Dis" & linsight_scorebin == "score_low") %>% count()
print(test2_FN)

#Calculate Uncertain (event no call)
test2_UNP = df_test2 %>% filter(INFO=="Dis" & linsight_scorebin == "score_no_call") %>% count()

#Calculate Uncertain (no event event no call )
test2_UNN = df_test2 %>% filter(INFO=="Cont" & linsight_scorebin == "score_no_call") %>% count()

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
test2_varLRpos = (1/test2_TN-1/test2_cont_n)+(1/test2_FN-1/test2_dis_n)
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


###############################################################
#reg_variant_tool_evaluation - LRs
#import summary table of all LR and CIs
#thresholds - defined based on the independent calibration via wiggins score
#graph the LRs of each group

#data sourced from Supplemental table 10_Combined_score_evaluation
cisreg_eval = read.table("Data/reg_variant_tool_evaluation - summary.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
colnames(cisreg_eval)
cisregLRs = read.table("Data/reg_variant_tool_evaluation - LRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

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

legend = get_legend(ggplot(cisregANNOT_c, aes(CADD_PHRED, fill = INFO))+
                      geom_density(position = "identity", alpha = 0.2)+
                      labs(x = "CADD", y = "density")+
                      theme(legend.position = "bottom", legend.title = element_blank())+
                      geom_vline(xintercept = CADD_lower_threshold) + geom_vline(xintercept = CADD_upper_threshold)+
                      scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value="#999999"))

ggdraw(plot_grid(cadd_density, cadd_LRplot,
                 remm_density, remm_LRplot,
                 fathmmkl_density, fathmmkl_LRplot,
                 fathmmxf_density, fathmmxf_LRplot,
                 eigen_density, eigen_LRplot,
                 linsight_density, linsight_LRplot,
                 align = "hv", axis = "b", rel_widths = c(0.75, 1.5),
                 labels = c("A","B", "C", "D","E", "F", "G", "H", "I", "J", "K", "L" ),
                 ncol = 2, nrow = 6))

ggsave(path = "D:/cisreg_manuscript/figures", filename = "Fig4.png", width = 10, height = 12, device='png', dpi=1200)

ggdraw(plot_grid(legend))
ggsave(path = "D:/cisreg_manuscript/figures", filename = "Fig4_legend.png", width = 10, height = 12, device='png', dpi=1200)

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
