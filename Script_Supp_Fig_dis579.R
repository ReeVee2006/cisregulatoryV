##cis reg disease variant annotation
##R code complemented by base annotation datasets
##last executed 25.10.23
##FIGURE 3 - cisreg manuscript

library(tidyverse)
library(GenomicRanges)
library(cowplot)

#--------------------------------------------
#setwd 
setwd("D:/cisreg_manuscript/SUPP")

## Load datasets to annotate. Cisreg variants, control and disease to combine dataset and create into VCF file
#load disease-associated variant set
cisregDraw = read.table("Tables/SUPP_Table3_disvar.txt", header= TRUE, sep = "\t")

#create a list of genes represented in disease variant set
testset_genes_all = unique(cisregDraw$Source_Reported_GeneName)

##--------------------------------------------------------

#select disease vcf columns, IDs with project identifiers CRD
cisregDis = cisregDraw[ ,c("Chr_GRCh38","Location_GRCh38","Testset_ID","Reference_Allele","Alternative_Allele")]
colnames(cisregDis) = c("CHROM", "POS", "ID", "REF", "ALT")

#add vcf columns to disease set
cisregDis$QUAL = "."
cisregDis$FILTER = "."
cisregDis$INFO = "Dis"

#sort by chromosome and then location
cisregDis = arrange(cisregDis, CHROM, POS)
cisregDis_vcf = cisregDis

# save a copy of the vcf file. 
#write.table(cisregDis_vcf, "temp/cisreg_dis_grc38_vcf.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# add the header to align with required vcf file format
##fileformat=VCFv4.2
##fileDate= <date>
##source script=Fig3_dis579
##reference=cisregDis_vcf
#

#annotate with VEP web GUI

#---------------------------------------------------
#ANNOTATION adding base information

## ANNOTATION exp.group - check the exp.group column is filled 
cisregDis$INFO %>% as.factor() %>% summary()

#ANNOTATION location38 and variant38 - columns giving location and variant specific values
cisregDis$location38 = paste(cisregDis$CHROM, cisregDis$POS, sep = ":")
cisregDis$variant38ref = paste(cisregDis$location38, cisregDis$REF, sep = "")
cisregDis$variant38 = paste(cisregDis$variant38ref, cisregDis$ALT, sep = ">")

#calculating end position - support determination of 1based including all variant types
cisregDis$ALTcount = nchar(cisregDis$ALT)
cisregDis$REFcount = nchar(cisregDis$REF)
cisregDis$allelesize = cisregDis$ALTcount - cisregDis$REFcount
cisregDis$width = cisregDis$allelesize +1
cisregDis$end.field = cisregDis$POS + cisregDis$REFcount

#--------------------------------------------------------
#create GRanges object from cisreg  variant list for input to annotations

#converting dataframe into GRanges object - from vcf like file
cisregDis$chr_38 = paste("chr",cisregDis$CHROM,sep="")
grcisregDis =  makeGRangesFromDataFrame(cisregDis,
                                          keep.extra.columns=FALSE,
                                          ignore.strand=TRUE,
                                          seqinfo=NULL,
                                          seqnames.field=c("chr_38"),
                                          start.field=c("POS"),
                                          end.field=c("end.field"),
                                          starts.in.df.are.0based=FALSE)

#-------------------------------------------------------------------------
##ANNOTATIONS - 1 -  bespoke additions from the original MAIN data

#adding additional columns from the original raw variant list dataset - depending on dataset modify the column headings
#first select the columns to add
#cisregDraw_tsin = select(cisregDraw)
colnames(cisregDraw)
colsraw = cisregDraw[ ,c("Testset_ID",
                         "GRCh38variant_B",
                         "Selected_gene_name","Source_Reported_GeneName",
                         "GRCh37variant")]

#adding datasource to dataframe header
colnames(colsraw) = paste("raw", colnames(colsraw),sep="_")

cisregANNOT1 = left_join(cisregDis, colsraw, by = c("ID" = "raw_Testset_ID"))

##--------ANNOTATIONS--2---VEP-----########################## 
# VEP annotation - currently input above file into web GUI, transcript database Ensembl, dont include protein annotation, phenotype yes to all, add spliceAI and conservation - currently only AA and BLUOSUM

#VEP annotation performed - 
#remove # at start of VEP annotation to read
VEPfull = read.delim("Data/dis579_VEP_KnYrlYkXXEctn0hs.txt", sep = "\t", header= TRUE)
colnames(VEPfull)

#count the number of unique variants based on project identified "CR#"
VEPfull$Uploaded_variation %>% as.factor() %>% n_distinct()
VEPfull$MANE_SELECT %>% as.factor() %>% n_distinct()
VEPfull$Gene %>% as.factor() %>% n_distinct()
VEPfull$testsetID = VEPfull$Uploaded_variation

#Analysis of only Refseq annotations - annotation gene and transcripts
VEPrefseq = VEPfull[VEPfull$SYMBOL_SOURCE == "EntrezGene", ]
VEP = VEPfull[VEPfull$SYMBOL_SOURCE == "HGNC", ]

VEP$Uploaded_variation %>% as.factor() %>% n_distinct()
VEP$MANE_SELECT %>% as.factor() %>% n_distinct()
VEP$MANE_PLUS_CLINICAL %>% as.factor() %>% n_distinct()
VEP$Gene %>% as.factor() %>% n_distinct()
VEP$SYMBOL %>% as.factor() %>% n_distinct()
VEP$Feature %>% as.factor() %>% n_distinct()

#variants annotating to multiple transcripts
VEP_mulitrans = VEP %>% group_by(Uploaded_variation) %>% filter(n() >1)
VEP_mulitrans$Uploaded_variation %>% as.factor() %>% n_distinct()
VEP_mulitrans$SYMBOL %>% as.factor() %>% n_distinct()

#variants annotating to sngl transcript
VEP_sng = VEP %>% group_by(Uploaded_variation) %>% filter(n() ==1)
VEP_sng$Uploaded_variation %>% as.factor() %>% n_distinct()

# looking at the number of selected genes represented in the disease set based on ALL annotations
VEP$SYMBOL %>% as.factor() %>% n_distinct()
genecount = data.frame(table(VEP$SYMBOL))
colnames(genecount)= c("SYMBOL", "count")
transcriptcount = data.frame(table(VEP$Feature))
summary(as.factor(transcriptcount$Var1))

#count the number of transcripts per variant
VEPrs_trans =  VEP %>% group_by(Uploaded_variation) %>% count()

ggplot(VEPrs_trans, aes(n)) +
  geom_histogram(binwidth = 1)+ 
  labs(x="# transcripts",  y = "# variants")

#-------------------------------------------------------------
#selecting the VEP annotation according to the MAIN dataset selected gene
#ensure that "Selected_gene_name" column in present on the ANNOT file - currently sourced from manual annotation

#create a table with the selected gene symbol and gene information for each variant
colnames(cisregDraw)
var_geneinfo = cisregDraw[ ,c("Testset_ID","GRCh38variant_B","Selected_gene_name")]

# annotate the VEP file with the selected gene name
VEPfull = left_join(x = VEPfull, y= var_geneinfo, by = c("Uploaded_variation" = "Testset_ID"))

#select the VEP annotation rows that match the selected gene 
VEP_select = VEPfull %>% filter(VEPfull$SYMBOL == VEPfull$Selected_gene_name)
VEPfull$SYMBOL %>% as.factor() %>% n_distinct()
VEP_select$SYMBOL %>% as.factor() %>% n_distinct()
VEP_select$Uploaded_variation %>% as.factor() %>% n_distinct()
#select MANEselect transcripts only
VEP_select_m = VEP_select[VEP_select$MANE_SELECT != "-", ]
VEP_select_m$Uploaded_variation %>% n_distinct()
VEP_select_m$SYMBOL %>% as.factor() %>% n_distinct()

#join with MANE plus 
VEP_select_c = VEP_select[VEP_select$MANE_PLUS_CLINICAL != "-", ]
VEP_select_mc = rbind(VEP_select_m,VEP_select_c)
VEP_select_mc$Uploaded_variation %>% n_distinct()

#selecting only the refseq variants, however potentially not necessary if only ensembl transcript annotation - selection in web gui upload
VEP_select_mc = VEP_select_mc[VEP_select_mc$SYMBOL_SOURCE == "HGNC", ]
VEP_select_mc = VEP_select_mc %>% distinct(Uploaded_variation, .keep_all = TRUE)
VEP_select_mc$Uploaded_variation %>% n_distinct()

#identify the variants without a VEP annotation line and select their VEP annotations
VEPselectvar = VEP_select_mc$Uploaded_variation
testsetID = cisregDis$ID
nonmc_var = setdiff(cisregDis$ID,VEP_select_mc$Uploaded_variation)

#step by step inclusion by logical priory
#subset out the nonincluded variants (17 not identifed)
VEP_rem= subset(VEPfull, (Uploaded_variation %in% nonmc_var))
as.factor(VEP_rem$Uploaded_variation) %>% n_distinct()

##DEV TST - try selecting the canonical transcript annotation for the 17 remaining, # does that solve

#selecting the remaining variants by gene and protein coding transcript
VEP_rem_select = VEP_rem %>% filter(VEP_rem$SYMBOL == VEP_rem$Selected_gene_name)
VEP_rem_select$Uploaded_variation %>% n_distinct()
VEP_rem_select = VEP_rem_select[VEP_rem_select$SYMBOL_SOURCE == "HGNC", ]
VEP_rem_select$Uploaded_variation %>% n_distinct()

#select the variants with only a single annotation (3)
VEP_rem_sng = VEP_rem_select %>% group_by(Uploaded_variation) %>% filter(n() ==1)
VEP_select_mcr = rbind(VEP_select_mc,VEP_rem_sng)
VEP_select_mcr$Uploaded_variation %>% n_distinct()

VEP_rem_dup = VEP_rem_select %>% group_by(Uploaded_variation) %>% filter(n() >1)
VEP_rem_dup$Uploaded_variation %>% n_distinct()
VEP_rem_dup2 = subset(VEP_rem_dup, BIOTYPE == "protein_coding")
VEP_rem_dup2$Uploaded_variation %>% n_distinct()

#the final transcripts/annotations to be selected by clinically relevant, ie standard transcript in clinvar
#APC = ENST00000507379.5 or "NM_001127511.3"
VEP_select_rem2 = filter(VEP_rem_dup2, Feature == "ENST00000507379.6")

#SHOX = ENST00000381578.6 "NM_001127511.3" (CR473, 474, 475)
VEP_select_rem3 = filter(VEP_rem_dup2, Feature == "ENST00000381578.6")

VEP_select_mcr4 = rbind(VEP_select_mcr,VEP_select_rem2, VEP_select_rem3)
VEP_select_mcr4$Uploaded_variation %>% n_distinct()

#remove from the VEP_select_rem list (2 to go)
VEP_rem4= subset(VEPfull, !(Uploaded_variation %in% VEP_select_mcr4$Uploaded_variation))
VEP_rem4$Uploaded_variation %>% n_distinct()
VEP_rem4 = VEP_rem4[VEP_rem4$MANE_SELECT != "-", ]
VEP_rem4 = VEP_rem4[VEP_rem4$SYMBOL_SOURCE == "HGNC", ]
VEP_rem4$Uploaded_variation %>% n_distinct()

VEP_select_mcr5 = rbind(VEP_select_mcr4,VEP_rem4)

VEP_select_mcr5$Uploaded_variation %>% n_distinct()
VEP_select_mcr5 = distinct(VEP_select_mcr5,Uploaded_variation, .keep_all = TRUE)

#saving a copy of the VEP annotations selected by the testset selected gene
write.table(VEP_select_mcr5, "temp/VEPsngannot_dis579.txt", sep = "\t", col.names = TRUE, quote = FALSE)

#--------------------------------------------------------------
#select vep columns to add to annotation base file and add a header distinguishing

#rename colnames so that clear sourced from closest gene
colnames(VEP_select_mcr5) = paste("VEP", colnames(VEP_select_mcr5),sep="_")

#add columns to annotation
cisregANNOT2 = dplyr::left_join(cisregANNOT1, VEP_select_mcr5, by = c("ID" = "VEP_Uploaded_variation"), keep = TRUE)

#add columns that indicate VEP selected consequence is splicing or not splicing
# if string contains 'splic*' "VEPc_splice", "NA"
cisregANNOT2 = cisregANNOT2 %>% mutate(VEP_splice_con = if_else(str_detect(VEP_Consequence, "splic"), "vep_splice", "No"))
cisregANNOT2$ID %>% as.factor() %>% n_distinct()
cisregANNOT2$VEP_splice_con %>% as.factor() %>% summary()

###---VEP-closest-MANE-transcript-annotation----------------------------------------------------------

#selecting the annotation by closest MANE transcripts, first replace the "-" which occur when 'within' a gene to 0 
VEP_c0 = VEPfull
VEP_c0$DISTANCE = str_replace(VEP_c0$DISTANCE, "-", "0")
VEP_c0$DISTANCE = as.numeric(VEP_c0$DISTANCE)
# select the row which has the lowest number for each variant
VEP_closest = VEP_c0 %>%
  filter(SYMBOL_SOURCE == "HGNC" & MANE_SELECT != "-") %>%
  group_by(Uploaded_variation) %>% 
  slice_min(DISTANCE)  %>%
  ungroup()

VEP_closest_dup = VEP_closest %>% group_by(Uploaded_variation) %>% filter(n() >1)
VEP_closest = VEP_closest[!duplicated(VEP_closest$Uploaded_variation), ]
VEP_closest$Uploaded_variation %>% n_distinct()

#graph the distance to 'feature' in full vs 'closest selected' dataset
hist(as.numeric(VEP_closest$DISTANCE))

#select columns to add to annotation base file and add a header distinguishing
colnames(VEP_closest)
cols2b = c("Uploaded_variation","Consequence","IMPACT","SYMBOL","Gene","STRAND", "DISTANCE")
VEP_cols_close = VEP_closest[ ,c(cols2b)]
summary(VEP_cols_close)

#rename colnames so that clear sourced from closest gene
colnames(VEP_cols_close) = paste("VEPclose", colnames(VEP_cols_close),sep="_")

#optional add columns to annotation
cisregANNOT3 = dplyr::left_join(cisregANNOT2, VEP_cols_close, by = c("ID" = "VEPclose_Uploaded_variation"), keep = TRUE)

write.table(cisregANNOT3, "temp/ANNOT3_dis579.txt", sep = "\t", col.names = TRUE, quote = FALSE)


#--ANNOTATION-4----CADD-------##########################################

#CADD annotation addition - using grc38 v1.6 collected via INPUT.R
CADD = read.table("Data/dis579_CADD_GRCh38-v1.6_anno_f2c052d18d502f786e25fd3c720767cb.tsv", sep = "\t", header= TRUE)
colnames(CADD)

#make a variant key (to match to ANNOT set as not variant ID in CADD annotation)
CADD$location38 = paste(CADD$Chrom, CADD$Pos, sep = ":")
CADD$variant38ref = paste(CADD$location38, CADD$Ref, sep = "")
CADD$variant38 = paste(CADD$variant38ref, CADD$Alt, sep = ">")
CADD$variant38 %>% as.factor() %>% n_distinct()

CADD_unscored= subset(cisregANNOT3, !(variant38 %in% CADD$variant38))

#check
CADD$variant38 %>% n_distinct()

#first select the columns to add
CADD_cols = unique(CADD[ ,c("variant38","minDistTSS","minDistTSE","GC","CpG",
                            "tOverlapMotifs","EnsembleRegulatoryFeature","RawScore","PHRED")])

#adding datasource to dataframe header
colnames(CADD_cols) = paste("CADD", colnames(CADD_cols),sep="_")
cisregANNOT4 = left_join(cisregANNOT3, CADD_cols, by = c("variant38" = "CADD_variant38"))

cisregANNOT4$ID %>% n_distinct()

write.table(cisregANNOT4, "temp/ANNOT4_dis580.txt", sep = "\t", col.names = TRUE, quote = FALSE)

#########################################
#-final-analysis--------------------------
#profile of variants in set
cisregANNOT_final = cisregANNOT4

#number of SNVs in full dataset
sum(cisregANNOT_final$allelesize==0)

#number of indels in full dataset
sum(cisregANNOT_final$allelesize!=0)

#count number with literature
summary(as.factor(cisregDraw$literature_evidence))

#graph summary of the consequence annotation #########
#setting graphing settings
theme_set(theme_gray(base_size = 16))

# looking at the consequences of the variants - based on the selected transcript
cisregANNOT_final$VEP_Consequence = as.factor(cisregANNOT_final$VEP_Consequence)
cisregANNOT_final$VEP_Consequence %>% n_distinct()

ggplot(cisregANNOT_final, aes(x = VEP_Consequence)) + geom_bar() + 
  theme(axis.text.x = element_text(angle=45, size = 11))+ labs(title = "VEP_Consequence")

#adding a column with splice consequences grouped
cisregANNOT_final$VEP_Consequence = as.character(cisregANNOT_final$VEP_Consequence)
summary(as.factor(str_detect(cisregANNOT_final$VEP_Consequence, "splice")))

cisregANNOT_final$VEP_Consequence_adj = ifelse((str_detect(cisregANNOT_final$VEP_Consequence, "splice"))== TRUE, "splicing", cisregANNOT_final$VEP_Consequence)
summary(as.factor(cisregANNOT_final$VEP_Consequence_adj))

cisregANNOT_final$VEP_Consequence_adj = as.factor(cisregANNOT_final$VEP_Consequence_adj)

ggplot(cisregANNOT_final, aes(x = VEP_Consequence_adj)) + geom_bar() + 
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  labs(x = "",  y = "# variants")+
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.3, nudge_y = 1)

#-export-graphs-------------------------------------------------------------------------------------

svg("temp/dis580_trasnoV.svg")
ggplot(VEPrs_trans, aes(n)) +
  geom_histogram(binwidth = 1)+ 
  labs(x="# transcripts",  y = "# variants")
dev.off()

svg("temp/dis580_VEP_Consequence.svg")
ggplot(cisregANNOT_final, aes(x = VEP_Consequence_adj)) + geom_bar() + 
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  labs(x = "",  y = "# variants")+
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.3, nudge_y = 1)
dev.off()

ggplot(VEPrs_trans, aes(n)) + 
  geom_histogram()+ 
  labs(x="# transcripts",  y = "# variants")

#-----------------------------------------------------------------
#figure generation

transct = ggplot(VEPrs_trans, aes(n)) +
  geom_histogram(binwidth = 1)+ 
  labs(x="# transcripts",  y = "# variants")+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value= "darkslategray")

consplot = ggplot(cisregANNOT_final, aes(x = VEP_Consequence_adj)) + geom_bar() + 
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  labs(x = "",  y = "# variants")+
  ylim(0,350)+
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.3, nudge_y = 1)+
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value= "darkslategray")

plot_grid(transct, consplot, align = "h", axis = "b", rel_widths = c(1, 1.3),labels = c("A", "B"))

ggsave(path = "D:/cisreg_manuscript/figures", filename = "Fig_dis580.png", width = 9, height = 4, device='png', dpi=600)
