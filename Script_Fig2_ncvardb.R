##cis reg variant - ncdbvar analysis
#R code - using data from https://github.com/Gardner-BinfLab/ncVarDB/tree/master/data
#Version 1 - 21.09.23

#set up
library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(cowplot)
library(grid)

#--------------------------------------------
#ANNOTATION FILE BASE - INPUT
#SET UP FILEPATH OR upload the 'standard data for annotation'
#setwd 

#-----------------------------------------------------------------
#####################################################################
#Load starting dataset - will be translated ito a vcf like format
ncdbvar_path = read.delim("Data/ncVar_pathogenic.tsv", sep = "\t", header= TRUE)
ncdbvar_ben = read.delim("Data/ncVar_benign.tsv", sep = "\t", header= TRUE)
colnames(ncdbvar_ben)
colnames(ncdbvar_path)

#label with a column for exp.goup
ncdbvar_ben$exp.group = "benign"
ncdbvar_path$exp.group = "pathogenic"

#-----------------------------------------------------------------
#####################################################################
#Load starting dataset - will be translated ito a vcf like format

#ensure annotation base file in vcf format grc38 CHROM, LOC, ID (snpID or dataset ID), REF ( allele vcf format AG>A), ALT, QUAL, FILT and INFO columns
#select columns
ncdbvar_b = ncdbvar_ben[ ,c("chr","pos","ID","ref","alt","exp.group")]
ncdbvar_p = ncdbvar_path[ ,c("chr","pos","ID","ref","alt","exp.group")]
ncdbvar = rbind(ncdbvar_b,ncdbvar_p)

#create and fill vcf columns if not present 
ncdbvar = ncdbvar %>% 
  add_column(QUAL = ".") %>%
  add_column(FILTER = ".") %>%
  relocate(exp.group, .after = last_col())

#rename colnames if required
colnames(ncdbvar) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
head(ncdbvar)
summary(as.factor(ncdbvar$INFO))
summary(as.factor(ncdbvar$ID))

#replace blank cells with a ".", for e.g. in snpID with .
#ncdbvar$ID[ncdbvar$ID==""] <- NA
#ncdbvar$ID[is.na(ncdbvar$ID)] <- "."

#remove chr from CHROM column - as per vcf 
#ncdbvar$CHROM = gsub("chr","",ncdbvar$CHROM)

#sort by chromosome and then location
ncdbvar = arrange(ncdbvar, CHROM, POS)

# save a copy of the vcf file. 
write.table(ncdbvar, "temp/ncdbvar_VCF.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# add the header to align with required vcf file format
##fileformat=VCFv4.2
##fileDate= <date>
##source= .../Data/ncvardb/Script_Fig2_ncdbvar.R
##reference=ncdbvar
#

#---------------------------------------------------
#adding a variant key column into ncdbvar set to support joining and graphing 
colnames(ncdbvar)
ncdbvar$location38 = paste(ncdbvar$CHROM, ncdbvar$POS, sep = ":")
ncdbvar$variant38ref = paste(ncdbvar$location38, ncdbvar$REF, sep = "")
ncdbvar$variant38 = paste(ncdbvar$variant38ref, ncdbvar$ALT, sep = ">")

###################################################---------------------------------------
#ANNOTATIONS - 2
#VEP annotation - currently input above file into web GUI, transcript database Ensembl, don't include protein annotation, phenotype yes to all, add spliceAI and conservation - currently only AA and BLUOSUM
#----------------------------------------

##VEP annotation (includes) ANNOTATION created a vcf file as input for online GUI VEP see INPUT script, 
# TEMPLATE add date and details of most VEP annotation - 23.03.23 full set of ~13K
#remove # at start of VEP annotation to read
VEPfull = read.delim("Data/ncvardb_VEP_M9UhtJ4urDJmLWPj.txt", sep = "\t", header= TRUE)
colnames(VEPfull)

#formatting VEP to support downstream process - if required- data has # direct from vep web gui, parses in via delm as X. can remove manual and import at dt also
#colnames(VEPfull)[which(names(VEPfull) == "X.Uploaded_variation")] <- "Uploaded_variation"

#count the number of unique variants based on project identified "CR#"
VEPfull$Uploaded_variation %>% as.factor() %>% n_distinct()
VEPfull$MANE_SELECT %>% as.factor() %>% n_distinct()
VEPfull$Gene %>% as.factor() %>% n_distinct()

#count the number in each sample group
ncPath = ncdbvar_path$ID
VEPfull$exp.group = ifelse(VEPfull$Uploaded_variation %in% ncPath, "Path", "Ben")
summary(VEPfull$exp.group)
VEPfull %>% group_by(exp.group) %>% summarise(n_distinct(Uploaded_variation))

#adding a column with splice consequences grouped
VEPfull$Consequence = as.character(VEPfull$Consequence)
summary(as.factor(str_detect(VEPfull$Consequence, "splice")))

VEPfull$VEP_Consequence_adj = ifelse((str_detect(VEPfull$Consequence, "splice"))== TRUE, "splicing", VEPfull$Consequence)
summary(as.factor(VEPfull$VEP_Consequence_adj))

VEPfull$VEP_Consequence_adj = as.factor(VEPfull$VEP_Consequence_adj)
summary(VEPfull$exp.group)
VEPfull %>% group_by(exp.group) %>% summarise(n_distinct(Uploaded_variation))


#add column grouping consequences outside of standard into 1 group
#separate into 'main consequences' and 'other to combine'
VEP1 = VEPfull %>% filter(VEP_Consequence_adj == "splicing"| 
                                          VEP_Consequence_adj == "3_prime_UTR_variant" | 
                                          VEP_Consequence_adj == "5_prime_UTR_variant" |
                                          VEP_Consequence_adj == "downstream_gene_variant"|
                                          VEP_Consequence_adj == "intron_variant"|
                                          VEP_Consequence_adj == "missense_variant"|
                                          VEP_Consequence_adj == "upstream_gene_variant"|
                                          VEP_Consequence_adj == "NA"
)
summary(VEP1$VEP_Consequence_adj)
VEP1$VEP_Consequence_adj2 = paste(VEP1$VEP_Consequence_adj)
summary(as.factor(VEP1$VEP_Consequence_adj2))

VEP2 = subset(VEPfull, !(Uploaded_variation %in% VEP1$ID))
summary(VEP2$VEP_Consequence_adj)
VEP2$VEP_Consequence_adj2 = "other"
summary(VEP2$VEP_Consequence_adj2)

VEP = rbind(VEP1, VEP2)
VEP$VEP_Consequence_adj2 = as.factor(VEP$VEP_Consequence_adj2)
summary(VEP)
VEP %>% group_by(exp.group) %>% summarise(n_distinct(Uploaded_variation))

#select path variants
VEP_path = subset(VEP, (Uploaded_variation %in% ncPath))
VEP_path$Uploaded_variation %>% as.factor() %>% n_distinct()

#looking at MANE transcript annotations
#select MANEselect transcripts only
VEP_m1 = VEP[VEP$MANE_SELECT != "-", ]
VEP_m2 = VEP[VEP$MANE_PLUS_CLINICAL != "-", ]
VEP_m = rbind(VEP_m1, VEP_m2)
VEP_m_u = distinct(VEP_m, Uploaded_variation, .keep_all = TRUE)

VEP_m$Uploaded_variation %>% n_distinct()
VEP_m_u$Uploaded_variation %>% n_distinct()
VEP_m %>% group_by(exp.group) %>% summarise(n_distinct(Uploaded_variation))
VEP_m$SYMBOL %>% as.factor() %>% n_distinct()

VEP_m$Consequence %>% as.factor() %>% summary()

#looking at canonical transcript annotations
#select canonical transcripts only
VEP_c = VEP[VEP$CANONICAL == "YES", ]
VEP_c_u = distinct(VEP_c, Uploaded_variation, .keep_all = TRUE)
VEP_c$Uploaded_variation %>% n_distinct()
VEP_c %>% group_by(exp.group) %>% summarise(n_distinct(Uploaded_variation))
VEP_c$SYMBOL %>% as.factor() %>% n_distinct()

VEP_c$Consequence %>% as.factor() %>% summary()
VEP_c$VEP_Consequence_adj %>% as.factor() %>% summary()
VEP$VEP_Consequence_adj2 %>% as.factor() %>% summary()

#choosing to proceed with VEP canonical transcript annotations
VEP_cols = VEP_c_u %>% dplyr::select(
  "Uploaded_variation","Location","Allele","Consequence" ,"VEP_Consequence_adj", "VEP_Consequence_adj2", "exp.group") 
VEP_cols$Uploaded_variation %>% as.factor() %>% n_distinct()

ncdbvar_VEP = left_join(ncdbvar,VEP_cols, join_by(ID ==Uploaded_variation))

#---------------------------------------------
#joining in the CADD annotation - CADD web GUI v1.6 information colleted 3.6.23
#using the CADD annotation to look at some of the dataset features
CADD = read.table("Data/ncvardb_CADD_GRCh38-v1.6_anno_f96061336ede37af302eb852e8210d10.tsv", sep = "\t", header= TRUE)
colnames(CADD)

#ANNOTATION location38 and variant38 - columns giving location and variant specific values
CADD$location38 = paste(CADD$Chrom, CADD$Pos, sep = ":")
CADD$variant38ref = paste(CADD$location38, CADD$Ref, sep = "")
CADD$variant38 = paste(CADD$variant38ref, CADD$Alt, sep = ">")

CADD_cols = CADD %>% dplyr::select(
  "Chrom","Pos","Ref","Alt","Type","Length","RawScore", "PHRED" ,"location38" , "variant38") 
CADD_cols = unique(CADD_cols)
CADD_cols$variant38 %>% as.factor() %>% n_distinct()

ncdbvar_cadd = left_join(ncdbvar_VEP, CADD_cols, by = "variant38")
ncdbvar_cadd %>% group_by(variant38) %>% summarise(n_distinct(variant38))
summary(as.factor(ncdbvar_cadd$INFO))
ncdbvar_cadd$PHRED = as.numeric(ncdbvar_cadd$PHRED)
summary(ncdbvar_cadd$PHRED)

#graph CADD ncdbvar_cadd benign versus pathogenic PHRED - pejaver indicated
ggplot(ncdbvar_cadd, aes(PHRED, fill = INFO)) + geom_boxplot()

#exporting the graph
svg("231016_ncdbvar_CADDdensity2.svg")
ggplot(ncdbvar_cadd, aes(PHRED, fill = INFO)) + 
  geom_density(position = "identity", alpha = 0.2) + 
  geom_vline(xintercept = 22.7) + geom_vline(xintercept = 25.3)+
  theme(legend.position="none", text = element_text(size= 16)) +
  labs(x = "CADD",  y = "density")
dev.off()

colnames(ncdbvar_cadd)
ncdbvar_cadd %>% group_by(as.factor(INFO)) %>% summarise(mean = mean(PHRED),median = median(PHRED))

#looking at the number of variants that scores - total # scores, and # number unique variants scores (some scored twice based on alt genes, score same though)
ncdbvar_cadd %>% group_by(INFO) %>%  summarise()
colSums(!is.na(ncdbvar_cadd))
length(unique(ncdbvar_cadd$PHRED))

group_by(ncdbvar_cadd, INFO) %>%
  summarise(
    count = n(),
    mean = mean(PHRED, na.rm = TRUE),
    sd = sd(PHRED, na.rm = TRUE),
    median = median(PHRED, na.rm = TRUE),
    IQR = IQR(PHRED, na.rm = TRUE)
  )

#-------------------------------------------------------
#ncVarDB annotation - final analysis
#select
ncVarDBannot = ncdbvar_cadd
summary(ncVarDBannot$VEP_Consequence_adj)
ncVarDBannot$VEP_Consequence_adj2 = as.factor(ncVarDBannot$VEP_Consequence_adj2)
summary(ncVarDBannot$VEP_Consequence_adj2)


#analysing the results of the VEP - ensembl consequences
colnames(ncVarDBannot)
ncVarDBannot$PHRED %>% summary()
ncVarDBannot$VEP_Consequence_adj = as.factor(ncVarDBannot$VEP_Consequence_adj)
ncVarDBannot$exp.group = as.factor(ncVarDBannot$exp.group)
ncVarDBannot$PHRED = as.numeric(ncVarDBannot$PHRED)
ncVarDBannot$RawScore = as.numeric(ncVarDBannot$RawScore)

table(ncVarDBannot$VEP_Consequence_adj, ncVarDBannot$exp.group)
count = count(group_by(ncVarDBannot,exp.group),VEP_Consequence_adj)
ncVarDBannot$PHRED %>% summary()
count(group_by(ncVarDBannot,exp.group),VEP_Consequence_adj2)

ggplot(ncVarDBannot, aes(Consequence, PHRED, fill = exp.group))+ 
  geom_boxplot( alpha = 0.2) + theme(axis.text.x = element_text(angle=45), axis.text = element_text(size = 6))

ncConsplot = ncVarDBannot %>% filter(VEP_Consequence_adj2 != "NA")

svg("231017_ncVarDBannotvConsequence.svg")
ggplot(ncConsplot, aes(VEP_Consequence_adj2, PHRED, fill = exp.group))+ 
  geom_boxplot(alpha = 0.2) + theme(axis.text.x = element_text(angle=45))+  
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.title = element_blank(),text = element_text(size= 12)) +
  labs(x = "",  y = "CADD")

dev.off()

#-----------------------------------------------------------------
#figure generation

caddplot = ggplot(ncdbvar_cadd, aes(PHRED, fill = INFO)) + 
  geom_density(position = "identity", alpha = 0.2) + 
  geom_vline(xintercept = 22.7) + geom_vline(xintercept = 25.3)+
  theme(legend.position="none", text = element_text(size= 12)) +
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value= "darkslategray")+
  labs(x = "CADD",  y = "density")+
  annotation_custom(grobTree(textGrob("BP4", x=0.4,  y=0.9, hjust=0,gp=gpar(col="#666666", fontsize=10))))+
  annotation_custom(grobTree(textGrob("PP3", x=0.65,  y=0.9, hjust=0,gp=gpar(col="#666666", fontsize=10))))

ngroup = ncConsplot %>% group_by(VEP_Consequence_adj2) %>% count(exp.group)

consplot = ggplot(ncConsplot, aes(VEP_Consequence_adj2, PHRED, fill = exp.group))+ 
  geom_boxplot(alpha = 0.2) + theme(axis.text.x = element_text(angle=45))+  
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.title = element_blank(),text = element_text(size= 12)) +
  scale_fill_manual(values=c("#00BFC4","#F8766D","#9999CC"), na.value= "darkslategray")+
  labs(x = "",  y = "CADD")

plot_grid(caddplot, consplot, align = "h", axis = "b", rel_widths = c(1, 1.3),labels = c("A", "B"))
