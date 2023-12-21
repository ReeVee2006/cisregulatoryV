##cis reg variants - FILTER
##web based input files - annotation of variant list for filter
##base FILTERation datasets in data folder - 
##Creating input files with web/GUI based annotations 
##last executed 26.10.23 RV

library(tidyverse)
library(GenomicRanges)
#library(liftOver)

#set up files--------------------------------------------
#setwd 
setwd("D:/cisreg_manuscript/SUPP")

## Load datasets-Cisreg variants, control and disease to combine dataset and create into VCF file
#load disease-associated variant set
cisregDraw = read.table("Data/Supplemental_DisVar.txt", header= TRUE, sep = "\t")
#load control variant set
cisregCraw = read.table("Data/Supplemental_ContVar.txt", header= TRUE, sep = "\t")

#create list of genes represented in disease variant set
n_distinct(cisregDraw$Source_Reported_GeneName)
n_distinct(cisregDraw$Selected_gene_name)

#--base-file-annotation-tool-input-vcf---------------------------------------------------------------------
## create a combined dataset with C and D variants labelled in INFO column
#select Control vcf columns
colnames(cisregCraw)
cisregCraw
cisregCont = cisregCraw[ ,c("CHROM","POS","contsetid","REF","ALT","QUAL","FILTER")]
cisregCont = cisregCont %>% filter(FILTER == "PASS")
cisregCont$contsetid = paste("CRcont",cisregCont$contsetid, sep = "")
names(cisregCont)[names(cisregCont) == "contsetid"] <- "ID"
cisregCont$INFO = "Cont"

#select disease vcf columns, IDs with project identifiers CRD
cisregDis = cisregDraw[ ,c("Chr_GRCh38","Location_GRCh38","Testset_ID","Reference_Allele","Alternative_Allele")]
colnames(cisregDis) = c("CHROM", "POS", "ID", "REF", "ALT")

#add vcf columns to disease set
cisregDis$QUAL = "."
cisregDis$FILTER = "."
cisregDis$INFO = "Dis"

#add chr to the CHROM column of disease for vcf ,depending on input file required
#cisregDis$CHROM = paste("chr",cisregDis$CHROM,sep ="")
#remove chr from CHROM column - as per vcf, or depending on input file required
cisregCont$CHROM = str_replace(cisregCont$CHROM, "chr", "")

#combine cont and disease
cisregAll38 = rbind(cisregCont,cisregDis)

#sort by chromosome and then location
cisregAll38 = arrange(cisregAll38, CHROM, POS)
cisregAll38_vcf = cisregAll38
cisregAll38_vcf$QUAL = "."
cisregAll38_vcf$FILTER = "."
cisregAll38_vcf$INFO = "."

#--generate--vcf-for-annotation-tool-input---------------------------------
# save a copy of the vcf file.
write.table(cisregAll38_vcf, "temp/cisreg_all38_grc38_vcf.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# add the header to align with required vcf file format
##fileformat=VCFv4.2
##fileDate= <date>
##source script=ANNOT_FILTER.R
##reference=cisregAll38_vcf
#

#--building-base-file_ for-FILTER------------------------------------------
#add location and variant columns 
cisregAll38$location38 = paste(cisregAll38$CHROM, cisregAll38$POS, sep = ":")
cisregAll38$variant38ref = paste(cisregAll38$location38, cisregAll38$REF, sep = "")
cisregAll38$variant38 = paste(cisregAll38$variant38ref, cisregAll38$ALT, sep = ">")
cisregall38loc = cisregAll38$location38
cisregall38var = cisregAll38$variant38

#create a list of control and disease location and variants
cisreg_dis = cisregAll38[ cisregAll38$INFO =="Dis", ]
cisreg_cont = cisregAll38[cisregAll38$INFO =="Cont", ] 

#####annotation ROI onto base file -----------
#--------------annotate variants with the gene ROI
#ROI row by row creation of full dataset
#upload input file with ROI regions and associated gene name 
mane_crroi = read.table("Data/maneroi_annot.txt", header= TRUE, sep = "\t")

#create a column for the ROI region location
colnames(mane_crroi)
mane_crroi$loc38regionstart = paste(mane_crroi$Chromosome.scaffold.name, mane_crroi$start, sep = ":")
mane_crroi$loc38region = paste(mane_crroi$loc38regionstart, mane_crroi$end, sep = "-")

#list of genes represented in disease variant set= crdisgenes_all
#selecting the cisreg disease variant set ROIs
mane_crroi$crdis_gene = ifelse(mane_crroi$Gene.name %in% cisregDraw$ROI_gene_name, "crdis_gene", "NA")
mane_crroi$crdis_gene = as.factor(mane_crroi$crdis_gene) 
summary(mane_crroi$crdis_gene)

roi = mane_crroi[mane_crroi$crdis_gene == "crdis_gene", ]
row.names(roi) = 1:nrow(roi)
colnames(roi)
roitest_cols = roi[ ,c("Gene.name","loc38region","Strand","width")]

#creating ROI file for Disease variants
#match the disease variants ROI information by the 'selected gene'
cisregDraw_ROI = left_join(cisregDraw, roitest_cols, by = c("ROI_gene_name" = "Gene.name"), keep = TRUE)

#select columns to add to annotation base file and add a header distinguishing
ROI_cols_dis = unique(cisregDraw_ROI[, c("Testset_ID","Gene.name","loc38region","Strand","width")])
unique(ROI_cols_dis$Testset_ID) %>% n_distinct()
colnames(ROI_cols_dis) = c("ID","ROIgene","ROIloc","ROIstrand","ROIwidth")

#creating ROI file for Control variants
cisregCraw_ROI = cisregAll38[cisregAll38$INFO =="Cont", ]
colnames(cisregCraw_ROI)
#provide column names
columns = c(colnames(cisregCraw_ROI))
#create empty file testsetCV
cisregC_ROI = data.frame(matrix(nrow = 0,ncol = length(columns)))
colnames(cisregC_ROI) = columns

for (i in 1:nrow(roi)) {
  row = roi[i, ]
  testset = filter(cisregCraw_ROI, CHROM == roi[[i,1]], POS >= roi[[i,2]] , POS <= roi[[i,3]])
  { if (nrow(testset != 0)) {
    testset$ROIgene = roi[[i,"Gene.name"]]
    testset$ROIloc = roi[[i,"loc38region"]]
    testset$ROIstrand = roi[[i,"Strand"]]
    testset$ROIwidth = roi[[i,"width"]]
    cisregC_ROI = rbind(cisregC_ROI,testset)}}}

#check
cisregC_ROI$ID %>% n_distinct()
cisregC_ROI$ROIgene %>% as.factor() %>% n_distinct()

#select columns to add to annotation base file and add a header distinguishing
ROI_cols_cont = unique(cisregC_ROI[, c("ID","ROIgene","ROIloc","ROIstrand","ROIwidth")])
ROI_cols_cont %>% distinct("ID", .keep_all = T)
ROI_cols_cont = ROI_cols_cont[!duplicated(ROI_cols_cont$ID),]
unique(ROI_cols_cont) %>% n_distinct()

ROI_cols = distinct(rbind(ROI_cols_cont,ROI_cols_dis))
ROI_cols = ROI_cols[!duplicated(ROI_cols$ID),]
n_distinct(ROI_cols$ID)
subset(ROI_cols,duplicated(ID))

#add ROI columns to the base annotation file
cisregAll38 = left_join(cisregAll38, ROI_cols, by = "ID")
n_distinct(cisregAll38$ROIgene)
summary(as.factor(cisregAll38$ROIgene))

#add raw data information on selected gene to the cisregdataset (for later ROI filter)
colnames(cisregDraw)
cisregDraw_genecols = cisregDraw[ ,c("Testset_ID","Selected_gene_name")]

#adding datasource to dataframe header
cisregAll38 = left_join(cisregAll38, cisregDraw_genecols, by = c("ID" = "Testset_ID"))
n_distinct(cisregAll38$Selected_gene_name)
summary(as.factor(cisregAll38$Selected_gene_name))
n_distinct(cisregAll38$ROIgene, na.rm = TRUE)
summary(as.factor(cisregAll38$ROIgene))

#check the variants with no ROI gene allocated - verify no viable ROI 
cisregAll38_noroi = cisregAll38 %>% filter(is.na(ROIgene))
#one NA remains, this is as variant overlaps the TSS, variant will be filtered out when ROI filter applied

#write.table(cisregAll38, "Data/cisreg_all38_grc38.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#--load VEP annotation------------

#load and check VEPannotation
#nb. remove # at start of VEP annotation to read
cisregAll_VEP = read.delim("Data/OaVt9eEsWgLnckRd.txt", sep = "\t", header= TRUE)
#formatting VEP to support downstream process
cisregAll_VEP$Uploaded_variation %>% as.factor() %>% n_distinct()
colnames(cisregAll_VEP)

#annotate the VEP file with exp group - noting that need to use the original raw MAIN sheet variant format
disID = cisreg_dis$ID
cisregAll_VEP$INFO <- ifelse(cisregAll_VEP$Uploaded_variation %in% disID, "Dis", "Cont")
cisregAll_VEP$Uploaded_variation %>% n_distinct()
cisregAll_VEP$INFO %>% as.factor() %>% summary()
#select the unique FILTERations, based on variant
VEPvariants = unique(cisregAll_VEP$Uploaded_variation)
summary(VEPvariants)
cisregAll_VEPu = dplyr::distinct(cisregAll_VEP, Uploaded_variation, .keep_all = TRUE)
cisregAll_VEPu$INFO = as.factor(cisregAll_VEPu$INFO)
summary(cisregAll_VEPu$INFO)

#formatting-VEP-data-to-separate-location-for-overlap
cisregAll_VEP$Location2 = cisregAll_VEP$Location
cisregAll_VEP = separate(data = cisregAll_VEP, col = Location2, into = c("CHROM","POS38"), sep = ":")
cisregAll_VEP = separate(data = cisregAll_VEP, col = POS38, into = c("POS38","POS38b"), sep = "-")
cisregAll_VEP$Location2 = cisregAll_VEP$Location
cisregAll_VEP = separate(data = cisregAll_VEP, col = Location2, into = c("location38","POS38b"), sep = "-")

#------identifying variants in a genomic coding region----------------------------------------------------
#ie overlap between the genomic coding intervals and the variant list (including all genomic coding region - non coding intronic variants also)
#load the file of all gene widest genomic coding intervals, created from ensembl biomart 27.02.23
codingregionsmane = read.table("Data/230228_manecodingregions.txt", header = TRUE, sep = "\t")

#creating a bed/granges of the coding regions
codingregionsmane$seqnames = paste("chr",codingregionsmane$Chromosome.scaffold.name,sep ="")
codingregionsmane$grstrand = gsub("1","",as.character(codingregionsmane$Strand))
codingregionsmane["grstrand"][codingregionsmane["grstrand"] == ''] <- "+"

gr_codingregions =  makeGRangesFromDataFrame(codingregionsmane,
                                             keep.extra.columns=FALSE,
                                             ignore.strand=FALSE,
                                             strand.field = "grstrand",
                                             seqinfo=NULL,
                                             seqnames.field=c("seqnames"),
                                             start.field=c("Genomic.coding.start"),
                                             end.field=c("Genomic.coding.end"),
                                             starts.in.df.are.0based=FALSE)

#create GRanges object from cisreg  variant list
#reformat-cisreg-list-forgranges-input
cisregAll38$grchr_38 = paste("chr",cisregAll38$CHROM,sep="")
#calculating end position - support determination of 1based including all variant types
cisregAll38$ALTcount = nchar(cisregAll38$ALT)
cisregAll38$REFcount = nchar(cisregAll38$REF)
cisregAll38$allelesize = cisregAll38$ALTcount - cisregAll38$REFcount
cisregAll38$width = cisregAll38$allelesize +1
cisregAll38$end.field = cisregAll38$POS + cisregAll38$REFcount

#converting dataframe into GRanges object - from vcf like file
gr_cisregAll38 =  makeGRangesFromDataFrame(cisregAll38,
                                           keep.extra.columns=FALSE,
                                           ignore.strand=TRUE,
                                           seqinfo=NULL,
                                           seqnames.field=c("grchr_38"),
                                           start.field=c("POS"),
                                           end.field=c("end.field"),
                                           starts.in.df.are.0based=FALSE)

#overlap genomic coding regions and test variants to see if test set are located within the BED file
all38coding_gr = findOverlaps(gr_codingregions, gr_cisregAll38)
all38coding_df = as(all38coding_gr, "data.frame")

#overlap coding regions and cisregdisease-set to see if test set are located within the BED file
codregvar_gr = subsetByOverlaps(gr_cisregAll38, gr_codingregions, ignore.strand=TRUE)
codregvar_df = as(codregvar_gr, "data.frame")
head(codregvar_df)

#create a list of all coding region variants
codregvar_df$CHROM = gsub("chr","",codregvar_df$seqnames)
codregvar_df$loc38 = paste(codregvar_df$CHROM, codregvar_df$start, sep = ":")
codregvar_df$loc38 %>% as.factor() %>% n_distinct()
cisregAll38$location38 %>% as.factor() %>% n_distinct()
codingregionLoc = unique(codregvar_df$loc38)

#annotate variant file with coding region 
cisregAll38$codingregion <- as.factor(ifelse(cisregAll38$location38 %in% codingregionLoc, "codingregion", "no"))
summary(cisregAll38)
codingregions_df = cisregAll38 %>% filter(codingregion == "codingregion")
codingregions_df %>% group_by(INFO) %>% summarise(n_distinct(variant38))

#------identifying protein coding variants----------------------------------------------------------------

#create a data.table of variants in coding region, create a list of coding locations
#select VEP annotations that overlap with an amino acid - any transript
coding = cisregAll_VEP[cisregAll_VEP$Amino_acids != "-", ]
coding$Amino_acid = as.factor(coding$Amino_acid)
summary(coding$Amino_acid)
codingLoc = unique(coding$location38)
codingVar = unique(coding$Uploaded_variation)
coding$INFO = as.factor(coding$INFO)
summary(coding)

#summary of coding variants
summary(as.factor(coding$SYMBOL))
coding$SYMBOL %>% n_distinct()

#annotate variant file to identify if aa annotation overlap
cisregAll38$codingaa <- as.factor(ifelse(cisregAll38$ID %in% codingVar, "codingaa", "no"))
summary(cisregAll38)

#------identify gwas variants-----
###------exclude-Alsheikh-2022-validated-GWAS variants--

#load dataset of control variants listed in Alsheikh et al 2022  
Alsh22 = read.csv("Data/12920_2022_1216_MOESM3_ESM.csv", header = TRUE)

#create a list of locations from Alsh22 - should be 306
colnames(Alsh22)
Alsh22$location38 = paste(Alsh22$Chr,sep = ":", Alsh22$Loc..hg38.)
Alsh22$Loc38 = Alsh22$location38
Alsh22$Loc38 <-gsub("chr","",as.character(Alsh22$Loc38))

#identify overlap in publication locations with our functional disease set (rare disease)
Alsh22Loc38 = Alsh22$Loc38
intersect(Alsh22Loc38, cisregAll38$location38)
overlap_ncvAlsh22 = intersect(Alsh22Loc38, cisregall38loc)

#annotate variant file to identify if a vaidated GWAS variant overlaps with Alsheikh et al 2022 variant location
cisregAll38$Alsh22 <- as.factor(ifelse(cisregAll38$location38 %in% Alsh22Loc38, "GWAS", "no"))
summary(cisregAll38)

#------identify variants for literature exclusion---------------------------------
#identify all manually excluded variants (combine - yes=clinical record queried and yes pc=literature suggests protein coding)
Draw_manexclude = cisregDraw[cisregDraw$Literature_exclude == "yes"| cisregDraw$Literature_exclude == "yes pc", ]
manualexclude = Draw_manexclude$Testset_ID
manualexclude_var =  Draw_manexclude$GRCh38variant_B

#annotate variant file with manual exclusion tag
cisregAll38$manexclude <- as.factor(ifelse(cisregAll38$ID %in% manualexclude, "manexclude", "no"))
summary(cisregAll38)

#------identify common variants based on AF--------------------------------------------
#excluding common variants, AF based on VEP annotation
cisreg_All38_vepAF = cisregAll_VEP
colnames(cisreg_All38_vepAF)
summary(cisreg_All38_vepAF$gnomADg_AF)

#selecting highest alternative AF from (non-founder) populations, Non-Finnish European, South-Asian, African-American/African ancestry, Latino, East Asian populations. 
#create new column with the calculated maxAF 
pops = c("gnomADg_NFE_AF","gnomADg_AFR_AF","gnomADg_SAS_AF","gnomADg_EAS_AF","gnomADg_AMR_AF")
cisreg_All38_vepAF <- cisreg_All38_vepAF %>% mutate_at(pops, as.numeric)
cisreg_All38_vepAF$maxAF = apply(cisreg_All38_vepAF[,c(pops)],1,max)
summary(cisreg_All38_vepAF$maxAF)

#create a new column listing the population with the max AF
cisreg_All38_vepAF$maxAFpop = as.factor(colnames(cisreg_All38_vepAF[ ,pops])[max.col(cisreg_All38_vepAF[ ,pops])])
summary(cisreg_All38_vepAF$maxAFpop)

cisreg_All38_vepAF_cols = cisreg_All38_vepAF[ ,c( "Uploaded_variation","gnomADg_AF",pops, "maxAF", "maxAFpop")]
cisreg_All38_vepAF_cols = distinct(cisreg_All38_vepAF_cols)
colnames(cisreg_All38_vepAF_cols)
summary(cisreg_All38_vepAF_cols)

#select variants by ID with maxAF >0.01
cisreg_All38_vepAF_cols$maxAF = as.numeric(cisreg_All38_vepAF_cols$maxAF)
commonVar = cisreg_All38_vepAF_cols %>% filter(., maxAF >= 0.01)
commonID = commonVar$Uploaded_variation

#annotate variant file if commonAF or no
cisregAll38$commonAF <- as.factor(ifelse(cisregAll38$ID %in% commonID, "commonAF", "no"))
summary(cisregAll38)


#------identify conflicting ClinVar classifications---------------------

#investigate disease variants and how they align with Clinvar classification, as annotated by VEP
cisregAll_VEP$CLIN_SIG = as.character(cisregAll_VEP$CLIN_SIG)
summary(as.factor(cisregAll_VEP$CLIN_SIG))
cisregAll_VEP$anyben = str_detect(cisregAll_VEP$CLIN_SIG, "benign")
cisregAll_VEP$anypath = str_detect(cisregAll_VEP$CLIN_SIG, "pathogenic")
cisregAll_VEP$anyben %>% as.factor %>% summary()
cisregAll_VEP$anypath %>% as.factor %>% summary()

#select the disease set with benign Clinvar in VEPu FILTERation
cisreg_dis_VEP = cisregAll_VEP[cisregAll_VEP$INFO == "Dis", ]
summary(as.factor(cisreg_dis_VEP$anyben))
cisreg_dis_VEP_anyben = cisreg_dis_VEP[cisreg_dis_VEP$anyben == "TRUE", ]
cisreg_dis_VEP_anyben_var38 = str_unique(cisreg_dis_VEP_anyben$variant38_vcfB)
cisreg_dis_VEP_anyben_ID = str_unique(cisreg_dis_VEP_anyben$Uploaded_variation)

#graph and export "bargraph_cisregDis_benignclinvar"..., add date and save in relevant dataset folder
e = ggplot(cisregAll_VEP, aes(x=INFO, fill=anyben)) 
e + geom_bar(position = "fill") + labs(x = "Experimental group", y = "proportion contain benign")

#count of disease variants with a benign classification (of any sort) in clinvar
cisreg_dis_VEP_anyben$Uploaded_variation %>% n_distinct()

#investigate control variants with a pathogenic classification

#select the control set with a pathogenic clinvar FILTERation in  VEPu FILTERation file
cisreg_cont_VEP = cisregAll_VEP[cisregAll_VEP$INFO == "Cont", ]
summary(as.factor(cisreg_cont_VEP$anypath))
cisreg_cont_VEP_anypath = cisreg_cont_VEP[cisreg_cont_VEP$anypath == "TRUE", ]
cisreg_cont_VEP_anypath_var38 = str_unique(cisreg_cont_VEP_anypath$variant38_vcfB)
cisreg_cont_VEP_anypath_ID = str_unique(cisreg_cont_VEP_anypath$Uploaded_variation)

#graph and export "bargraph_cisregCont_pathclinvar"..., add date and save in relevant dataset folder
c = ggplot(cisregAll_VEP, aes(x=INFO, fill=anypath)) 
c + geom_bar(position = "fill") + labs(x = "Experimental group", y = "proportion pathogenic")

#count of control variants with a disease classification (of any sort) in clinvar
cisreg_cont_VEP_anypath$Uploaded_variation %>% n_distinct()

#annotate variant file with conflicting clinvar classification - disease and cont separate
cisregAll38$contPathinCVconflict <- as.factor(ifelse(cisregAll38$ID %in% cisreg_cont_VEP_anypath_ID, "ContPathinCV", "no"))
cisregAll38$disBeninCVconflict <- as.factor(ifelse(cisregAll38$ID %in% cisreg_dis_VEP_anyben_ID, "DisBeninCV", "no"))
summary(cisregAll38$contPathinCVconflict)
summary(cisregAll38$disBeninCVconflict)

#------identify Splicing-annotation---------------------------------------------
#load spliceAI results sourced 22.04.2023 through  HPC cmd line spliceAI tool
splice_raw = read.table("Data/all38_scored_variants_annotatedv2_220423.tsv", sep = "\t", header= TRUE, stringsAsFactors = TRUE)
summary(splice_raw)

splice = tidyr::separate(splice_raw, col = ID, into = c("Testset_ID", "exp.group"), sep = "_", remove = FALSE )
splice$exp.group = str_replace(splice$exp.group, "ExpG=", "")
splice$Testset_ID = as.factor(splice$Testset_ID)
splice$exp.group = as.factor(splice$exp.group)

#how many variants have FILTERated with splice AI
splice$Testset_ID %>% as.factor() %>% n_distinct()

#how many variants have FILTERated with a score indicating they are 'potential splice variants'
#calculate a max spliceAI score, create new column, which includes the max score from the selected relevant splice AI columns
#select the max value from a range of columns.
colnames(splice)
spcols = c("DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")
splice = splice %>% rowwise() %>% mutate(max_spAI = max(DS_AG,DS_AL,DS_DG,DS_DL))
summary(splice$max_spAI)

#select the variants that have a maxAI greater than 0.2
splice_impact02 = splice %>% filter(max_spAI >= 0.2)
splice_impact02$Testset_ID %>% as.factor() %>% n_distinct()
splice_impact02_IDs = unique(splice_impact02$Testset_ID)

#annotation of the splicing file with yes or no, for potential splicing impact. 
splice$spliceAI_impact02 <- ifelse(splice$Testset_ID %in% splice_impact02_IDs, "splicing", "NA")
as.factor(splice$spliceAI_impact02) %>% summary()

#looking a little at what the splice impacts are - starting with the 'potential clinical splicing impact >0.2)
summary(splice_impact02)

cisregAll38$spliceAI_impact02 <- as.factor(ifelse(cisregAll38$ID %in% splice_impact02_IDs, "PredSplicing", "no"))
summary(cisregAll38$spliceAI_impact02)

##select the final reference set###################################
#select the final reference set- disease first
colnames(cisregAll38)
summary(cisregAll38)
#identify all the filter exclusion annotations
filters = c("codingregion","codingaa","Alsh22","manexclude","commonAF","contPathinCVconflict","disBeninCVconflict","spliceAI_impact02")
#select reference set - ie remove variants that with an exclusion annotation
cisregall38_filter = cisregAll38 %>% filter(codingregion == "no",
  codingaa == "no",
  Alsh22 == "no",
  manexclude == "no",
  commonAF== "no",
  contPathinCVconflict == "no",
  disBeninCVconflict == "no",
  spliceAI_impact02 == "no"
)
#select the disease variants that remain in reference set
cisregall38_dis_filter = cisregall38_filter %>% filter(INFO == "Dis")
summary(cisregall38_filter)
cisregall38_dis_filter$INFO %>% as.factor() %>% summary()

#identify the genes within the reference set
refset_genes = cisregall38_dis_filter$Selected_gene_name

#select the control variants that remain in reference set
cisregall38_cont_filter = cisregall38_filter %>% filter(INFO == "Cont")
summary(cisregall38_cont_filter)
cisregall38_cont_filter$INFO %>% as.factor() %>% summary()

#filter the control variants so any variants in filtered out ROIs are removed
cisregall38_cont_filter_ROIrem= subset(cisregall38_cont_filter, (ROIgene %in% c(refset_genes)))

#join the filtered control and disease variants - final reference set
refset = rbind(cisregall38_cont_filter_ROIrem, cisregall38_dis_filter)
refset = arrange(refset, CHROM, POS)

#annotate the all38 set with in and out reference set
cisregAll38$refsetROI <- as.factor(ifelse(cisregAll38$ROIgene %in% refset_genes, "ROIgene", "no"))
summary(cisregAll38$refsetROI)
cisregAll38$refsetV <- as.factor(ifelse(cisregAll38$ID %in% refset$ID, "refset", "no"))
summary(cisregAll38$refsetV)

##final save #################################

#save FILTER annotation - last executed 26.10.23
write.table(cisregAll38, "temp/cisregall38_FILTER.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Supplemental Table 8

colnames(cisregAll38)
supp8 = cisregAll38 %>% select("CHROM","POS","ID","REF","ALT","INFO",
                                      "ROIgene","ROIloc","ROIstrand", "ROIwidth",
                                      "Selected_gene_name", 23:32)
names(supp8)[names(supp8) == "INFO"] = "Group"

write.table(supp8, "Data/cisreg_allgrc38_supp8.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


#check and save final FILTER annotation
summary(refset)
refset$INFO %>% as.factor() %>% summary()
refset$ROIgene %>% as.factor() %>% n_distinct()
refset$variant38 %>% as.factor() %>% n_distinct()
refset$ID %>% as.factor() %>% n_distinct()

#save reference set - last executed 26.10.23
write.table(refset, "temp/cisreg_refset_grc38.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#--generate--vcf-for-annotation-tool-input---------------------------------
#select columns
refset_vcf = refset[ ,c("CHROM","POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")]

#rename colnames if required
#colnames(refset_vcf) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

#create and fill vcf columns if not present 
#refset_vcf$FILT = "."

#replace blank cells with a ".", for e.g. in snpID with .
#refset_vcf$ID[refset_vcf$ID==""] <- NA
#refset_vcf$ID[is.na(refset_vcf$ID)] <- "."

#sort by chromosome and then location
refset_vcf = arrange(refset_vcf, CHROM, POS)

# save a copy of the vcf file - last execute to save 25.10.23
write.table(refset_vcf, "temp/cisreg_refset_grc38_vcf.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# add the header to align with required vcf file format
##fileformat=VCFv4.2
##fileDate= <date>
##source script=ANNOT_FILTER.R
##reference=cisreg_refset_vcf
#

#####################################

