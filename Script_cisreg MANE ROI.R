##cis reg variant - control selection, from MANE start sites
##R code complemented by base annotation datasets in corresponding data folder - including source files for project and sources files for annotation, create a folder for the specific annotation set which is directd within specific folder in code 
## Creating input ROI files for web/GUI based annotation - last update 14.04.23 RV


#install required packages
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
#BiocManager::install("liftOver")
library(tidyverse)
library(GenomicRanges)
#library(liftOver)

#--------------------------------------------
#setwd 
setwd("D:/cisreg_manuscript/SUPP")

## Load datasets to annotate. Cisreg variants, control and disease to combine dataset and create into VCF file
#load disease-associated variant set
cisregDraw = read.table("Data/Supplemental_Table3_disvar.txt", header= TRUE, sep = "\t")
#load control variant set
cisregCraw = read.table("Data/Supplemental_Table6_contvar.txt", header= TRUE, sep = "\t")
#create list of genes represented in disease variant set
crdis_genes = unique(cisregDraw$Source_Reported_GeneName)

#---------------------------------------------------------------------------------------------

#load ensemble biomart dump of MANE gene features. 
MANE = read.delim("Data/martquery_0228035716_130.txt", header= TRUE, sep = "\t")
colnames(MANE)

#creating a MANE gene coding regions bed
#select the chr, genomic start and genomic end of coding. create a bed file
start = aggregate(Genomic.coding.start ~ Gene.stable.ID.version, MANE, function(x) min(x))
end = aggregate(Genomic.coding.end ~ Gene.stable.ID.version, MANE, function(x) max(x))
regions = left_join(start,end, by = "Gene.stable.ID.version")

#add mane gene info
maneinfo = MANE %>% dplyr::select(c("Gene.stable.ID.version","Gene.name","Chromosome.scaffold.name","Strand",
                                           "Transcription.start.site..TSS.","Transcript.end..bp.","Transcript.count")) %>% unique()
                                           
maneregions = left_join(regions,maneinfo, by = "Gene.stable.ID.version")

#sort by chromosome and then location
maneregions = arrange(maneregions, Chromosome.scaffold.name, Genomic.coding.start)

write.table(maneregions, "temp/manegregions.txt",col.names= TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#-------------------------------------------------------------
#creating a granges file for mane gene coding regions
colnames(maneregions)
maneregions$Strand = gsub( "1","", maneregions$Strand)
maneregions$Strand[maneregions$Strand== ""] = "+"
maneregions$width = maneregions$Genomic.coding.end - maneregions$Genomic.coding.start
gr_roi_testgenes =  makeGRangesFromDataFrame(maneregions,
                                             keep.extra.columns=FALSE,
                                             ignore.strand=FALSE,
                                             strand.field = "Strand",
                                             seqinfo=NULL,
                                             seqnames.field=c("Chromosome.scaffold.name"),
                                             start.field=c("Genomic.coding.start"),
                                             end.field=c("Genomic.coding.end"),
                                             starts.in.df.are.0based=TRUE)
gr_roi_testgenes$gene.name = maneregions$Gene.name

####################################################
#working on the cr disease gene set
#######################################################
#--------------------------------------------------------------
#selecting the disease set variant gene regions, mane genes. 
colnames(maneregions)
maneregions_crdis = subset(maneregions, (Gene.name %in% c(crdis_genes)))
maneregions_crdis$Gene.stable.ID.version %>% as.factor() %>% n_distinct()
manegenes = maneregions_crdis$Gene.name
setdiff(crdis_genes, manegenes)

#----------------------------------------------------------------
#creating a cisreg region ROI dataframe - all genes
#separate into strands, two sep files for creating a MANE gene upstream ROI regions bed,
mane_pos = maneregions[maneregions$Strand =="+", ]
mane_neg = maneregions[maneregions$Strand =="-", ]

#postive strand - new start column -5kb from TranscriptSS
#pos strand - new end column -5kb from TranscriptSS coding start 
mane_pos$start = as.numeric(mane_pos$Transcription.start.site..TSS.)-5000
mane_pos$end = mane_pos$Genomic.coding.start

#negative strand - new start column coding end site
mane_neg$start = mane_neg$Genomic.coding.end
#neg strand - new end column +5kb from Transcript end site
mane_neg$end = as.numeric(mane_neg$Transcript.end..bp.)+5000

#join pos and negative - select pos chr, start, end and 
posroi = mane_pos %>% dplyr::select(c("Chromosome.scaffold.name","start", "end", "Gene.name","Strand"))
#select neg chr, start end
negroi = mane_neg %>% dplyr::select(c("Chromosome.scaffold.name","start", "end", "Gene.name", "Strand"))
#bind dataframes and arrange
mane_crroi = rbind(posroi,negroi)
mane_crroi = arrange(mane_crroi, Chromosome.scaffold.name, start)

write.table(mane_crroi, "temp/maneroi_annot", quote = FALSE, sep = "\t", col.names = TRUE)
#save table for use in annotating variants with roi
write.table(mane_crroi, "Data/maneroi_annot", quote = FALSE, sep = "\t", col.names = TRUE)

#select the TSS and new end site for both positive and negative strands, then join file so pos and neg in one datatabase, covert to bed ROI file
#check that all our cisreg disease genes are in the MANE gene list
#-------------------------------------------------------------
#identifying overlapping regions of interest across all mane genes
#creating a granges file for mane gene promoter regions of interest
#formatting file to parse to Granges
colnames(mane_crroi)
mane_crroi$width = mane_crroi$end - mane_crroi$start

#creating a granges file
gr_roi_manegenes =  makeGRangesFromDataFrame(mane_crroi,keep.extra.columns=FALSE,
  ignore.strand=FALSE,
  strand.field = "Strand",
  seqinfo=NULL,
  seqnames.field=c("Chromosome.scaffold.name"),
  start.field=c("start"),
  end.field=c("end"),
  starts.in.df.are.0based=TRUE)

#save a bed file of all mane gene reg ROI
export(gr_roi_manegenes, "temp/maneregroi.bed")

#identifying overlapping regions
reduce(gr_roi_manegenes)
GenomicRanges::reduce(do.call(c, gr_roi_manegenes))
redroimane = as.data.frame(d)

#------------------------------------------------------------

#selecting the test set variant gene regions, mane genes. 
colnames(maneregions)
mane_crroi_testgenes = subset(mane_crroi, (Gene.name %in% c(crdis_genes)))
mane_crroi_testgenes$Gene.name %>% as.factor() %>% n_distinct()

#----------------------------------------------------------------
#converting dataframe into GRanges object - from vcf like file
colnames(mane_crroi_testgenes)
mane_crroi_testgenes$Strand = gsub( "1","", mane_crroi_testgenes$Strand)
mane_crroi_testgenes$Strand[mane_crroi_testgenes$Strand== ""] = "+"
mane_crroi_testgenes$width = mane_crroi_testgenes$end - mane_crroi_testgenes$start

write.table(mane_crroi_testgenes, "temp/maneroi_crdis__annot", quote = FALSE, sep = "\t", col.names = TRUE)

#------------------------------------------------------------
#looking at the size of the regions
colnames(mane_crroi_testgenes)
a = ggplot(mane_crroi_testgenes, aes(x = reorder(Gene.name, width), y = width))
a + geom_col() + labs(x = "Gene", y = "ROI width")


#-------------------------------------------------------------
#creating a granges file from the crdis gene region of interest
gr_roi_testgenes =  makeGRangesFromDataFrame(mane_crroi_testgenes,
                                                   keep.extra.columns=TRUE,
                                                   ignore.strand=FALSE,
                                                   strand.field = "Strand",
                                                   seqinfo=NULL,
                                                   seqnames.field=c("Chromosome.scaffold.name"),
                                                   start.field=c("start"),
                                                   end.field=c("end"),
                                                   starts.in.df.are.0based=TRUE)

#------------------------------------------------------------------------
#creating a granges file from the crdis gene region of interest
gr_roi_testgenes =  makeGRangesFromDataFrame(mane_crroi_testgenes,
                                             keep.extra.columns=TRUE,
                                             ignore.strand=FALSE,
                                             strand.field = "Strand",
                                             seqinfo=NULL,
                                             seqnames.field=c("Chromosome.scaffold.name"),
                                             start.field=c("start"),
                                             end.field=c("end"),
                                             starts.in.df.are.0based=TRUE)

granges(gr_roi_testgenes)
reduced_gr = GenomicRanges::reduce(gr_roi_testgenes, min.gapwidth=1L, ignore.strand=TRUE)
length(gr_roi_testgenes)
length(reduced_gr)

#convert to a bed file

export(reduced_gr, "temp/cisregROI_crdisgenes.txt")
export.bed(reduced_gr, "temp/cisregROI_crdisgenes.bed",format = "bed")

ts23_roi = as.data.frame(reduced_gr)
ts23_roi_bed = ts23_roi %>% dplyr::select(c("seqnames","start","end"))
ts23_roi_bed = arrange(ts23_roi_bed, seqnames, start, end)
export.bed(ts23_roi_bed, con = "temp/crdis23_ROI.bed")
ts23_roi_bed$seqnames = paste0("chr",ts23_roi_bed$seqnames,"")
write.table(ts23_roi_bed, "temp/crdis23_ROI_2.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#------------------------------------------------------------------------------------------
#creating a list of the genes, selected refseq transcript IDs
#gene names - 191 in crdis
ENSgenesymbol = data.frame(mane_crroi_testgenes[,"Gene.name"])
write.table(ENSgenesymbol,"temp/230414_ENSgenesymbol.txt",col.names= TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#refseq transcripts used. - 191 MANE transcripts
#select the 
colnames(MANE)
manegenes = MANE %>% dplyr::select(c("Gene.stable.ID.version","Gene.name","Transcript.stable.ID","Transcript.stable.ID.version")) %>% unique()
testgenes = subset(manegenes, (Gene.name %in% c(crdis_genes)))
testgenes$Gene.name %>% as.factor() %>% n_distinct()
write.table(testgenes,"temp/crdisgenes.txt",col.names= TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


#-------------------------------------------------------------------------------------------
#creating a data.frame for use in ANNOT process to enable annotating variants (by location) with their corresponding 'reg region'
#starting with a list of the regions of interest = mane_crroi
colnames(mane_crroi)
# add a region identifier column
mane_crroi$regionID =  paste("regregion", mane_crroi$Gene.name, sep = "_")
# add a region column in bed format (chr:start-end)
mane_crroi$loc38start = paste(mane_crroi$Chromosome.scaffold.name, mane_crroi$start, sep = ":")
mane_crroi$loc38region = paste(mane_crroi$loc38start, mane_crroi$end, sep = "-")
# annotate the file with a column to indicate if in the crdis gene list
mane_crroi$crdis280323_gene <- ifelse(mane_crroi$Gene.name %in% crdis_genes, "yes", "no")

#check
mane_crroi$Gene.name %>% as.factor() %>% n_distinct()
mane_crroi$crdis280323_gene %>% as.factor() %>% summary()
write.table(mane_crroi,"Data/230414_manegenes_annot.txt",col.names= TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

