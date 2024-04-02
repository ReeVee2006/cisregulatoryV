#working to redo the control graph without the
##cis reg variant - control variant analysis
#fig3
#Version 1 - 17.12.23
##GERP score(Conservation) and phyloP 100V (c_phyloP100) both annotated via VEP cmd line - score headers maintained as per VEP plugin

library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(cowplot)

#setwd 

#load variants - 316425 variants - from cont partA script
afcons = read.table("testset_vepafcons.txt",sep = "\t", header = TRUE)

#load disease-associated variant set
cisregDraw = read.table("D:/cisreg_manuscript/SUPP/Data/Supplemental_DisVar.txt", header= TRUE, sep = "\t")

#prepare the data for analysis and graphing
colnames(afcons)
afcons %>% n_distinct()
summary(afcons)

#select only the rows with filter category PASS
nrow(afcons[afcons$FILTER == "PASS", ] )
afcons = afcons[afcons$FILTER =="PASS", ]
#count number of variants with maxAF not calculated ie "."
nrow(afcons[afcons$maxAF == ".", ] )
#exclude variants with maxAF not calculated ie "."
afcons = afcons[afcons$maxAF != ".", ]
#remove rows with NA in maxAF column, rows with phylop or cons NA will remain, and excluded in plots
afcons = afcons[!is.na(afcons$maxAF),]
#convert to numeric class
afcons$maxAF = as.numeric(afcons$maxAF)
afcons$Conservation = as.numeric(afcons$Conservation)
afcons$c_phyloP100 = as.numeric(afcons$c_phyloP100)
summary(afcons)

###filter so only variants in ROIs matched to final cis-reg disease variant set remain--------------------
#select the MANE gene regions of interest associated with disease variants
#load mane gene rois
#upload input file with ROI regions and associated gene name 
mane_crroi = read.table("D:/cisreg_manuscript/SUPP/Data/maneroi_annot.txt", header= TRUE, sep = "\t")

#create a column for the ROI region location
colnames(mane_crroi)
mane_crroi$loc38regionstart = paste(mane_crroi$Chromosome.scaffold.name, mane_crroi$start, sep = ":")
mane_crroi$loc38region = paste(mane_crroi$loc38regionstart, mane_crroi$end, sep = "-")

#list disease genes
#create list of genes represented in disease variant set
crdis_genes = unique(cisregDraw$ROI_gene_name)

#selecting the cisreg disease variant set ROIs
mane_crroi$crdis_gene = ifelse(mane_crroi$Gene.name %in% cisregDraw$ROI_gene_name, "crdis_genes", "NA")
mane_crroi$crdis_gene = as.factor(mane_crroi$crdis_gene) 
summary(mane_crroi$crdis_gene)

#select disease mane gene rois
colnames(mane_crroi)
matchedroi = mane_crroi[mane_crroi$crdis_gene =="crdis_genes", ]

#selecting the disease rois - disease set variant mane gene cisreg regions 
matchedroi$Gene.stable.ID.version %>% as.factor() %>% n_distinct()
manegenes = matchedroi$Gene.name
setdiff(crdis_genes, manegenes)
setdiff(matchedroi$Gene.name,crdis_genes)

#select control variants in mane gene rois
#make granges from matchedroi
#converting dataframe into GRanges object - from vcf like file
colnames(matchedroi)
summary(as.factor(matchedroi$Chromosome.scaffold.name))
matchedroi$CHROM = paste("chr",matchedroi$Chromosome.scaffold.name,sep = "")
summary(as.factor(matchedroi$CHROM))
gr_matchedroi =  makeGRangesFromDataFrame(matchedroi,
                                          keep.extra.columns=TRUE,
                                          ignore.strand=FALSE,
                                          strand.field = "Strand",
                                          seqinfo=NULL,
                                          seqnames.field=c("CHROM"),
                                          start.field=c("start"),
                                          end.field=c("end"),
                                          starts.in.df.are.0based=TRUE)

#make Granges from afconsvariants
colnames(afcons)
summary(as.factor(afcons$CHROM))
#calculating end position - support determination of 1based including all variant types
afcons$ALTcount = nchar(afcons$ALT)
afcons$REFcount = nchar(afcons$REF)
afcons$allelesize = afcons$ALTcount - afcons$REFcount
afcons$width = afcons$allelesize +1
afcons$end.field = afcons$POS + afcons$REFcount

#converting dataframe into GRanges object - from vcf like file
gr_afcons =  makeGRangesFromDataFrame(afcons,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=TRUE,
                                           seqinfo=NULL,
                                           seqnames.field=c("CHROM"),
                                           start.field=c("POS"),
                                           end.field=c("end.field"),
                                           starts.in.df.are.0based=FALSE)

#overlap dis ROIs and test set to see if test set are located within the BED file
gr_disroi_contV = subsetByOverlaps(gr_afcons, gr_matchedroi)
disroi_contV = as(gr_disroi_contV, "data.frame")
head(disroi_contV)

#save matched control variants
write.table(disroi_contV, "temp/disROImatchedcontV.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

####ANALYSIS OF MATCHED CONTROLS------------------------------------
theme_set(theme_gray(base_size = 16))

#analysis of correlation across entire range of maxAF

## repeat control analysis but with disease Roi matched control variants 306000
cor.test(disroi_contV$maxAF,disroi_contV$c_phyloP100, method = 'spearman')
cor.test(disroi_contV$maxAF,disroi_contV$Conservation, method = 'spearman')

#max AF v phyloP jitter dot plot (A)
pdot = ggplot(disroi_contV, aes(x = maxAF, y = c_phyloP100))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_log10()+labs(x = "maxAF log10", y = "phyloP 100V")

# max AF v gerp jitter dot plot (B)
#gdot
gdot = ggplot(disroi_contV, aes(x = maxAF, y = Conservation))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_log10()+labs(x = "maxAF log10", y = "GERP") 

#analysis of the binned summary stats
print("summary stats maxAF")
group_by(disroi_contV, bin) %>%
  summarise(
    count = n(),
    mean = mean(maxAF, na.rm = TRUE),
    sd = sd(maxAF, na.rm = TRUE),
    median = median(maxAF, na.rm = TRUE),
    IQR = IQR(maxAF, na.rm = TRUE),
    min = min(maxAF),
    max = max(maxAF)
  )

print("summary stats phylop")
phylopsumm = group_by(disroi_contV, bin) %>%
  summarise(
    count = n(),
    mean = mean(c_phyloP100, na.rm = TRUE),
    sd = sd(c_phyloP100, na.rm = TRUE),
    se = sd / sqrt(count),na.rm = TRUE,
    median = median(c_phyloP100, na.rm = TRUE),
    IQR = IQR(c_phyloP100, na.rm = TRUE)
  )
phylopsumm$bin_no = c(1:7)
print(phylopsumm)

#source data for supplementary table 4
write.table(phylopsumm, "temp/phyloP_summarystats.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#psum
psum = ggplot(data = phylopsumm, aes(x=bin_no, y = mean, ymin = mean - (IQR/2), ymax = mean + (IQR/2)))+
  geom_pointrange()+
  scale_x_continuous(breaks = c(1:7))+
  labs(x = "maxAF bin", y = "phyloP 100V")

print("summary stats Conservation")

gerpsumm = group_by(disroi_contV, bin) %>%
  summarise(
    count = n(),
    mean = mean(Conservation, na.rm = TRUE),
    sd = sd(Conservation, na.rm = TRUE),
    se = sd / sqrt(count),na.rm = TRUE,
    median = median(Conservation, na.rm = TRUE),
    IQR = IQR(Conservation, na.rm = TRUE)
  )
gerpsumm$bin_no = c(1:7)
print(gerpsumm)

#source data for supplementary table 5
write.table(gerpsumm, "temp/gerp_summarystats.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#gsum
gsum = ggplot(data = gerpsumm, aes(x=bin_no, y = mean, ymin = mean - (IQR/2), ymax = mean + (IQR/2)))+
  geom_pointrange()+
  scale_x_continuous(breaks = c(1:7))+
  labs(x = "maxAF bin", y = "GERP")

ungroup(disroi_contV)

kruskal.test(c_phyloP100 ~ bin, data = disroi_contV)
kruskal.test(Conservation ~ bin, data = disroi_contV)

######-------------------------------
#Fig3 correlation of maxAF with conservation scores

plot_grid(pdot, gdot, psum, gsum,
                    labels = c("A","B","C","D"),
                    align = "v", rel_widths = c(1, 1, 1,1),
                    rel_heights = c(1,1,1,1),
                    ncol = 2, nrow = 2)

ggsave(path = "Data/", filename = "Fig3_maxAFcons_3IQR.png", width = 12, height = 10, device='png', dpi=1200)

#--------------------------------------------------------
#select control set
#select variants with appropriate maxAF
colnames(afcons)
ts_select= afcons %>% filter(maxAF > 0.00002 & maxAF <= 0.0001)
ts_selectcont2= ts_select %>% filter(row_number() %% 10 == 1)

#filter to ensure that are in disease ROIs
#converting dataframe into GRanges object - from vcf like file
gr_ts =  makeGRangesFromDataFrame(ts_selectcont2,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=TRUE,
                                  seqinfo=NULL,
                                  seqnames.field=c("CHROM"),
                                  start.field=c("POS"),
                                  end.field=c("end.field"),
                                  starts.in.df.are.0based=FALSE)

#overlap dis ROIs and test set to see if test set are located within the BED file
gr_disroi_ts = subsetByOverlaps(gr_ts, gr_matchedroi)
disroi_contV_ts = as(gr_disroi_ts, "data.frame")
head(disroi_contV_ts)

#SUPP table 6 - variants selects from maxAF bin 0.00002 - 0.0001
colnames(disroi_contV_ts)
disroi_contV_ts$Contset_ID = paste("CRcont",disroi_contV_ts$contsetid,sep ="")
contV_ts = disroi_contV_ts %>% select("seqnames","start","ID","REF","ALT","QUAL","FILTER","Conservation","c_phyloP100","maxAF","bin", "Contset_ID")
names(contV_ts)[names(contV_ts) == "start"] = "POS"
names(contV_ts)[names(contV_ts) == "seqnames"] = "CHROM"
