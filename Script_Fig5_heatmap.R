####INFORMATION ####
##19.03.2024 
##Author Maddison Mackenzie

# load packages
library(dplyr)
library(readr)
library(pheatmap)

##Adding binary to df testset

testset = read.table("Data/cisreg_testset_ANNOTcategories.txt", stringsAsFactors = TRUE , sep = "\t", header = TRUE)

library(dplyr)

testset <- testset %>%
  mutate(CADD = case_when(
    is.na(CADD_PHRED) ~ 3,
    CADD_PHRED >= 10 ~ 2,
    CADD_PHRED >= 8 ~ 1,
    TRUE ~ 0
  ),
  LINSIGHT = case_when(
    is.na(linsight_Lscore_hg19) ~ 3,
    linsight_Lscore_hg19 > 0.24 ~ 2,
    linsight_Lscore_hg19 >= 0.16 ~ 1,
    TRUE ~ 0
  ),
  `FATHMM_MKL` = case_when(
    is.na(fthmmkl_Non.Coding.Score_grc37) ~ 3,
    fthmmkl_Non.Coding.Score_grc37 > 0.59 ~ 2,
    fthmmkl_Non.Coding.Score_grc37 >= 0.39 ~ 1,
    TRUE ~ 0
  ),
  `FATHMM_XF` = case_when(
    is.na(fthmxf_Non.Coding.Score_grc38) ~ 3,
    fthmxf_Non.Coding.Score_grc38 > 0.14 ~ 2,
    fthmxf_Non.Coding.Score_grc38 >= 0.12 ~ 1,
    TRUE ~ 0
  ),
  Eigen = case_when(
    is.na(eigen_Eigen.raw) ~ 3,
    eigen_Eigen.raw > 0.594 ~ 2,
    eigen_Eigen.raw >= 0.394 ~ 1,
    TRUE ~ 0
  ),
  REMM = case_when(
    is.na(REMM_score) ~ 3,
    REMM_score > 0.86 ~ 2,
    REMM_score >= 0.8 ~ 1,
    TRUE ~ 0
  ))


####making the two DF for heatmaps 
# Create disheatmap dataframe for rows with "Dis" in INFO column
disheatmap <- subset(testset, INFO == "Dis")

# Create contheatmap dataframe for rows with "Cont" in INFO column
contheatmap <- subset(testset, INFO == "Cont")

##disease heatmap
# Define your color palette
colors <- c("#0066CC", "#FFD700", "#FF3333", "#B0B0B0")

# Select the columns of interest from your dataset (disheatmap)
heatmap_data <- disheatmap %>% 
  select(variant38, CADD, LINSIGHT, `FATHMM_MKL`, FATHMM_XF, Eigen, REMM)

# Order the data based on CADD and REMM columns
heatmap_data <- heatmap_data[order(heatmap_data$CADD, heatmap_data$REMM), ]

# Extract the relevant columns for the heatmap
heatmap_matrix <- as.matrix(heatmap_data[, c("CADD", "REMM", "FATHMM_MKL", "FATHMM_XF", "Eigen", "LINSIGHT")])

# Remove row names for Y-axis
rownames(heatmap_matrix) <- rep("", nrow(heatmap_matrix))

# Create custom labels for x-axis
colnames(heatmap_matrix) <- c("CADD", "REMM", "FATHMM-MKL", "FATHMM-XF", "Eigen", "LINSIGHT")

# Create the heatmap using pheatmaps
heatmap_dis <- pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE, 
  cluster_cols = FALSE,  
  color = colors,       
  main = "Disease",
  fontsize_row = 8,      
  fontsize_col = 10,     
  angle_col = 45,       
  scale = "none",        
  border_color = NA,     
  display_numbers = FALSE,  
  legend = TRUE,         
  YAxis = FALSE          
)

#final figure

png(file="hm_dis.png")
heatmap_dis
dev.off()


#####Control heatmap 

# Define your color palette
colors <- c("#0066CC", "#FFD700", "#FF3333", "#B0B0B0")

# Select the columns of interest from your dataset
heatmap_data <- contheatmap%>% 
  select(variant38, CADD, LINSIGHT, `FATHMM_MKL`, FATHMM_XF, Eigen, REMM)

# Order the data based on CADD and REMM columns
heatmap_data <- heatmap_data[order(heatmap_data$CADD, heatmap_data$REMM), ]

# Extract the relevant columns for the heatmap
heatmap_matrix <- as.matrix(heatmap_data[, c("CADD", "REMM", "FATHMM_MKL", "FATHMM_XF", "Eigen", "LINSIGHT")])

# Remove row names for Y-axis
rownames(heatmap_matrix) <- rep("", nrow(heatmap_matrix))

# Create custom labels for x-axis
colnames(heatmap_matrix) <- c("CADD", "REMM", "FATHMM-MKL", "FATHMM-XF", "Eigen", "LINSIGHT")

# Create the heatmap using pheatmaps
heatmap_cont <- pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  color = colors,        
  main = "Control",
  fontsize_row = 8,      
  fontsize_col = 10,    
  angle_col = 45,       
  scale = "none",       
  border_color = NA,     
  display_numbers = FALSE,  
  legend = TRUE,         
  YAxis = FALSE          
)

#final figure

png(file="hm_cont.png")
heatmap_cont
dev.off()
