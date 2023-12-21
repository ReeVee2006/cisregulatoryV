#### making graphs with thresholds cont/dis for each tool 

##CADD_PHRED

library(ggplot2)
library(dplyr)

# Adjust thresholds
threshold1 <- 8
threshold2 <- 10

# Set colors for Control and Disease groups
control_color <- "blue"
disease_color <- "red"

# Create the density plot
ggplot(data = df, aes(x = CADD_PHRED, group = INFO, fill = INFO, color = INFO)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(x = "CADD_PHRED", y = "Density", title = "Density Plot") +
  scale_fill_manual(values = c("Control" = control_color, "Dis" = disease_color), name = "Group") +
  scale_color_manual(values = c("Control" = control_color, "Dis" = disease_color), name = "Group") +
  geom_vline(xintercept = c(threshold1, threshold2), linetype = "dotted", color = "black") +
  annotate("text", x = threshold1, y = -0.01, label = as.character(threshold1), vjust = 1) +
  annotate("text", x = threshold2, y = -0.01, label = as.character(threshold2), vjust = 1)


ggplot(data=df, aes(x=CADD_PHRED, group=INFO, fill=INFO)) +
  geom_density(adjust=1, alpha=.4)


# Create the density plot
ggplot(data = df, aes(x = CADD_PHRED, group = INFO, fill = INFO)) +
  geom_density(adjust = 1, alpha = 0.4) +
  theme_minimal() +
  labs(x = "CADD_PHRED", y = "Density", title = "CADD_PHRED Control and Disease density plot") +
  scale_fill_manual(values = c("Cont" = "blue", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  scale_color_manual(values = c("Cont" = "blue", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  geom_vline(xintercept = c(8, 10), linetype = "dotted", color = "black") +
  annotate("text", x = 8, y = -0.01, label = "8", vjust = 0.5) +
  annotate("text", x = 10, y = -0.01, label = "10", vjust = 0.5)

# Create the density plot
ggplot(data = df, aes(x = CADD_PHRED, group = INFO, fill = INFO)) +
  geom_density(adjust = 1.3, alpha = 0.4, na.rm = TRUE) +
  theme_minimal() +
  labs(x = "CADD_PHRED", y = "Denstiy", title = "CADD_PHRED Control and Disease Density Plot") +
  scale_fill_manual(values = c("Cont" = "blue", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  scale_color_manual(values = c("Cont" = "blue", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  geom_vline(xintercept = c(8, 10), linetype = "dashed", color = "black") +
  geom_text(data = df[!is.na(df$CADD_PHRED) & df$CADD_PHRED == 8, ], x = 8, y = -0.01, label = "8", vjust = 1) +
  geom_text(data = df[!is.na(df$CADD_PHRED) & df$CADD_PHRED == 10, ], x = 10, y = -0.01, label = "10", vjust = 1) +
  ylim(0, max(density(df$CADD_PHRED, na.rm = TRUE)$y))

#########################################################
######Density plot START CODE CADD
ggplot(data = df, aes(x = CADD_PHRED, group = INFO, fill = INFO)) +
  geom_density(adjust = 1.3, alpha = 0.4, na.rm = TRUE) +
  theme_minimal() +
  labs(x = "CADD_PHRED", y = "Density", title = "CADD_PHRED Control and Disease Density Plot") +
  scale_fill_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  scale_color_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  geom_vline(xintercept = c(8, 10), linetype = "dashed", color = "black") +
  geom_text(data = df[!is.na(df$CADD_PHRED) & df$CADD_PHRED == 8, ], x = 8, y = -0.01, label = "8", vjust = 1) +
  geom_text(data = df[!is.na(df$CADD_PHRED) & df$CADD_PHRED == 10, ], x = 10, y = -0.01, label = "10", vjust = 1) +
  ylim(0, max(density(df$CADD_PHRED, na.rm = TRUE)$y))


##Linsight


ggplot(data = df, aes(x = linsight_Lscore_hg19, group = INFO, fill = INFO)) +
  geom_density(adjust = 1.3, alpha = 0.4, na.rm = TRUE) +
  theme_minimal() +
  labs(x = "LINSIGHT Score", y = "Density", title = "LINSIGHT Scores Control and Disease Density Plot") +
  scale_fill_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  scale_color_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  geom_vline(xintercept = c(0.16, 0.24), linetype = "dashed", color = "black") +
  geom_text(data = df[!is.na(df$linsight_Lscore_hg19) & df$linsight_Lscore_hg19 == 0.16, ], x = 0.16, y = -0.01, label = "0.16", vjust = 1) +
  geom_text(data = df[!is.na(df$linsight_Lscore_hg19) & df$linsight_Lscore_hg19 == 0.24, ], x = 0.24, y = -0.01, label = "0.24", vjust = 1) +
  scale_y_continuous(limits = c(0, NA))


##fathmm-mkl


ggplot(data = df, aes(x = fthmmkl_Non.Coding.Score_hg19, group = INFO, fill = INFO)) +
  geom_density(adjust = 1.3, alpha = 0.4, na.rm = TRUE) +
  theme_minimal() +
  labs(x = "FATHMM-MKL NonCoding Score", y = "Density", title = "FATHMM-MKL NonCoding Score Control and Disease Density Plot") +
  scale_fill_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  scale_color_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  geom_vline(xintercept = c(0.39, 0.59), linetype = "dashed", color = "black") +
  geom_text(data = df[!is.na(df$fthmmkl_Non.Coding.Score_hg19) & df$fthmmkl_Non.Coding.Score_hg19 == 0.39, ], x = 0.39, y = -0.01, label = "0.39", vjust = 1) +
  geom_text(data = df[!is.na(df$fthmmkl_Non.Coding.Score_hg19) & df$fthmmkl_Non.Coding.Score_hg19 == 0.59, ], x = 0.59, y = -0.01, label = "0.59", vjust = 1) +
  ylim(0, max(density(df$fthmmkl_Non.Coding.Score_hg19, na.rm = TRUE)$y))

##FATHMM-XF 

ggplot(data = df, aes(x = fthmxf_Non.Coding.Score_hg19, group = INFO, fill = INFO)) +
  geom_density(adjust = 1.3, alpha = 0.4, na.rm = TRUE) +
  theme_minimal() +
  labs(x = "FATHMM-XF NonCoding Score", y = "Density", title = "FATHMM-XF NonCoding Score Control and Disease Density Plot") +
  scale_fill_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  scale_color_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  geom_vline(xintercept = c(0.12, 0.14), linetype = "dashed", color = "black") +
  geom_text(data = df[!is.na(df$fthmxf_Non.Coding.Score_hg19) & df$fthmxf_Non.Coding.Score_hg19 == 0.12, ], x = 0.12, y = -0.01, label = "8", vjust = 1) +
  geom_text(data = df[!is.na(df$fthmxf_Non.Coding.Score_hg19) & df$fthmxf_Non.Coding.Score_hg19 == 0.14, ], x = 0.14, y = -0.01, label = "10", vjust = 1) +
  ylim(0, max(density(df$fthmxf_Non.Coding.Score_hg19, na.rm = TRUE)$y))



##Eigen 
ggplot(data = df, aes(x = eigen_Eigen.raw, group = INFO, fill = INFO)) +
  geom_density(adjust = 1.3, alpha = 0.4, na.rm = TRUE) +
  theme_minimal() +
  labs(x = "Eigen raw Score", y = "Density", title = "Eigen raw Scores Control and Disease Density Plot") +
  scale_fill_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  scale_color_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  geom_vline(xintercept = c(0.394, 0.594), linetype = "dashed", color = "black") +
  geom_text(data = df[!is.na(df$eigen_Eigen.raw) & df$eigen_Eigen.raw == 0.394, ], x = 0.394, y = -0.01, label = "8", vjust = 1) +
  geom_text(data = df[!is.na(df$eigen_Eigen.raw) & df$eigen_Eigen.raw == 0.594, ], x = 0.594, y = -0.01, label = "10", vjust = 1) +
  ylim(0, max(density(df$eigen_Eigen.raw, na.rm = TRUE)$y))


##REMM


ggplot(data = df, aes(x = REMM_score, group = INFO, fill = INFO)) +
  geom_density(adjust = 1.3, alpha = 0.4, na.rm = TRUE) +
  theme_minimal() +
  labs(x = "REMM Score", y = "Density", title = "REMM Score Control and Disease Density Plot") +
  scale_fill_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  scale_color_manual(values = c("Cont" = "#00AABB", "Dis" = "red"), labels = c("Control", "Disease"), name = "Group") +
  geom_vline(xintercept = c(0.8, 0.86), linetype = "dashed", color = "black") +
  geom_text(data = df[!is.na(df$REMM_score) & df$REMM_score == 0.8, ], x = 0.8, y = -0.01, label = "8", vjust = 1) +
  geom_text(data = df[!is.na(df$REMM_score) & df$REMM_score == 0.86, ], x = 0.86, y = -0.01, label = "10", vjust = 1) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0))




