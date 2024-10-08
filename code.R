setwd("your directory")

library(tidyverse)
library(mixOmics)
library(ggpubr)
library(vegan)
library(caret)
library(limma)
library(patchwork)
library(rstatix)
library(Spectra)
library(MsBackendMgf)
library(UpSetR)
library(lmerTest)
library(pheatmap)

# Read data
feature_table <- read_csv("gnps_quant.csv")

metadata <- read.delim("metadata.txt") %>% 
  dplyr::filter(str_starts(pattern = "A", anonymized_name)) %>% 
  distinct(anonymized_name, .keep_all = TRUE) %>%
  dplyr::select(2:6,8:9,22,24,26,33,35,39,41,46,48,53,54)
  
annotations <- read.delim("fbmn.tsv") %>%
  dplyr::filter(!str_detect(pattern = "REFRAME", LibraryName)) # remove drug library
annotations$X.Scan. <- as.character(annotations$X.Scan.)

rev_cosine <- read_delim("candidate_BA_rev_cos.tsv") %>% 
  dplyr::select(1,2,35)

info_feature <- feature_table %>% dplyr::select(1:3,7)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)

info_feature_complete <- info_feature %>% 
  left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  dplyr::select(1:5,18,24) %>% 
  left_join(rev_cosine, by = c("SpectrumID" = "query_id"))


# Feature table
data <- feature_table %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)

data$SampleID <- gsub(".mzML Peak area", "", data$SampleID)

# Metadata
metadata_metabolomics <- data_frame(SampleID = data$SampleID) %>% 
  left_join(metadata, by = c("SampleID" = "metabolomics_sample_id"))

# Investigate total peak area
metadata_TIC <- data.frame(TIC = rowSums(data %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

metadata_TIC %>% 
  dplyr::filter(!str_detect(SampleID, "(?i)qc|blank|blk")) %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order cause it was not provided
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 3e9) +
  stat_cor() 
# there two sample with low TIC. Remove them (A19_103_20_Animal_Died, A19_103_15_Animal_Died)

# Check sample type
sample_tic <- metadata_TIC %>% dplyr::filter(!str_detect(pattern = "(?i)qc|blank|blk|Died", SampleID)) %>% summarise(mean(TIC))
blank_tic <- metadata_TIC %>% dplyr::filter(str_detect(pattern = "Blank", SampleID)) %>% summarise(mean(TIC))

ratio_tic_sb <- sample_tic/blank_tic

# Check internal standard
fbmn_IS <- annotations %>% 
  dplyr::filter(str_detect(Compound_Name, regex("Amitriptyline|sulfa", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract Amitriptyline (IS) from the table
table_IS <- data %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("ID") %>% dplyr::filter(ID %in% fbmn_IS$X.Scan.) %>% 
  column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics) %>%
  dplyr::filter(!(str_detect(SampleID, "(?i)qc|blank|blk|died"))) %>%
  dplyr::select(SampleID, `40131`)

colnames(table_IS)[2] <- "Amitriptyline"

table_IS %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order cause it was not provided
  ggscatter(x = "Order", y = "Amitriptyline", add = c("reg.line")) + ylim(0, 5e5) +
  stat_cor()

table_IS %>%
  dplyr::mutate(Order = seq_len(n())) %>% # fake order cause it is not provided
  ggbarplot(x = "Order", y = "Amitriptyline", xlab = "Order", 
            ylab = "Peak Area Amitriptyline", title = "Internal Standard Acquisition") +
  geom_hline(yintercept = mean(table_IS$Amitriptyline, na.rm = TRUE), linetype = "dashed", color = "blue")

cv_is <- sd(table_IS$Amitriptyline)/mean(table_IS$Amitriptyline)

# Remove two samples with low TIC
data_final <- data %>% 
  dplyr::filter(!(SampleID %in% c("A19_103_20_Animal_Died", "A19_103_15_Animal_Died")))

# Check features in blanks and samples
data_blank <- data_final %>% dplyr::filter(str_detect(pattern = "lank", SampleID))
data_sample <- data_final %>% dplyr::filter(!(str_detect(pattern = "lank", SampleID)))

# Blanks
blanks_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                  Mean_blank = data_blank %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_blank =  data_blank %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_blank, SD_blank, CV_blank) %>% 
  dplyr::filter(Mean_blank > 0) %>% arrange(desc(Mean_blank))

# Samples
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_sample =  data_sample %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_sample, SD_sample, CV_sample) %>% 
  dplyr::filter(Mean_sample > 0) %>% arrange(desc(Mean_sample))


# Features to be removed - Sample/Blank < 5
feature_to_remove <- blanks_feature_info %>% left_join(sample_feature_info) %>%
  dplyr::filter(Mean_blank > 0) %>% 
  dplyr::mutate(Sample_Blank = Mean_sample/Mean_blank) %>% 
  dplyr::filter(Sample_Blank < 5 | is.na(Sample_Blank))

# Data with blank removal
data_clean <- data_final %>% dplyr::select(-c(feature_to_remove$Feature)) %>%
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# Remove feature before 0.2 minutes and after 7 minutes
feature_to_remove_rt <- info_feature_complete %>% dplyr::filter(RT < 0.2 | RT > 7) %>%
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))

# Final cleaned table
data_clean2 <- data_clean %>% dplyr::select(-c(feature_to_remove_rt$Feature))

# PCA on raw data
PCA_raw <- mixOmics::pca(data_clean2 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

i <- "diet"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Study 1 - Fecal Metabolome", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Keep only samples
data_sample <- data_clean2 %>% 
  dplyr::filter(!(str_detect(pattern = "(?i)qc|blank|blk", SampleID)))

# RCLR transformation
data_sample_clr <- decostand(data_sample %>% column_to_rownames("SampleID"), method = "rclr")

# Fix metadata - At the first timepoint all animals were receiving RC
metadata_fix <- metadata_metabolomics %>% 
  dplyr::mutate(diet = case_when(host_age == 10 ~ "Regular chow",
                                 TRUE ~ diet))

# PCA on samples
PCA_whole <- mixOmics::pca(data_sample_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_fix)

PCA_study_plots <- list()

for (i in c("exposure_type", "diet", "cage_location", 
            "cage_number", "host_age", "host_subject_id")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Study 1 - Fecal Metabolome", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_study_plots[[i]] <- PCA_plot
  
}

PCA_study_plots_final <- wrap_plots(PCA_study_plots, nrow = 3)

# PERMANOVA
dist_metabolites <- vegdist(data_sample_clr, method = "euclidean")
disper_diet <- betadisper(dist_metabolites, PCA_whole_scores$diet)
anova(disper_diet)
permanova <- adonis2(dist_metabolites ~ diet * exposure_type * host_age, 
                     PCA_whole_scores, na.action = na.omit)


# Save full table for joint-RPCA integration
table_metabolites <- data_sample %>% 
  column_to_rownames("SampleID") %>% 
  t() %>% as.data.frame() %>% rownames_to_column("Feature")

#write_csv(x = table_metabolites, file = "table_metabolomics_jointRPCA.csv")


#############
# Figure 3A #
#############
data_sample_no_base <- data_sample %>% 
  dplyr::filter(!(str_detect(pattern = "_01", SampleID))) # remove first timepoint

# RCLR transformation
data_sample_no_base_clr <- decostand(data_sample_no_base %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_sample <- mixOmics::pca(data_sample_no_base_clr %>%
                             select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_sample_scores <- data.frame(PCA_sample$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_fix)

# Plot
PCA_sample_plot <- PCA_sample_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "diet_exp", alpha = 0.6,
            title = "PCA - Fecal Metabolome", palette = c("#404fa2","#ed2a26","#050708"),
            xlab = paste("PC1 (", round(PCA_sample$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sample$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sample_scores %>% group_by(diet_exp) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = diet_exp), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PERMANOVA
dist_metabolites <- vegdist(data_sample_no_base_clr, method = "euclidean")
disper_diet <- betadisper(dist_metabolites, PCA_sample_scores$diet_exp)
anova(disper_diet)
permanova <- adonis2(dist_metabolites ~ diet_exp * host_age + 
                       host_subject_id, PCA_sample_scores, na.action = na.omit)

#ggsave(plot = PCA_sample_plot, filename = "PCA_all_samples.svg", device = "svg", dpi = "retina", height = 3, width = 4)


####################
# HFHC Diet Effect #
####################
# Extract sample with only air exposure
sample_air <- metadata_fix %>% 
  dplyr::filter(exposure_type == "AIR") %>%
  dplyr::filter(host_age != 10)

data_air <- data_sample %>% dplyr::filter(SampleID %in% sample_air$SampleID)

# RCLR transformation
data_air_clr <- decostand(data_air %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_air <- mixOmics::pca(data_air_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                         ncomp = 2, center = TRUE, scale = TRUE)
PCA_air_scores <- data.frame(PCA_air$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_fix)

PCA_air_plots <- list()

for (i in c("exposure_type", "diet", "cage_location", 
            "cage_number", "host_age", "host_subject_id")) {
  
  PCA_plot <- PCA_air_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Study 1 - Air exposure", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_air$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_air$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_air_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_air_plots[[i]] <- PCA_plot
  
}

PCA_air_plots_final <- wrap_plots(PCA_air_plots, nrow = 2)

# PLSDA - High Fat Diet
PLSDA_diet <- mixOmics::plsda(data_air_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                              PCA_air_scores$diet, ncomp = 2, scale = TRUE)
PLSDA_diet_scores <- data.frame(PLSDA_diet$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_fix)

PLSDA_diet_plot <- PLSDA_diet_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "diet", alpha = 0.6, title = "PLSDA - Study 1 Fecal Metabolome - Diet",
            xlab = paste("Component 1 (", round(PLSDA_diet$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_diet$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_diet_scores %>% group_by(diet) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = diet), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_diet <- plotLoadings(PLSDA_diet, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_diet <- perf(PLSDA_diet, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_diet, legend = FALSE)

VIPs_diet <- as.data.frame(mixOmics::vip(PLSDA_diet))
VIPs_diet_filter <- dplyr::filter(VIPs_diet, VIPs_diet$comp1 > 1)
VIPs_diet_filter$ID <- rownames(VIPs_diet_filter)
VIPs_diet_select <- VIPs_diet_filter %>% dplyr::select(ID, comp1)
VIPs_diet_Load <- VIPs_diet_select %>% 
  left_join(Loadings_diet, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_diet_Load, file = "SupplementaryTable1.csv")


#############
# Figure 3B #
#############
diet_rc <- VIPs_diet_Load %>% dplyr::filter(GroupContrib == "Regular chow")
diet_hfhc <- VIPs_diet_Load %>% dplyr::filter(GroupContrib == "High Fat Diet")

diet_hfhc_effect <- data_sample %>%
  dplyr::select("SampleID", diet_rc$ID, diet_hfhc$ID) %>%
  dplyr::mutate(HFHC = rowSums(select(., diet_hfhc$ID))) %>%
  dplyr::mutate(RC = rowSums(select(., diet_rc$ID))) %>%
  dplyr::mutate(Ratio = log(RC/HFHC)) %>%
  left_join(metadata_fix) %>%
  dplyr::mutate(diet_exp = factor(diet_exp, levels = c("RC_AIR", "HFD_AIR", "HFD_IHC")))

# Calculate the mean per group and age
df_mean <- diet_hfhc_effect %>%
  group_by(diet_exp, host_age) %>%
  summarise(Mean_Measurement = mean(Ratio))

# Plot
ratio_hfhc_plot <- ggplot(diet_hfhc_effect, aes(x = host_age, y = Ratio,)) +
  geom_line(aes(alpha = 0.3, colour = diet_exp, group = host_subject_id), linewidth = 0.3) +
  geom_line(data = df_mean, aes(x = host_age, y = Mean_Measurement, color = diet_exp), 
            linewidth = 2.5) +
  scale_color_manual(values = c("black", "red", "blue")) +
  theme_bw() +
  labs(x = "Host Age (weeks)", y = "Log(Top RC/Top HFHC") +
  theme(legend.position = "none")

#ggsave(plot = ratio_hfhc_plot, filename = "Figure3B.svg", device = "svg", dpi = "retina", height = 3, width = 5)

# Linear mixed effect model to test separation
model_ratio <- lmerTest::lmer(Ratio ~ host_age + diet_exp + (1 | mouse_number), data = diet_hfhc_effect)

# Summarize the model to see the results
summary(model_ratio)


# Separate samples based on diet (removing first timepoint)
sample_rc <- metadata_fix %>% dplyr::filter(diet == "Regular chow") %>%
  dplyr::filter(host_age != 10)
sample_hfhc <- metadata_fix %>% dplyr::filter(diet == "High Fat Diet") %>%
  dplyr::filter(host_age != 10)


##############
# AIR in HFD #
##############
data_hfhc <- data_sample %>% dplyr::filter(SampleID %in% sample_hfhc$SampleID)

# RCLR transformation
data_hfhc_clr <- decostand(data_hfhc %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_hfd <- mixOmics::pca(data_hfhc_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                         ncomp = 2, center = TRUE, scale = TRUE)
PCA_hfd_scores <- data.frame(PCA_hfd$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

PCA_hfd_plots <- list()

for (i in c("exposure_type", 
            "cage_number", "host_age", "host_subject_id")) {
  
  PCA_plot <- PCA_hfd_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Study 1", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_hfd$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_hfd$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_hfd_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_hfd_plots[[i]] <- PCA_plot
  
}

PCA_hfd_plots_final <- wrap_plots(PCA_hfd_plots, nrow = 2)

# PERMANOVA
dist_metabolites <- vegdist(data_hfhc_clr, method = "euclidean")
disper_air <- betadisper(dist_metabolites, PCA_hfd_scores$exposure_type)
anova(disper_air)
permanova <- adonis2(dist_metabolites ~ exposure_type * host_age, PCA_hfd_scores, na.action = na.omit)

# PLSDA - HFD Air
PLSDA_air <- mixOmics::plsda(data_hfhc_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                             PCA_hfd_scores$exposure_type, ncomp = 3, scale = TRUE)
PLSDA_air_scores <- data.frame(PLSDA_air$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_fix)

PLSDA_air_plot <- PLSDA_air_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "exposure_type", alpha = 0.6, title = "PLSDA - Study 1 Fecal Metabolome HFD - Air Exposure",
            xlab = paste("Component 1 (", round(PLSDA_air$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_air$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_air_scores %>% group_by(exposure_type) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = exposure_type), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_air <- plotLoadings(PLSDA_air, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_air <- perf(PLSDA_air, validation = "Mfold", folds = 5, nrepeat = 99, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_air, legend = FALSE)

VIPs_air <- as.data.frame(mixOmics::vip(PLSDA_air))
VIPs_air_filter <- dplyr::filter(VIPs_air, VIPs_air$comp1 > 1)
VIPs_air_filter$ID <- rownames(VIPs_air_filter)
VIPs_air_select <- VIPs_air_filter %>% dplyr::select(ID, comp1)
VIPs_air_Load <- VIPs_air_select %>% 
  left_join(Loadings_air, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#write_csv(x = VIPs_air_Load, file = "SupplementaryTable2.csv")


# Extract features associated with Air or IHC
exposure_air <- VIPs_air_Load %>% dplyr::filter(GroupContrib == "AIR") 
exposure_ihc <- VIPs_air_Load %>% dplyr::filter(GroupContrib == "IHC") 


##################################################
# Feature affected by HFD and exacerbated by IHC #
##################################################
feature_interest <- VIPs_diet_Load %>% dplyr::select(ID, comp1, GroupContrib) %>%
  dplyr::rename(vip_hfhc = "comp1", group_diet = "GroupContrib") %>% 
  full_join(VIPs_air_Load) %>% 
  dplyr::rename(vip_air = "comp1", group_air = "GroupContrib") %>%
  dplyr::mutate(combined_vip = rowSums(cbind(vip_hfhc, vip_air), na.rm = TRUE)) %>%
  dplyr::mutate(Group = paste(group_diet, group_air, sep = "_")) %>%
  dplyr::mutate(Outcome = case_when(group_diet == "Regular chow" & is.na(group_air) ~ "Decreased HFHC",
                                    group_diet == "High Fat Diet" & is.na(group_air) ~ "Increased HFHC",
                                    group_diet == "High Fat Diet" & group_air == "AIR"  ~ "Increased HFHC and Decresead IHC", 
                                    group_diet == "High Fat Diet" & group_air == "IHC"  ~ "Increased HFHC and Increased IHC",
                                    is.na(group_diet) & group_air == "AIR"  ~ "Decreased IHC",
                                    is.na(group_diet) & group_air == "IHC"  ~ "Increased IHC",
                                    group_diet == "Regular chow" & group_air == "AIR"  ~ "Decreased HFHC and Decresead IHC",
                                    TRUE ~ "remove")) %>% dplyr::filter(Outcome != "remove")

summary(factor(feature_interest$Outcome))

feature_rc_air <- feature_interest %>% 
  dplyr::filter(Outcome == "Decreased HFHC and Decresead IHC")

feature_hfd_ihc <- feature_interest %>% 
  dplyr::filter(Outcome == "Increased HFHC and Increased IHC")

# Extract exacerbated features
final_features <- feature_interest %>% 
  dplyr::filter(Outcome %in% c("Increased HFHC and Increased IHC", 
                               "Decreased HFHC and Decresead IHC")) %>%
  arrange(desc(combined_vip))

#write_csv(x = final_features, file = "SupplementaryTable3.csv")


#############
# Figure 3C #
#############
data_hfhc_ihc <- data_sample %>%
  dplyr::select("SampleID", feature_rc_air$ID, feature_hfd_ihc$ID) %>%
  dplyr::mutate(HFHC_IHC = rowSums(select(., feature_hfd_ihc$ID))) %>%
  dplyr::mutate(RC_Air = rowSums(select(., feature_rc_air$ID))) %>%
  dplyr::mutate(Ratio = log(RC_Air/HFHC_IHC)) %>%
  left_join(metadata_fix) %>%
  dplyr::mutate(diet_exp = factor(diet_exp, levels = c("RC_AIR", "HFD_AIR", "HFD_IHC")))

data_hfhc_ihc_plot <- data_hfhc_ihc %>% dplyr::filter(host_age > 15) %>%
  ggboxplot(x = "diet_exp", y = "Ratio", legend = "none",
            title = "Driving HFHC and IHC metabolic features", 
            ylab = "Log(Top RC/Top HFHC IHC)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) + 
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif", ) + 
  theme(plot.title = element_text(size = 10), axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6), axis.text.x = element_blank(), axis.ticks.x = element_blank())

anova <- aov(Ratio ~ diet_exp, data_hfhc_ihc)
tukey_hsd(anova)

#ggsave(plot = data_hfhc_ihc_plot, filename = "Figure3C.svg", device = "svg", dpi = "retina", width = 2.5, height = 2.5)

#write_csv(data_hfhc_ihc %>% dplyr::select("SampleID", "Ratio"), file = "metabolome_log_ratio.csv")


#############
# Figure 3E #
#############
mean_peak_long <- data_sample %>%
  left_join(metadata_fix) %>%
  dplyr::filter(host_age != 10) %>%
  dplyr::select(feature_interest$ID, "host_subject_id") %>%
  group_by(host_subject_id) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  left_join(metadata_fix %>% dplyr::select(host_subject_id, diet_exp) %>% distinct(host_subject_id, .keep_all = TRUE)) %>%
  dplyr::mutate(diet_exp = factor(diet_exp, levels = c("RC_AIR", "HFD_AIR", "HFD_IHC"))) %>%
  pivot_longer(cols = 2:3743, names_to = "Feature", values_to = "PeakArea") %>%
  dplyr::mutate(log10_PeakArea = log10(PeakArea)) %>%
  left_join(feature_interest %>% dplyr::select(ID, Compound_Name), by = c("Feature" = "ID")) %>%
  dplyr::mutate(Compound_Name = gsub("Candidate ", "", Compound_Name)) %>%
  dplyr::mutate(Compound_Name = case_when(str_detect(pattern = "bile acid", Compound_Name) & 
                                            str_detect(pattern = "spectral library match", Compound_Name) ~ gsub(".*;(.*)", "\\1", Compound_Name),
                                          TRUE ~ Compound_Name))  %>%
  dplyr::mutate(Compound_Name = case_when(str_detect(pattern = "spectral library match", Compound_Name) ~ gsub("^[^,]*,", "", Compound_Name),
                                          TRUE ~ Compound_Name)) %>%
  dplyr::mutate(Compound_Name = gsub("Spectral Match to ", "", Compound_Name)) %>%
  dplyr::mutate(Compound_Name = gsub(" from NIST14", "", Compound_Name)) %>% 
  dplyr::mutate(Compound_Name = paste(Compound_Name, Feature, sep = "_"))


# Plot linoleic acid
linoleic_plot <- mean_peak_long %>% dplyr::filter(Feature == "30197") %>%
  ggboxplot(x = "diet_exp", y = "log10_PeakArea", legend = "none",
            title = "Linoleic Acid", 
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) + ylim(2.5,7.5) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif", ) + 
  theme(plot.title = element_text(size = 10), axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6), axis.text.x = element_blank(), axis.ticks.x = element_blank())

anova <- aov(log10_PeakArea ~ diet_exp, (mean_peak_long %>% dplyr::filter(Feature == "30197")))
tukey_hsd(anova)

# Plot Diacetylspermidine
sperm_plot <- mean_peak_long %>% dplyr::filter(Feature == "2005") %>%
  ggboxplot(x = "diet_exp", y = "log10_PeakArea", legend = "none",
            title = "Diacetylspermidine", 
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) + ylim(2.5,7.5) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif", ) + 
  theme(plot.title = element_text(size = 10), axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6), axis.text.x = element_blank(), axis.ticks.x = element_blank())

anova <- aov(log10_PeakArea ~ diet_exp, (mean_peak_long %>% dplyr::filter(Feature == "2005")))
tukey_hsd(anova)

# Plot Leucine-C16:0
leu_plot <- mean_peak_long %>% dplyr::filter(Feature == "33139") %>%
  ggboxplot(x = "diet_exp", y = "log10_PeakArea", legend = "none",
            title = "Isoleucine/Leucine-C16:0", 
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) + ylim(2.5,7.5) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif", ) + 
  theme(plot.title = element_text(size = 10), axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6), axis.text.x = element_blank(), axis.ticks.x = element_blank())

anova <- aov(log10_PeakArea ~ diet_exp, (mean_peak_long %>% dplyr::filter(Feature == "33139")))
tukey_hsd(anova)

# Plot PC(20:4)
pc_plot <- mean_peak_long %>% dplyr::filter(Feature == "29177") %>%
  ggboxplot(x = "diet_exp", y = "log10_PeakArea", legend = "none",
            title = "PC(20:4)", 
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) + ylim(2.5,7.5) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif", ) + 
  theme(plot.title = element_text(size = 10), axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6), axis.text.x = element_blank(), axis.ticks.x = element_blank())

anova <- aov(log10_PeakArea ~ diet_exp, (mean_peak_long %>% dplyr::filter(Feature == "29177")))
tukey_hsd(anova)

# Plot Taurocholic acid
tau_plot <- mean_peak_long %>% dplyr::filter(Feature == "20406") %>%
  ggboxplot(x = "diet_exp", y = "log10_PeakArea", legend = "none",
            title = "Taurine-(OH)3", 
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) + ylim(2.5,7.5) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif", ) + 
  theme(plot.title = element_text(size = 10), axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6), axis.text.x = element_blank(), axis.ticks.x = element_blank())

anova <- aov(log10_PeakArea ~ diet_exp, (mean_peak_long %>% dplyr::filter(Feature == "20406")))
tukey_hsd(anova)

# Combined plots
plot_combined <- ggarrange(linoleic_plot, sperm_plot, leu_plot, pc_plot, tau_plot, nrow = 1)

#ggsave(plot = plot_combined, filename = "Figure3E.svg", device = "svg", dpi = "retina", width = 7, height = 2)


##############
# Joint-RPCA #
##############

# Outputs from RPCA (Axis 1 separates Diet, Axis 2 separates Exposure)
rpca_table_1_axis1 <- read.delim("rpca/selected_features_1_percent_axis1.tsv")
rpca_table_10_axis1 <- read.delim("rpca/selected_features_10_percent_axis1.tsv")
rpca_table_1_axis2 <- read.delim("rpca/selected_features_1_percent_axis2.tsv")
rpca_table_10_axis2 <- read.delim("rpca/selected_features_10_percent_axis2.tsv")

colnames(rpca_table_1_axis1)[1] <- "ID"
rpca_table_1_axis1$ID <- as.character(rpca_table_1_axis1$ID)
colnames(rpca_table_10_axis1)[1] <- "ID"
rpca_table_10_axis1$ID <- as.character(rpca_table_10_axis1$ID)
colnames(rpca_table_1_axis2)[1] <- "ID"
rpca_table_1_axis2$ID <- as.character(rpca_table_1_axis2$ID)
colnames(rpca_table_10_axis2)[1] <- "ID"
rpca_table_10_axis2$ID <- as.character(rpca_table_10_axis2$ID)

# Check overlap with features identified by stratified PLS-DA models
plsda_rpca_overlap_1 <- final_features %>% 
  left_join(rpca_table_1_axis1) %>% 
  left_join(rpca_table_1_axis2, by = c("ID" = "ID"))
plsda_rpca_overlap_10 <- final_features %>%
  left_join(rpca_table_10_axis1) %>% 
  left_join(rpca_table_10_axis2, by = c("ID" = "ID"))

# Check features only found by RPCA
rpca_1_plsda_miss_axis1 <- rpca_table_1_axis1 %>% anti_join(final_features) %>% 
  left_join(VIPs_diet_Load)
rpca_10_plsda_miss_axis1 <- rpca_table_10_axis1 %>% anti_join(final_features) %>% 
  left_join(VIPs_diet_Load)
rpca_1_plsda_miss_axis2 <- rpca_table_1_axis2 %>% anti_join(final_features) %>% 
  left_join(VIPs_air_Load)
rpca_10_plsda_miss_axis2 <- rpca_table_10_axis2 %>% anti_join(final_features) %>% 
  left_join(VIPs_air_Load)

# Identified features by RPCA are either only associated to diet or air exposure
# and are recovered by general plsda models

# Examples picked by RPCA
data_sample %>% dplyr::select("SampleID", "1882") %>%
  left_join(metadata_fix) %>% 
  dplyr::filter(host_age != 10) %>%
  dplyr::select("host_subject_id", "1882") %>%
  group_by(host_subject_id) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  left_join(metadata_fix %>% dplyr::select(host_subject_id, diet_exp) %>% 
              distinct(host_subject_id, .keep_all = TRUE)) %>%
  dplyr::mutate(LogPeak = log(`1882`+1)) %>%
  ggboxplot(x = "diet_exp", y = "LogPeak", legend = "none",
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif")

data_sample %>% dplyr::select("SampleID", "6114") %>%
  left_join(metadata_fix) %>% 
  dplyr::filter(host_age != 10) %>%
  dplyr::select("host_subject_id", "6114") %>%
  group_by(host_subject_id) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  left_join(metadata_fix %>% dplyr::select(host_subject_id, diet_exp) %>% 
              distinct(host_subject_id, .keep_all = TRUE)) %>%
  dplyr::mutate(LogPeak = log(`6114`+1)) %>%
  ggboxplot(x = "diet_exp", y = "LogPeak", legend = "none",
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif")

data_sample %>% dplyr::select("SampleID", "14002") %>%
  left_join(metadata_fix) %>% 
  dplyr::filter(host_age != 10) %>%
  dplyr::select("host_subject_id", "14002") %>%
  group_by(host_subject_id) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  left_join(metadata_fix %>% dplyr::select(host_subject_id, diet_exp) %>% 
              distinct(host_subject_id, .keep_all = TRUE)) %>%
  dplyr::mutate(LogPeak = log(`14002`+1)) %>%
  ggboxplot(x = "diet_exp", y = "LogPeak", legend = "none",
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif")

data_sample %>% dplyr::select("SampleID", "26036") %>%
  left_join(metadata_fix) %>% 
  dplyr::filter(host_age != 10) %>%
  dplyr::select("host_subject_id", "26036") %>%
  group_by(host_subject_id) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  left_join(metadata_fix %>% dplyr::select(host_subject_id, diet_exp) %>% 
              distinct(host_subject_id, .keep_all = TRUE)) %>%
  dplyr::mutate(LogPeak = log(`26036`+1)) %>%
  ggboxplot(x = "diet_exp", y = "LogPeak", legend = "none",
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif")

# Examples NOT picked by RPCA
data_sample %>% dplyr::select("SampleID", "29881") %>%
  left_join(metadata_fix) %>% 
  dplyr::filter(host_age != 10) %>%
  dplyr::select("host_subject_id", "29881") %>%
  group_by(host_subject_id) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  left_join(metadata_fix %>% dplyr::select(host_subject_id, diet_exp) %>% 
              distinct(host_subject_id, .keep_all = TRUE)) %>%
  dplyr::mutate(LogPeak = log(`29881`+1)) %>%
  ggboxplot(x = "diet_exp", y = "LogPeak", legend = "none",
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif")

data_sample %>% dplyr::select("SampleID", "33139") %>%
  left_join(metadata_fix) %>% 
  dplyr::filter(host_age != 10) %>%
  dplyr::select("host_subject_id", "33139") %>%
  group_by(host_subject_id) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  left_join(metadata_fix %>% dplyr::select(host_subject_id, diet_exp) %>% 
              distinct(host_subject_id, .keep_all = TRUE)) %>%
  dplyr::mutate(LogPeak = log(`33139`+1)) %>%
  ggboxplot(x = "diet_exp", y = "LogPeak", legend = "none",
            ylab = "Log10(Peak Area)", xlab = FALSE,
            add = "jitter", palette = c("#050708","#404fa2","#ed2a26"),
            add.params = list(color = "diet_exp", alpha = 0.6)) +
  stat_compare_means(comparisons = list(c("RC_AIR", "HFD_AIR"), c("RC_AIR", "HFD_IHC"),
                                        c("HFD_AIR", "HFD_IHC")), label = "p.signif")


# Output joint-RPCA
joint_rpca_interest <- read_csv("rpca/sz-subset_jointRPCA.csv")
asv_annotation <- read_csv("rpca/shortlist_taxonomy.csv") %>%
  dplyr::mutate(Taxon = gsub(".*f__", "", Taxon))

joint_rpca_interest$featureid <- as.character(joint_rpca_interest$featureid)

final_features_joint <- final_features %>% 
  dplyr::mutate(name = paste(Outcome, ID, round(combined_vip, digits = 2), sep = "_")) %>%
  dplyr::mutate(name = gsub(" ", "_", name))

joint_output <- joint_rpca_interest %>% 
  left_join(final_features_joint %>% dplyr::select(ID, name), by = c("featureid" = "ID")) %>%
  dplyr::select(-featureid) %>% dplyr::filter(!(is.na(name))) %>%
  column_to_rownames("name") %>% t() %>% as.data.frame() %>% rownames_to_column("ASV") %>% 
  left_join(asv_annotation) %>% dplyr::select(-ASV) %>% dplyr::mutate(Taxon = paste(Taxon, seq(1:38), sep = "_")) %>%
  column_to_rownames("Taxon")

# Create the heatmap
heatmpa_rpca <- pheatmap(joint_output,
         #scale = "row",          
         clustering_distance_rows = "euclidean",  # Distance method for clustering rows
         clustering_distance_cols = "euclidean",  # Distance method for clustering columns
         clustering_method = "complete",          # Clustering method (hierarchical)
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))  # Color palette

# Save as a PDF
pdf("heatmap_rpca.pdf", width = 50, height = 20)
heatmpa_rpca
dev.off()
