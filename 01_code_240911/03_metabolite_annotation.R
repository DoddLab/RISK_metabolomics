################################################################################
# Metabolite annotation --------------------------------------------------------

# the metabolites were annoated by match MS1, RT, and MS2 with in-house and public database
# Level 1: MS1, RT, and MS2 match with in-house library
# Level 2.1: MS1 and RT match with in-house library
# Level 2.2: MS1 and MS2 match with public library

# In-house library: DoddLib
# Public library: MSDial_Lib, GNPS_Bile_Acid_Lib, GNPS_Acyl_Amides_Lib, GNPS_Acyl_Esters, Fiehn_Peptide_Lib 

# The annotated metabolites were merged and reduced the replication by confidence level and manually check the annotation.

# Packages:
# DoddLabMetID: v0.1.19, https://github.com/DoddLab/DoddLabMetID 
# DoddLabDatabase: v0.2.6, https://github.com/DoddLab/DoddLabDatabase


################################################################################
# Libraries match ---------------------------------------------------------------

# C18 pos ----------------------------------------------------------------------
setwd('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/')


# library(DoddLabPackages)
library(DoddLabMetID)
library(DoddLabDatabase)

load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/01_input_data_cleaning/05_object_c18_pos_outlier_removal.RData')

# level 1: mz + RT + ms2
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/',
                                                                  lib = 'dodd',
                                                                  column = 'c18',
                                                                  polarity = 'positive',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_c18_pos, parameter_set_annotation)

file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')



# mz + RT
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/01_input_data_cleaning/05_object_c18_pos_outlier_removal.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/',
                                                                  lib = 'dodd',
                                                                  column = 'c18',
                                                                  polarity = 'positive',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = FALSE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_c18_pos, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt')

# mz + ms2
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/01_input_data_cleaning/05_object_c18_pos_outlier_removal.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/',
                                                                  lib = 'all_public',
                                                                  column = 'c18',
                                                                  polarity = 'positive',
                                                                  is_rt_score = FALSE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.7
parameter_set_annotation@para_ms2_match$direction <- 'forward'
parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
annotate_metabolite(object = object_c18_pos, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_mz_ms2_public_db')



# C18 neg ----------------------------------------------------------------------

# library(DoddLabPackages)
library(DoddLabMetID)
library(DoddLabDatabase)

setwd('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_neg/')
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_neg/01_input_data_cleaning/05_object_c18_neg_outlier_removal.RData')

# mz + RT + ms2
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_neg/',
                                                                  lib = 'dodd',
                                                                  column = 'c18',
                                                                  polarity = 'negative',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_c18_neg, parameter_set_annotation)

file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')



# mz + RT
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_neg/01_input_data_cleaning/05_object_c18_neg_outlier_removal.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_neg/',
                                                                  lib = 'dodd',
                                                                  column = 'c18',
                                                                  polarity = 'negative',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = FALSE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_c18_neg, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt')

# mz + ms2
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_neg/01_input_data_cleaning/05_object_c18_neg_outlier_removal.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_neg/',
                                                                  lib = 'all_public',
                                                                  column = 'c18',
                                                                  polarity = 'negative',
                                                                  is_rt_score = FALSE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.7
parameter_set_annotation@para_ms2_match$direction <- 'forward'
parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
annotate_metabolite(object = object_c18_neg, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_mz_ms2_public_db')


# hilic pos ----------------------------------------------------------------------

# library(DoddLabPackages)
library(DoddLabMetID)
library(DoddLabDatabase)

setwd('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_pos/')
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_pos/01_input_data_cleaning/05_object_hilic_pos_outlier_removal.RData')

# mz + RT + ms2
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_pos/',
                                                                  lib = 'dodd',
                                                                  column = 'hilic',
                                                                  polarity = 'positive',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$tolerance_rt_range <- 30
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_hilic_pos, parameter_set_annotation)

file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')



# mz + RT
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_pos/01_input_data_cleaning/05_object_hilic_pos_outlier_removal.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_pos/',
                                                                  lib = 'dodd',
                                                                  column = 'hilic',
                                                                  polarity = 'positive',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = FALSE)

parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$tolerance_rt_range <- 30
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_hilic_pos, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt')

# mz + ms2
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_pos/01_input_data_cleaning/05_object_hilic_pos_outlier_removal.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_pos/',
                                                                  lib = 'all_public',
                                                                  column = 'hilic',
                                                                  polarity = 'positive',
                                                                  is_rt_score = FALSE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.7
parameter_set_annotation@para_ms2_match$direction <- 'forward'
parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
annotate_metabolite(object = object_hilic_pos, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_mz_ms2_public_db')



# hilic neg ----------------------------------------------------------------------

# library(DoddLabPackages)
library(DoddLabMetID)
library(DoddLabDatabase)

setwd('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_neg/')
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_neg/01_input_data_cleaning/05_object_hilic_neg_outlier_removal.RData')

# mz + RT + ms2
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_neg/',
                                                                  lib = 'dodd',
                                                                  column = 'hilic',
                                                                  polarity = 'negative',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$tolerance_rt_range <- 30
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_hilic_neg, parameter_set_annotation)

file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')



# mz + RT
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_neg/01_input_data_cleaning/05_object_hilic_neg_outlier_removal.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(object = object_hilic_neg,
                                                                  path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_neg/',
                                                                  lib = 'dodd',
                                                                  column = 'hilic',
                                                                  polarity = 'negative',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = FALSE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$tolerance_rt_range <- 30
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
annotate_metabolite(object = object_hilic_neg, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt')

# mz + ms2
load('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_neg/01_input_data_cleaning/05_object_hilic_neg_outlier_removal.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(object = object_hilic_neg,
                                                                  path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_neg/',
                                                                  lib = 'all_public',
                                                                  column = 'hilic',
                                                                  polarity = 'negative',
                                                                  is_rt_score = FALSE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.7
parameter_set_annotation@para_ms2_match$direction <- 'forward'
parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
annotate_metabolite(object = object_hilic_neg, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_mz_ms2_public_db')



################################################################################
# Merge annotation results -----------------------------------------------------
library(tidyverse)
library(DoddLabMetID)

merge_one_modes(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_pos/',column = 'c18', polarity = 'positive')
merge_one_modes(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/c18_neg/',column = 'c18', polarity = 'negative')
merge_one_modes(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_pos/',column = 'hilic', polarity = 'positive')
merge_one_modes(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/hilic_neg/',column = 'hilic', polarity = 'negative')


hilic_pos <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/HILIC_pos/annotation_table_merge_hilic_positive.xlsx')
hilic_neg <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/HILIC_neg/annotation_table_merge_hilic_negative.xlsx')
c18_pos <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/C18_pos/annotation_table_merge_c18_positive.xlsx')
c18_neg <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/C18_neg/annotation_table_merge_c18_negative.xlsx')

temp_merge <- hilic_pos %>% 
  dplyr::bind_rows(hilic_neg) %>% 
  dplyr::bind_rows(c18_pos) %>% 
  dplyr::bind_rows(c18_neg) %>% 
  dplyr::arrange(inchikey1, confidence_level, desc(msms_score_forward), rt_error)

writexl::write_xlsx(temp_merge, 
                    path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/merge_annotation_24626.xlsx', 
                    format_headers = FALSE)

# plot merged ms2 for manual check ---------------------------------------------

plot_merge_ms2(path = '~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/', 
               annot_merge_file = 'merge_annotation_24626.xlsx')

################################################################################
# Reduce the replication -------------------------------------------------------

# one feature - multiple metabolites 
result_manual <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/merge_annotation_manual_24626.xlsx', sheet = 2)

result_manual %>% 
  count(feature_name) %>% 
  filter(n > 1) %>% 
  pull(feature_name)


# one metabolite - multiple features
result_manual <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240524_B001_B033_data_processing/merge_annotation_manual_24626.xlsx', sheet = 3)

result_manual %>% 
  count(feature_name) %>% 
  filter(n > 1) %>% 
  pull(feature_name)

# one metabolite - multiple features
temp_inchikey <- result_manual %>% 
  count(inchikey1) %>% 
  filter(n > 1) %>% 
  pull(inchikey1)

result_manual %>%
  filter(inchikey1 %in% temp_inchikey)

