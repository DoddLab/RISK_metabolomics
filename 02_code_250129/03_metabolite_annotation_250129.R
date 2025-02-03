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
# Libraries match --------------------------------------------------------------
  # C18 pos --------------------------------------------------------------------
setwd('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/')


# library(DoddLabPackages)
library(DoddLabMetID)
library(DoddLabDatabase)

load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/01_input_data_cleaning/03_object_c18_pos_serrf.RData')

# mz + RT + ms2
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/',
                                                                  lib = 'dodd',
                                                                  column = 'c18',
                                                                  polarity = 'positive',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_c18_pos_serrf, parameter_set_annotation)

file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')



# mz + RT
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/01_input_data_cleaning/03_object_c18_pos_serrf.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/',
                                                                  lib = 'dodd',
                                                                  column = 'c18',
                                                                  polarity = 'positive',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = FALSE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_c18_pos_serrf, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt')

# mz + ms2
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/01_input_data_cleaning/03_object_c18_pos_serrf.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/',
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
annotate_metabolite(object = object_c18_pos_serrf, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_mz_ms2_public_db')



  # C18 neg ----------------------------------------------------------------------

# library(DoddLabPackages)
library(DoddLabMetID)
library(DoddLabDatabase)

setwd('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/')
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/01_input_data_cleaning/03_object_c18_neg_serrf.RData')

# mz + RT + ms2
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/',
                                                                  lib = 'dodd',
                                                                  column = 'c18',
                                                                  polarity = 'negative',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.8
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_c18_neg_serrf, parameter_set_annotation)

file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')



# mz + RT
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/01_input_data_cleaning/03_object_c18_neg_serrf.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/',
                                                                  lib = 'dodd',
                                                                  column = 'c18',
                                                                  polarity = 'negative',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = FALSE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
parameter_set_annotation@para_ms2_match$direction <- 'reverse'
annotate_metabolite(object = object_c18_neg_serrf, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt')

# mz + ms2
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/01_input_data_cleaning/03_object_c18_neg_serrf.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/',
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
annotate_metabolite(object = object_c18_neg_serrf, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_mz_ms2_public_db')


  # hilic pos ----------------------------------------------------------------------

# library(DoddLabPackages)
library(DoddLabMetID)
library(DoddLabDatabase)

setwd('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/')
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/01_input_data_cleaning/03_object_hilic_pos_serrf.RData')

# mz + RT + ms2
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/',
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
annotate_metabolite(object = object_hilic_pos_serrf, parameter_set_annotation)

file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')



# mz + RT
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/01_input_data_cleaning/03_object_hilic_pos_serrf.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/',
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
annotate_metabolite(object = object_hilic_pos_serrf, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt')

# mz + ms2
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/01_input_data_cleaning/03_object_hilic_pos_serrf.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/',
                                                                  lib = 'all_public',
                                                                  column = 'hilic',
                                                                  polarity = 'positive',
                                                                  is_rt_score = FALSE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.7
parameter_set_annotation@para_ms2_match$direction <- 'forward'
parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
annotate_metabolite(object = object_hilic_pos_serrf, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_mz_ms2_public_db')



  # hilic neg ----------------------------------------------------------------------

# library(DoddLabPackages)
library(DoddLabMetID)
library(DoddLabDatabase)

setwd('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/')
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/01_input_data_cleaning/03_object_hilic_neg_serrf.RData')

# mz + RT + ms2
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/',
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
annotate_metabolite(object = object_hilic_neg_serrf, parameter_set_annotation)

file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt_ms2')



# mz + RT
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/01_input_data_cleaning/03_object_hilic_neg_serrf.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(object = object_hilic_neg_serrf,
                                                                  path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/',
                                                                  lib = 'dodd',
                                                                  column = 'hilic',
                                                                  polarity = 'negative',
                                                                  is_rt_score = TRUE,
                                                                  is_ms2_score = FALSE)
parameter_set_annotation@para_ms1_match$mz_tol <- 10
parameter_set_annotation@para_ms1_match$tolerance_rt_range <- 30
parameter_set_annotation@para_ms1_match$mz_ppm_thr <- 200
annotate_metabolite(object = object_hilic_neg_serrf, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_dodd_mz_rt')

# mz + ms2
load('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/01_input_data_cleaning/03_object_hilic_neg_serrf.RData')
parameter_set_annotation <- initialize_annotation_parameter_class(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/',
                                                                  lib = 'all_public',
                                                                  column = 'hilic',
                                                                  polarity = 'negative',
                                                                  is_rt_score = FALSE,
                                                                  is_ms2_score = TRUE)
parameter_set_annotation@para_ms2_match$dp_cutoff <- 0.7
parameter_set_annotation@para_ms2_match$direction <- 'forward'
parameter_set_annotation@para_ms2_match$matched_frag_cutoff <- 2
annotate_metabolite(object = object_hilic_neg_serrf, parameter_set_annotation = parameter_set_annotation)
file.rename(from = './01_metabolite_annotation', to = './01_metabolite_annotation_mz_ms2_public_db')



################################################################################
# Merge annotation results for manual check ------------------------------------
library(tidyverse)
library(DoddLabMetID)

merge_one_modes(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_pos/',column = 'c18', polarity = 'positive')
merge_one_modes(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/c18_neg/',column = 'c18', polarity = 'negative')
merge_one_modes(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_pos/',column = 'hilic', polarity = 'positive')
merge_one_modes(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/hilic_neg/',column = 'hilic', polarity = 'negative')


hilic_pos <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/HILIC_pos/annotation_table_merge_hilic_positive.xlsx')
hilic_neg <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/HILIC_neg/annotation_table_merge_hilic_negative.xlsx')
c18_pos <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/C18_pos/annotation_table_merge_c18_positive.xlsx')
c18_neg <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/C18_neg/annotation_table_merge_c18_negative.xlsx')

temp_merge <- hilic_pos %>% 
  dplyr::bind_rows(hilic_neg) %>% 
  dplyr::bind_rows(c18_pos) %>% 
  dplyr::bind_rows(c18_neg) %>% 
  dplyr::arrange(inchikey1, confidence_level, desc(msms_score_forward), rt_error)

writexl::write_xlsx(temp_merge, 
                    path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/merge_annotation_240923.xlsx', 
                    format_headers = FALSE)

# plot merged ms2 for manual check ---------------------------------------------

plot_merge_ms2(path = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/', 
               annot_merge_file = 'merge_annotation_240923.xlsx')

################################################################################
# Manually check, dereplication, and modify the annotation result --------------
setwd('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/')

library(tidyverse)
library(tidymass)
library(DoddLabDatabase)

# load tidymass objects
# load('./hilic_pos/01_input_data_cleaning/04_object_hilic_pos_outlier_removal.RData')
load('~/Project/00_IBD_project/Data/20241010_manual_integration_leucine_isoleucine/04_object_hilic_pos_outlier_removal_241016.RData')
load('./hilic_neg/01_input_data_cleaning/04_object_hilic_neg_outlier_removal.RData')
load('./c18_pos/01_input_data_cleaning/04_object_c18_pos_outlier_removal.RData')
load('./c18_neg/01_input_data_cleaning/04_object_c18_neg_outlier_removal.RData')


annot_table_final <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/manual_check_annotation_table_with_class_info_manual_check_241218_2.xlsx')
annot_table_tidymass <- annot_table_final

# modify the annotation
annot_table_tidymass <- annot_table_tidymass %>% 
  dplyr::mutate(ms2_files_id = 'all_qc_dda', 
                ms2_spectrum_id = feature_name) %>% 
  dplyr::rename(variable_id = feature_name,
                Compound.name = name) %>% 
  dplyr::select(variable_id, ms2_files_id, ms2_spectrum_id, Compound.name, dplyr::everything())


# generate a new tidymass object for the statistics analysis 

# modify mass from each mode
  # hilic pos ------------------------------------------------------------------
expression_hilic_pos <- object_hilic_pos@expression_data
sample_info_hilic_pos <- object_hilic_pos@sample_info
variable_info_hilic_pos <- object_hilic_pos@variable_info
ms2_hilic_pos <- object_hilic_pos@ms2_data

# extract features in the final annot table
annot_table_hilic_pos <- annot_table_tidymass %>% 
  dplyr::filter(column == 'hilic' & polarity == 'positive') %>% 
  dplyr::arrange(mz, rt)
object_hilic_pos@annotation_table <- annot_table_hilic_pos

object_hilic_pos_annot <- object_hilic_pos %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  dplyr::filter(!is.na(confidence_level))


  # hilic neg ------------------------------------------------------------------
expression_hilic_neg <- object_hilic_neg@expression_data
sample_info_hilic_neg <- object_hilic_neg@sample_info
variable_info_hilic_neg <- object_hilic_neg@variable_info
ms2_hilic_neg <- object_hilic_neg@ms2_data

# extract features in the final annot table
annot_table_hilic_neg <- annot_table_tidymass %>% 
  dplyr::filter(column == 'hilic' & polarity == 'negative') %>% 
  dplyr::arrange(mz, rt)
object_hilic_neg@annotation_table <- annot_table_hilic_neg

object_hilic_neg_annot <- object_hilic_neg %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  dplyr::filter(!is.na(confidence_level))


  # c18 pos --------------------------------------------------------------------
expression_c18_pos <- object_c18_pos@expression_data
sample_info_c18_pos <- object_c18_pos@sample_info
variable_info_c18_pos <- object_c18_pos@variable_info
ms2_c18_pos <- object_c18_pos@ms2_data

# extract features in the final annot table
annot_table_c18_pos <- annot_table_tidymass %>% 
  dplyr::filter(column == 'c18' & polarity == 'positive') %>% 
  dplyr::arrange(mz, rt)
object_c18_pos@annotation_table <- annot_table_c18_pos

object_c18_pos_annot <- object_c18_pos %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  dplyr::filter(!is.na(confidence_level))


  # c18 neg --------------------------------------------------------------------
expression_c18_neg <- object_c18_neg@expression_data
sample_info_c18_neg <- object_c18_neg@sample_info
variable_info_c18_neg <- object_c18_neg@variable_info
ms2_c18_neg <- object_c18_neg@ms2_data

# extract features in the final annot table
annot_table_c18_neg <- annot_table_tidymass %>% 
  dplyr::filter(column == 'c18' & polarity == 'negative') %>% 
  dplyr::arrange(mz, rt)
object_c18_neg@annotation_table <- annot_table_c18_neg

object_c18_neg_annot <- object_c18_neg %>% 
  activate_mass_dataset(what = 'annotation_table') %>% 
  dplyr::filter(!is.na(confidence_level))


# Merge final data set -------------------------------------------------------
object_final <- merge_mass_dataset(x = object_hilic_pos_annot,
                                   y = object_hilic_neg_annot,
                                   sample_direction = "inner",
                                   variable_direction = "full", 
                                   sample_by = "sample_id", 
                                   variable_by = c("variable_id", "mz", "rt"))

object_final <- merge_mass_dataset(x = object_final,
                                   y = object_c18_pos_annot,
                                   sample_direction = "inner",
                                   variable_direction = "full", 
                                   sample_by = "sample_id", 
                                   variable_by = c("variable_id", "mz", "rt"))

object_final <- merge_mass_dataset(x = object_final,
                                   y = object_c18_neg_annot,
                                   sample_direction = "inner",
                                   variable_direction = "full", 
                                   sample_by = "sample_id", 
                                   variable_by = c("variable_id", "mz", "rt"))

colnames(object_final@sample_info)
object_final@sample_info <- object_final@sample_info %>% dplyr::select(sample_id:note)
object_final@variable_info <- object_final@variable_info %>% 
  dplyr::select(variable_id:rt)


object_final@sample_info_note <- object_final@sample_info_note %>% slice(1:21)
object_final@variable_info_note <- object_final@variable_info_note %>% slice(1:3)


save(object_final, 
     file = '~/Project/00_IBD_project/Data/20240919_IBD_B001_B044_analysis/00_manual_check_merge/object_final_241218.RData')