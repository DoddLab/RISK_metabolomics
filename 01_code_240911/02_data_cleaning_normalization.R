################################################################################
# Data cleaning and normalization  ---------------------------------------------

# Run all data cleaning on local computer, except SERRF normalization on Stanford SCG cluster

# Workflow:
#   1. create tidymass object by merging the peak area table, worklists.
#   2. applying the tidymass workflow to clean the data.
#     a. remove low-quality features (noise)
#     b. missing value imputation
#     c. normalization using SERRF

# DoddLabQC: v0.1.2; https://github.com/DoddLab/DoddLabQC
# DoddLabTool: v0.1.5; https://github.com/DoddLab/DoddLabTool




################################################################################
# Function: create_ibd_object ------------------------------------------------------------

create_ibd_object <- function(peak_table = 'Peak-table.csv',
                              merged_worklist = 'merged_worklist_c18_pos.xlsx',
                              dir_path = '.',
                              mode = c('c18_pos', 'c18_neg', 'hilic_pos', 'hilic_neg')) {
  # browser()
  mode <- match.arg(mode)
  
  cat('Modify the worklist for tidymass worklist ...\n')
  # modify the worklist format for object
  worklist_data <- readxl::read_xlsx(merged_worklist)
  sample_info <- worklist_data %>% 
    dplyr::filter(!(sample_type %in% c("blank", "poolQC_equilibrate", "procedureBlank", "poolQC_dda"))) %>% 
    dplyr::rename(sample_id = sample_name, 
                  injection.order = injection_order, 
                  class = sample_type, 
                  batch = batch_id, 
                  subject.id = deidentified_master_patient_id,
                  group = sample_type2) %>% 
    dplyr::select(sample_id, injection.order, class, batch, subject.id, group, everything()) %>% 
    dplyr::mutate(batch = stringr::str_extract(batch, pattern = '\\d+') %>% as.numeric()) %>% 
    dplyr::mutate(sample_id = stringr::str_replace(sample_id, pattern = '\\-', replacement = '_'))  %>% 
    dplyr::group_by(batch) %>% 
    dplyr::mutate(injection.order = seq(n())) %>% 
    dplyr::ungroup()
  
  # modify the peak table
  cat('Modify the expression_data & variable_info ...\n')
  sample_id <- sample_info$sample_id
  peak_table <- readr::read_csv(file.path(dir_path, peak_table), show_col_types = FALSE)
  colnames(peak_table) <- colnames(peak_table) %>% 
    stringr::str_replace(pattern = '-', replacement = '_')
  
  peak_table <- peak_table %>% 
    dplyr::select(name, mzmed, rtmed, sample_id) %>% 
    dplyr::rename(mz = mzmed,
                  rt = rtmed,
                  variable_id = name) %>% 
    dplyr::mutate(variable_id = paste0(variable_id, '_', mode))
  
  # expression_data_pos
  expression_data <- peak_table %>% 
    dplyr::select(-c(variable_id:rt)) %>% 
    as.data.frame()
  
  # variable_info_pos
  variable_info <- peak_table %>% 
    dplyr::select(variable_id:rt) %>% 
    as.data.frame()
  
  rownames(expression_data) <- variable_info$variable_id
  
  # check sample_id in sample_info whether same with expression_data, if not, reorganize
  colnames(expression_data) == sample_info$sample_id
  expression_data <- expression_data[, sample_info$sample_id]
  
  # creat mass_data object -------------------------------------------------------
  cat('Create the mass dataset object ...\n')
  object <- create_mass_dataset(expression_data = expression_data, 
                                sample_info = sample_info, 
                                variable_info = variable_info)
  
  dir.create(path = file.path(dir_path, '01_input_data_cleaning'), showWarnings = FALSE, recursive = TRUE)
  save(object, 
       file = file.path(dir_path, '/01_input_data_cleaning/object.RData'))
  
  export_mass_dataset(object = object,
                      file_type = "xlsx",
                      path = file.path(dir_path, "/01_input_data_cleaning/mass_data_tables"))
  
  cat('Done!\n')
}



################################################################################
# load packages ----------------------------------------------------------------
library(tidyverse)
library(tidymass)
library(sjmisc)
library(DoddLabQC)
library(DoddLabTool)

################################################################################
# C18 positive -----------------------------------------------------------------
  # 1. Create tidymass object --------------------------------------------------
setwd('~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_pos/')

create_ibd_object(peak_table = 'Peak-table.csv',
                  merged_worklist = '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/merged_worklist/merged_worklist_c18_pos.xlsx',
                  dir_path = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_pos/00_raw_data_processing/',
                  mode = 'c18_pos')



  # 2. Data cleaning ------------------------------------------------------------
load('./00_raw_data_processing/01_input_data_cleaning/object.RData')

object_c18_pos <- object
rm(object);gc()

# plot global PCA for all batches & samples before cleaning
temp_data_trans <- object_c18_pos@expression_data %>% 
  rotate_df()
pca_pos <- pca(temp_data_trans, center = TRUE, scale = TRUE) 


plotIndiv(pca_pos,
          group = object_c18_pos@sample_info$batch,
          ellipse = FALSE,
          legend = TRUE,
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = object_c18_pos@sample_info$sample_id,
          title = 'Before cleaning')

  #   2.a remove low-quality features ------------------------------------------

# object: replace 0 value as NA
object_c18_pos@expression_data <- object_c18_pos@expression_data %>% 
  mutate(across(where(is.numeric), ~na_if(., 0)))

# criteria of removal:
# MV number >= 80% pool QC
# MV Non-IBD >= 50%
# IBD (B1->B1, B1->B2, B1->B3) >= 50%, either group

pool_qc_id <- object_c18_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "poolQC") %>%
  pull(sample_id)

non_IDB_id <- object_c18_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == "Non-IBD") %>%
  pull(sample_id)

B1_id <- object_c18_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B1') %>%
  pull(sample_id)

B2_id <- object_c18_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B2') %>%
  pull(sample_id)

B3_id <- object_c18_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B3') %>%
  pull(sample_id)

object_c18_pos <- object_c18_pos %>%
  mutate_variable_na_freq(according_to_samples = pool_qc_id) %>% 
  mutate_variable_na_freq(according_to_samples = non_IDB_id) %>% 
  mutate_variable_na_freq(according_to_samples = B1_id) %>% 
  mutate_variable_na_freq(according_to_samples = B2_id) %>% 
  mutate_variable_na_freq(according_to_samples = B3_id)

# extract_variable_info(object_c18_pos) %>% head()
object_c18_pos <- 
  object_c18_pos %>% 
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2 & (na_freq.1 < 0.5 | na_freq.2 < 0.5 | na_freq.3 < 0.5 | na_freq.4 < 0.5))

dir.create('./00_raw_data_processing/01_input_data_cleaning', showWarnings = FALSE, recursive = TRUE)
save(object_c18_pos, file = './00_raw_data_processing/01_input_data_cleaning/01_object_c18_pos_noise_removal_24829.RData')

object_c18_pos@expression_data %>% dim()



  #   2.b missing value imputation ---------------------------------------------
get_mv_number(object_c18_pos)

object_c18_pos <- 
  impute_mv(object = object_c18_pos, method = "knn")

# check whether introduce zero values by KNN. If it exist, remove them
temp <- which(object_c18_pos@expression_data < 0.0001, arr.ind = TRUE)
if (length(temp) > 0) {
  peak_removal <- which(object_c18_pos@expression_data < 0.0001, arr.ind = TRUE) %>% rownames() %>% unique()
  temp_idx <- match(peak_removal, rownames(object_c18_pos@expression_data))
  object_c18_pos@expression_data[temp_idx,] %>% 
    rotate_df() %>% 
    rownames_to_column() %>% 
    tidyr::pivot_longer(cols = -rowname) %>% 
    dplyr::group_by(name) %>% 
    summarise(sum(value < 0.0001)) 
  
  object_c18_pos <- object_c18_pos %>% 
    activate_mass_dataset(what = 'variable_info') %>% 
    filter(!(variable_id %in% peak_removal))
}


object_c18_pos@sample_info <- object_c18_pos@sample_info %>% 
  mutate(batch = as.character(batch))

get_mv_number(object_c18_pos)
show_missing_values(object = object_c18_pos, show_column_names = FALSE, percentage = TRUE)

save(object_c18_pos, 
     file = './00_raw_data_processing/01_input_data_cleaning/02_object_c18_pos_imputation.RData')


  #   2.c normalization --------------------------------------------------------------

# evaluate before normalization 
load('./00_raw_data_processing/01_input_data_cleaning/02_object_c18_pos_imputation.RData')
object_c18_pos <- object_c18_pos %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(class != 'nistQC') %>% 
  mutate(class = case_when(class == 'poolQC' ~ 'QC',
                           class == 'sample' ~ 'Subject')) %>% 
  mutate(injection.order = seq(n()))

# original PCA plot
temp_object <- object_c18_pos
temp_data_trans <- temp_object@expression_data %>% 
  rotate_df()

pca_pos <- pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$class
pca_batch <- temp_object@sample_info$batch
pca_var_id <- temp_object@sample_info$sample_id

dir.create('./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots', showWarnings = FALSE, recursive = TRUE)
pdf(file = './00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/original_PCA.pdf', 
    width = 8, height = 6, onefile = TRUE)

plotIndiv(pca_pos, 
          group = pca_group, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

dev.off()

# show ISTD intensity plot 
istd_peak <- match_istd(object = object_c18_pos, mode = 'c18_pos', lib_version = 'v2', mz_ppm = 10, rt_second = 20)
order_injection <- object_c18_pos@sample_info$sample_id %>% stringr::str_detect('PoolQC06') %>% which()

temp_plot <- lapply(seq_along(istd_peak$variable_id), function(i){
  x <- istd_peak$variable_id[i]
  temp_data <- temp_object %>% 
    intensity_plot(variable_id = x,
                   color_by = 'class',
                   order_by = 'injection.order', 
                   interactive = FALSE)
  
  temp_plot <- temp_data$data %>% 
    arrange(batch, injection.order) %>% 
    mutate(x = seq(n())) %>%
    ggplot(aes(x = x, y = int, color = class)) +
    geom_point() +
    scale_x_discrete(labels = NULL) +
    ZZWtool::ZZWTheme() +
    theme(axis.text.y.left = element_text(hjust = 0.5)) +
    ggtitle(paste(istd_peak$name_istd[i])) +
    geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed') +
    xlab('Injection Order')
})

ggsave(
  filename = "./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/ISTD_plots_before_intrabatch_normalization_2.pdf", 
  plot = gridExtra::marrangeGrob(temp_plot, nrow=3, ncol=3), 
  width = 15, height = 9
)


# normalization via the SERRF
#  convert data for serrf normalization, and run on Stanford SCG cluster

load('./00_raw_data_processing/01_input_data_cleaning/02_object_c18_pos_imputation.RData')

library(tidyverse)
library(tidymass)
library(sjmisc)
library(DoddLabQC)
library(DoddLabTool)
library(DoddLabSERRF)


# convert data for serrf
table_serrf <- convert_object_for_serrf(object = object_c18_pos)

dir.create('~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_pos/00_raw_data_processing/01_input_data_cleaning/01_serrf', 
           showWarnings = FALSE, recursive = TRUE)

saveRDS(table_serrf, 
        file = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_pos/00_raw_data_processing/01_input_data_cleaning/01_serrf/raw_data_240829.RDS')

# run serrf normalization
normalize_serrf(data_file = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_pos/00_raw_data_processing/01_input_data_cleaning/01_serrf/raw_data_240829.RDS',
                dir_path = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_pos/00_raw_data_processing/01_input_data_cleaning/01_serrf/')


# merge SERRF result and update the tidymass object 
setwd('~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_pos/')

load('./00_raw_data_processing/01_input_data_cleaning/02_object_c18_pos_imputation.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/serrf_result/serrf_result.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/p.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/f.RData')

object_c18_pos_serrf <- merge_serrf_result(object = object_c18_pos, 
                                           serrf_result = serrf_result, 
                                           p = p,
                                           f = f)

save(object_c18_pos_serrf,
     file = './00_raw_data_processing/01_input_data_cleaning/03_object_c18_pos_serrf.RData')

rm(list = ls())


# evaluate the SERRF normalization result
temp_object <- object_c18_pos_serrf %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class %in% c('sample', 'poolQC'))


temp_data_trans <- temp_object@expression_data %>% 
  # activate_mass_dataset(what = 'expression_data') %>%
  # mutate(across(everything(), as.numeric)) %>%
  rotate_df()

pca_pos <- mixOmics::pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$class
pca_group2 <- temp_object@sample_info$group
pca_batch <- temp_object@sample_info$batch
pca_var_id <- temp_object@sample_info$sample_id

dir.create('./00_raw_data_processing/01_input_data_cleaning/plots', showWarnings = FALSE, recursive = TRUE)
pdf(file = './00_raw_data_processing/01_input_data_cleaning/plots/serrf_normalized_PCA.pdf', 
    width = 8, height = 6, onefile = TRUE)

plotIndiv(pca_pos, 
          group = pca_group, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_group2, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

dev.off()



# show ISTD intensity plot 
istd_peak <- match_istd(object = temp_object, 
                        mode = 'c18_pos', 
                        lib_version = 'v2', 
                        mz_ppm = 10, 
                        rt_second = 30)

order_injection <- temp_object@sample_info$sample_id %>% stringr::str_detect('PoolQC06') %>% which()

temp_plot <- lapply(seq_along(istd_peak$variable_id), function(i){
  x <- istd_peak$variable_id[i]
  temp_data <- temp_object %>% 
    intensity_plot(variable_id = x,
                   color_by = 'class',
                   order_by = 'injection.order', 
                   interactive = FALSE)
  
  temp_plot <- temp_data$data %>% 
    arrange(batch, injection.order) %>% 
    mutate(x = seq(n())) %>%
    ggplot(aes(x = x, y = int, color = class)) +
    geom_point() +
    scale_x_discrete(labels = NULL) +
    ZZWtool::ZZWTheme() +
    theme(axis.text.y.left = element_text(hjust = 0.5)) +
    ggtitle(paste(istd_peak$name_istd[i])) +
    geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed') +
    xlab('Injection Order')
})

ggsave(
  filename = "./00_raw_data_processing/01_input_data_cleaning/plots/serrf_normalized_ISTD_plots.pdf", 
  plot = gridExtra::marrangeGrob(temp_plot, nrow=3, ncol=3), 
  width = 15, height = 9
)


# RSD calculation 
temp_qc_rsd <- temp_object %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class == 'poolQC') %>% 
  extract_expression_data() %>%
  apply(., 1, function(x){sd(x)/mean(x)*100})

temp_all_rsd <- temp_object %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class != 'poolQC') %>%
  extract_expression_data() %>%
  apply(., 1, function(x){sd(x)/mean(x)*100})

temp_object <- temp_object %>%
  activate_mass_dataset(what = 'variable_info') %>%
  mutate(rsd_qc = temp_qc_rsd,
         rsd_all_sample = temp_all_rsd)

# all feature
temp_object %>%
  extract_variable_info() %>%
  dplyr::summarise(mean_qc_rsd = mean(rsd_qc, na.rm = TRUE),
                   mean_all_sample_rsd = mean(rsd_all_sample, na.rm = TRUE),
                   qc_rsd_less_30 = sum(rsd_qc <= 30)/n(),
                   all_sample_rsd_less_30 = sum(rsd_all_sample <= 30)/n())

# ISTD  
temp_object %>%
  extract_variable_info() %>%
  dplyr::filter(variable_id %in% istd_peak$variable_id) %>%
  dplyr::summarise(mean_qc_rsd = mean(rsd_qc, na.rm = TRUE),
                   mean_all_sample_rsd = mean(rsd_all_sample, na.rm = TRUE),
                   qc_rsd_less_30 = sum(rsd_qc <= 30)/n(),
                   all_sample_rsd_less_30 = sum(rsd_all_sample <= 30)/n())

rm(list = ls())


################################################################################
# C18 negative -----------------------------------------------------------------
  # 1. Create tidymass object --------------------------------------------------
rm(list = ls());gc()

setwd('~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_neg/')

create_ibd_object(peak_table = 'Peak-table.csv',
                  merged_worklist = '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/merged_worklist/merged_worklist_c18_neg.xlsx',
                  dir_path = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_neg/00_raw_data_processing/',
                  mode = 'c18_neg')



  # 2. Data cleaning -----------------------------------------------------------
load('./00_raw_data_processing/01_input_data_cleaning/object.RData')

object_c18_neg <- object
rm(object);gc()

# plot global PCA for all batches & samples before cleaning
temp_data_trans <- object_c18_neg@expression_data %>% 
  rotate_df()
pca_pos <- pca(temp_data_trans, center = TRUE, scale = TRUE) 


plotIndiv(pca_pos,
          group = object_c18_neg@sample_info$batch,
          ellipse = FALSE,
          legend = TRUE,
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = object_c18_neg@sample_info$sample_id,
          title = 'Before cleaning')

  #   2.a remove low-quality features ------------------------------------------

# object: replace 0 value as NA
object_c18_neg@expression_data <- object_c18_neg@expression_data %>% 
  mutate(across(where(is.numeric), ~na_if(., 0)))

# criteria of removal:
# MV number >= 80% pool QC
# MV Non-IBD >= 50%
# IBD (B1->B1, B1->B2, B1->B3) >= 50%, either group

pool_qc_id <- object_c18_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "poolQC") %>%
  pull(sample_id)

non_IDB_id <- object_c18_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == "Non-IBD") %>%
  pull(sample_id)

B1_id <- object_c18_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B1') %>%
  pull(sample_id)

B2_id <- object_c18_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B2') %>%
  pull(sample_id)

B3_id <- object_c18_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B3') %>%
  pull(sample_id)

object_c18_neg <- object_c18_neg %>%
  mutate_variable_na_freq(according_to_samples = pool_qc_id) %>% 
  mutate_variable_na_freq(according_to_samples = non_IDB_id) %>% 
  mutate_variable_na_freq(according_to_samples = B1_id) %>% 
  mutate_variable_na_freq(according_to_samples = B2_id) %>% 
  mutate_variable_na_freq(according_to_samples = B3_id)

# extract_variable_info(object_c18_neg) %>% head()
object_c18_neg <- 
  object_c18_neg %>% 
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2 & (na_freq.1 < 0.5 | na_freq.2 < 0.5 | na_freq.3 < 0.5 | na_freq.4 < 0.5))

dir.create('./00_raw_data_processing/01_input_data_cleaning', showWarnings = FALSE, recursive = TRUE)
save(object_c18_neg, file = './00_raw_data_processing/01_input_data_cleaning/01_object_c18_neg_noise_removal_240829.RData')

object_c18_neg@expression_data %>% dim()



  #   2.b missing value imputation ---------------------------------------------
get_mv_number(object_c18_neg)

object_c18_neg <- 
  impute_mv(object = object_c18_neg, method = "knn")

# check whether introduce zero values by KNN. If it exist, remove them
temp <- which(object_c18_neg@expression_data < 0.0001, arr.ind = TRUE)
if (length(temp) > 0) {
  peak_removal <- which(object_c18_neg@expression_data < 0.0001, arr.ind = TRUE) %>% rownames() %>% unique()
  temp_idx <- match(peak_removal, rownames(object_c18_neg@expression_data))
  object_c18_neg@expression_data[temp_idx,] %>% 
    rotate_df() %>% 
    rownames_to_column() %>% 
    tidyr::pivot_longer(cols = -rowname) %>% 
    dplyr::group_by(name) %>% 
    summarise(sum(value < 0.0001)) 
  
  object_c18_neg <- object_c18_neg %>% 
    activate_mass_dataset(what = 'variable_info') %>% 
    filter(!(variable_id %in% peak_removal))
}


object_c18_neg@sample_info <- object_c18_neg@sample_info %>% 
  mutate(batch = as.character(batch))

get_mv_number(object_c18_neg)
show_missing_values(object = object_c18_neg, show_column_names = FALSE, percentage = TRUE)

save(object_c18_neg, 
     file = './00_raw_data_processing/01_input_data_cleaning/02_object_c18_neg_imputation.RData')


  #   2.c normalization --------------------------------------------------------
# evaluate before normalization 
load('./00_raw_data_processing/01_input_data_cleaning/02_object_c18_neg_imputation.RData')
object_c18_neg <- object_c18_neg %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(class != 'nistQC') %>% 
  mutate(class = case_when(class == 'poolQC' ~ 'QC',
                           class == 'sample' ~ 'Subject')) %>% 
  mutate(injection.order = seq(n()))

# original PCA plot
temp_object <- object_c18_neg
temp_data_trans <- temp_object@expression_data %>% 
  rotate_df()

pca_pos <- pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$class
pca_batch <- temp_object@sample_info$batch
pca_var_id <- temp_object@sample_info$sample_id

dir.create('./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots', showWarnings = FALSE, recursive = TRUE)
pdf(file = './00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/original_PCA.pdf', 
    width = 8, height = 6, onefile = TRUE)

plotIndiv(pca_pos, 
          group = pca_group, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

dev.off()

# show ISTD intensity plot 
istd_peak <- match_istd(object = object_c18_neg, mode = 'c18_neg', lib_version = 'v2', mz_ppm = 10, rt_second = 15)
order_injection <- object_c18_neg@sample_info$sample_id %>% stringr::str_detect('PoolQC06') %>% which()



temp_plot <- lapply(seq_along(istd_peak$variable_id), function(i){
  x <- istd_peak$variable_id[i]
  temp_data <- temp_object %>% 
    intensity_plot(variable_id = x,
                   color_by = 'class',
                   order_by = 'injection.order', 
                   interactive = FALSE)
  
  temp_plot <- temp_data$data %>% 
    arrange(batch, injection.order) %>% 
    mutate(x = seq(n())) %>%
    ggplot(aes(x = x, y = int, color = class)) +
    geom_point() +
    scale_x_discrete(labels = NULL) +
    ZZWtool::ZZWTheme() +
    theme(axis.text.y.left = element_text(hjust = 0.5)) +
    ggtitle(paste(istd_peak$name_istd[i])) +
    geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed') +
    xlab('Injection Order')
})

ggsave(
  filename = "./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/ISTD_plots_before_intrabatch_normalization_2.pdf", 
  plot = gridExtra::marrangeGrob(temp_plot, nrow=4, ncol=3), 
  width = 15, height = 9
)



# normalization via the SERRF 

load('./00_raw_data_processing/01_input_data_cleaning/02_object_c18_neg_imputation.RData')

library(DoddLabSERRF)


# convert data for serrf
table_serrf <- convert_object_for_serrf(object = object_c18_neg)

dir.create('./00_raw_data_processing/01_input_data_cleaning/01_serrf', 
           showWarnings = FALSE, recursive = TRUE)

saveRDS(table_serrf, 
        file = './00_raw_data_processing/01_input_data_cleaning/01_serrf/raw_data_240829.RDS')

# run serrf normalization
normalize_serrf(data_file = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_pos/00_raw_data_processing/01_input_data_cleaning/01_serrf/raw_data_240829.RDS',
                dir_path = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_pos/00_raw_data_processing/01_input_data_cleaning/01_serrf/')


# merge SERRF result and update the tidymass object
setwd('~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/c18_neg/')

library(tidyverse)
library(tidymass)
library(sjmisc)
library(DoddLabQC)
library(DoddLabTool)
library(DoddLabSERRF)

load('./00_raw_data_processing/01_input_data_cleaning/02_object_c18_neg_imputation.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/serrf_result/serrf_result.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/p.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/f.RData')

object_c18_neg_serrf <- merge_serrf_result(object = object_c18_neg, 
                                           serrf_result = serrf_result, 
                                           p = p,
                                           f = f)

save(object_c18_neg_serrf,
     file = './00_raw_data_processing/01_input_data_cleaning/03_object_c18_neg_serrf.RData')

rm(list = c('object_c18_neg', 'serrf_result', 'p', 'f'));gc()


# evaluate the SERRF normalization result
temp_object <- object_c18_neg_serrf %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class %in% c('sample', 'poolQC'))


temp_data_trans <- temp_object@expression_data %>% 
  # activate_mass_dataset(what = 'expression_data') %>%
  # mutate(across(everything(), as.numeric)) %>%
  rotate_df()

pca_pos <- mixOmics::pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$class
pca_group2 <- temp_object@sample_info$group
pca_batch <- temp_object@sample_info$batch
pca_var_id <- temp_object@sample_info$sample_id

dir.create('./00_raw_data_processing/01_input_data_cleaning/plots', showWarnings = FALSE, recursive = TRUE)
pdf(file = './00_raw_data_processing/01_input_data_cleaning/plots/serrf_normalized_PCA.pdf', 
    width = 8, height = 6, onefile = TRUE)

plotIndiv(pca_pos, 
          group = pca_group, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_group2, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

dev.off()



# show ISTD intensity plot 
istd_peak <- match_istd(object = temp_object, 
                        mode = 'c18_neg', 
                        lib_version = 'v2', 
                        mz_ppm = 10, 
                        rt_second = 30)

order_injection <- temp_object@sample_info$sample_id %>% stringr::str_detect('PoolQC06') %>% which()

temp_plot <- lapply(seq_along(istd_peak$variable_id), function(i){
  x <- istd_peak$variable_id[i]
  temp_data <- temp_object %>% 
    intensity_plot(variable_id = x,
                   color_by = 'class',
                   order_by = 'injection.order', 
                   interactive = FALSE)
  
  temp_plot <- temp_data$data %>% 
    arrange(batch, injection.order) %>% 
    mutate(x = seq(n())) %>%
    ggplot(aes(x = x, y = int, color = class)) +
    geom_point() +
    scale_x_discrete(labels = NULL) +
    ZZWtool::ZZWTheme() +
    theme(axis.text.y.left = element_text(hjust = 0.5)) +
    ggtitle(paste(istd_peak$name_istd[i])) +
    geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed') +
    xlab('Injection Order')
})

ggsave(
  filename = "./00_raw_data_processing/01_input_data_cleaning/plots/serrf_normalized_ISTD_plots.pdf", 
  plot = gridExtra::marrangeGrob(temp_plot, nrow=3, ncol=3), 
  width = 15, height = 9
)


# RSD calculation 
temp_qc_rsd <- temp_object %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class == 'poolQC') %>% 
  extract_expression_data() %>%
  apply(., 1, function(x){sd(x)/mean(x)*100})

temp_all_rsd <- temp_object %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class != 'poolQC') %>%
  extract_expression_data() %>%
  apply(., 1, function(x){sd(x)/mean(x)*100})

temp_object <- temp_object %>%
  activate_mass_dataset(what = 'variable_info') %>%
  mutate(rsd_qc = temp_qc_rsd,
         rsd_all_sample = temp_all_rsd)

# all feature
temp_object %>%
  extract_variable_info() %>%
  dplyr::summarise(mean_qc_rsd = mean(rsd_qc, na.rm = TRUE),
                   mean_all_sample_rsd = mean(rsd_all_sample, na.rm = TRUE),
                   qc_rsd_less_30 = sum(rsd_qc <= 30)/n(),
                   all_sample_rsd_less_30 = sum(rsd_all_sample <= 30)/n())

# ISTD  
temp_object %>%
  extract_variable_info() %>%
  dplyr::filter(variable_id %in% istd_peak$variable_id) %>%
  dplyr::summarise(mean_qc_rsd = mean(rsd_qc, na.rm = TRUE),
                   mean_all_sample_rsd = mean(rsd_all_sample, na.rm = TRUE),
                   qc_rsd_less_30 = sum(rsd_qc <= 30)/n(),
                   all_sample_rsd_less_30 = sum(rsd_all_sample <= 30)/n())

rm(list = ls())



################################################################################
# HILIC positive ---------------------------------------------------------------
rm(list = ls());gc()

  # 1. Create tidymass object --------------------------------------------------
setwd('~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_pos/')

create_ibd_object(peak_table = 'Peak-table.csv',
                  merged_worklist = '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/merged_worklist/merged_worklist_hilic_pos.xlsx',
                  dir_path = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_pos/00_raw_data_processing/',
                  mode = 'hilic_pos')



  # 2. Data cleaning -----------------------------------------------------------
load('./00_raw_data_processing/01_input_data_cleaning/object.RData')

object_hilic_pos <- object
rm(object);gc()

# plot global PCA for all batches & samples before cleaning
temp_data_trans <- object_hilic_pos@expression_data %>% 
  rotate_df()
pca_pos <- pca(temp_data_trans, center = TRUE, scale = TRUE) 


plotIndiv(pca_pos,
          group = object_hilic_pos@sample_info$batch,
          ellipse = FALSE,
          legend = TRUE,
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = object_hilic_pos@sample_info$sample_id,
          title = 'Before cleaning')

  #   2.a remove low-quality features ------------------------------------------

# object: replace 0 value as NA
object_hilic_pos@expression_data <- object_hilic_pos@expression_data %>% 
  mutate(across(where(is.numeric), ~na_if(., 0)))

# criteria of removal:
# MV number >= 80% pool QC
# MV Non-IBD >= 50%
# IBD (B1->B1, B1->B2, B1->B3) >= 50%, either group

pool_qc_id <- object_hilic_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "poolQC") %>%
  pull(sample_id)

non_IDB_id <- object_hilic_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == "Non-IBD") %>%
  pull(sample_id)

B1_id <- object_hilic_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B1') %>%
  pull(sample_id)

B2_id <- object_hilic_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B2') %>%
  pull(sample_id)

B3_id <- object_hilic_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B3') %>%
  pull(sample_id)

object_hilic_pos <- object_hilic_pos %>%
  mutate_variable_na_freq(according_to_samples = pool_qc_id) %>% 
  mutate_variable_na_freq(according_to_samples = non_IDB_id) %>% 
  mutate_variable_na_freq(according_to_samples = B1_id) %>% 
  mutate_variable_na_freq(according_to_samples = B2_id) %>% 
  mutate_variable_na_freq(according_to_samples = B3_id)

# extract_variable_info(object_hilic_pos) %>% head()
object_hilic_pos <- 
  object_hilic_pos %>% 
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2 & (na_freq.1 < 0.5 | na_freq.2 < 0.5 | na_freq.3 < 0.5 | na_freq.4 < 0.5))

dir.create('./00_raw_data_processing/01_input_data_cleaning', showWarnings = FALSE, recursive = TRUE)
save(object_hilic_pos, file = './00_raw_data_processing/01_input_data_cleaning/01_object_hilic_pos_noise_removal_240829.RData')

object_hilic_pos@expression_data %>% dim()

# # show missing value plot again
# show_missing_values(object = object_hilic_pos, show_column_names = FALSE, percentage = TRUE)


  #   2.b missing value imputation ---------------------------------------------
get_mv_number(object_hilic_pos)

object_hilic_pos <- 
  impute_mv(object = object_hilic_pos, method = "knn")

# check whether introduce zero values by KNN. If it exist, remove them
temp <- which(object_hilic_pos@expression_data < 0.0001, arr.ind = TRUE)
if (length(temp) > 0) {
  peak_removal <- which(object_hilic_pos@expression_data < 0.0001, arr.ind = TRUE) %>% rownames() %>% unique()
  temp_idx <- match(peak_removal, rownames(object_hilic_pos@expression_data))
  object_hilic_pos@expression_data[temp_idx,] %>% 
    rotate_df() %>% 
    rownames_to_column() %>% 
    tidyr::pivot_longer(cols = -rowname) %>% 
    dplyr::group_by(name) %>% 
    summarise(sum(value < 0.0001)) 
  
  object_hilic_pos <- object_hilic_pos %>% 
    activate_mass_dataset(what = 'variable_info') %>% 
    filter(!(variable_id %in% peak_removal))
}


object_hilic_pos@sample_info <- object_hilic_pos@sample_info %>% 
  mutate(batch = as.character(batch))

get_mv_number(object_hilic_pos)
show_missing_values(object = object_hilic_pos, show_column_names = FALSE, percentage = TRUE)

save(object_hilic_pos, 
     file = './00_raw_data_processing/01_input_data_cleaning/02_object_hilic_pos_imputation.RData')


  #   2.c normalization --------------------------------------------------------
# before normalization 
load('./00_raw_data_processing/01_input_data_cleaning/02_object_hilic_pos_imputation.RData')
object_hilic_pos <- object_hilic_pos %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(class != 'nistQC') %>% 
  mutate(class = case_when(class == 'poolQC' ~ 'QC',
                           class == 'sample' ~ 'Subject')) %>% 
  mutate(injection.order = seq(n()))

# original PCA plot
temp_object <- object_hilic_pos
temp_data_trans <- temp_object@expression_data %>% 
  rotate_df()

pca_pos <- pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$class
pca_batch <- temp_object@sample_info$batch
pca_var_id <- temp_object@sample_info$sample_id

dir.create('./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots', showWarnings = FALSE, recursive = TRUE)
pdf(file = './00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/original_PCA.pdf', 
    width = 8, height = 6, onefile = TRUE)

plotIndiv(pca_pos, 
          group = pca_group, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

dev.off()

# show ISTD intensity plot 
istd_peak <- match_istd(object = object_hilic_pos, mode = 'hilic_pos', lib_version = 'v2', mz_ppm = 10, rt_second = 30)
order_injection <- object_hilic_pos@sample_info$sample_id %>% stringr::str_detect('PoolQC06') %>% which()

# temp_plot <- lapply(istd_peak$variable_id, function(x){
#   object_hilic_pos %>% 
#     intensity_plot(variable_id = x,
#                    color_by = 'class',
#                    order_by = 'injection.order', 
#                    interactive = FALSE) + 
#     scale_x_discrete(labels = NULL) +
#     ZZWtool::ZZWTheme() +
#     theme(axis.text.y.left = element_text(hjust = 0.5)) +
#     ggtitle(paste(x, 'Before intra-batch normalization')) +
#     geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed')
# })
# 
# ggsave(
#   filename = "./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/ISTD_plots_before_intrabatch_normalization.pdf", 
#   plot = gridExtra::marrangeGrob(temp_plot, nrow=1, ncol=1), 
#   width = 15, height = 9
# )

temp_plot <- lapply(seq_along(istd_peak$variable_id), function(i){
  x <- istd_peak$variable_id[i]
  temp_data <- temp_object %>% 
    intensity_plot(variable_id = x,
                   color_by = 'class',
                   order_by = 'injection.order', 
                   interactive = FALSE)
  
  temp_plot <- temp_data$data %>% 
    arrange(batch, injection.order) %>% 
    mutate(x = seq(n())) %>%
    ggplot(aes(x = x, y = int, color = class)) +
    geom_point() +
    scale_x_discrete(labels = NULL) +
    ZZWtool::ZZWTheme() +
    theme(axis.text.y.left = element_text(hjust = 0.5)) +
    ggtitle(paste(istd_peak$name_istd[i])) +
    geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed') +
    xlab('Injection Order')
})

ggsave(
  filename = "./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/ISTD_plots_before_intrabatch_normalization_2.pdf", 
  plot = gridExtra::marrangeGrob(temp_plot, nrow=3, ncol=3), 
  width = 15, height = 9
)


# normalization via the SERRF

load('./00_raw_data_processing/01_input_data_cleaning/02_object_hilic_pos_imputation.RData')

library(tidyverse)
library(tidymass)
library(sjmisc)
library(DoddLabQC)
library(DoddLabTool)
library(DoddLabSERRF)


# convert data for serrf
table_serrf <- convert_object_for_serrf(object = object_hilic_pos)

dir.create('./00_raw_data_processing/01_input_data_cleaning/01_serrf', 
           showWarnings = FALSE, recursive = TRUE)

saveRDS(table_serrf, 
        file = './00_raw_data_processing/01_input_data_cleaning/01_serrf/raw_data_240829.RDS')

# run serrf normalization
normalize_serrf(data_file = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_pos/00_raw_data_processing/01_input_data_cleaning/01_serrf/raw_data_240829.RDS',
                dir_path = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_pos/00_raw_data_processing/01_input_data_cleaning/01_serrf/')

# merge SERRF result and update the tidymass object 
setwd('~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_pos/')

library(tidyverse)
library(tidymass)
library(sjmisc)
library(DoddLabQC)
library(DoddLabTool)
library(DoddLabSERRF)

load('./00_raw_data_processing/01_input_data_cleaning/02_object_hilic_pos_imputation.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/serrf_result/serrf_result.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/p.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/f.RData')

object_hilic_pos_serrf <- merge_serrf_result(object = object_hilic_pos, 
                                             serrf_result = serrf_result, 
                                             p = p,
                                             f = f)

save(object_hilic_pos_serrf,
     file = './00_raw_data_processing/01_input_data_cleaning/03_object_hilic_pos_serrf.RData')

rm(list = c('object_hilic_pos', 'serrf_result', 'p', 'f'));gc()


# evaluate the SERRF normalization result
temp_object <- object_hilic_pos_serrf %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class %in% c('sample', 'poolQC'))


temp_data_trans <- temp_object@expression_data %>% 
  # activate_mass_dataset(what = 'expression_data') %>%
  # mutate(across(everything(), as.numeric)) %>%
  rotate_df()

pca_pos <- mixOmics::pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$class
pca_group2 <- temp_object@sample_info$group
pca_batch <- temp_object@sample_info$batch
pca_var_id <- temp_object@sample_info$sample_id

dir.create('./00_raw_data_processing/01_input_data_cleaning/plots', showWarnings = FALSE, recursive = TRUE)
pdf(file = './00_raw_data_processing/01_input_data_cleaning/plots/serrf_normalized_PCA.pdf', 
    width = 8, height = 6, onefile = TRUE)

plotIndiv(pca_pos, 
          group = pca_group, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_group2, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

dev.off()



# show ISTD intensity plot 
istd_peak <- match_istd(object = temp_object, 
                        mode = 'hilic_pos', 
                        lib_version = 'v2', 
                        mz_ppm = 10, 
                        rt_second = 30)

order_injection <- temp_object@sample_info$sample_id %>% stringr::str_detect('PoolQC06') %>% which()

temp_plot <- lapply(seq_along(istd_peak$variable_id), function(i){
  x <- istd_peak$variable_id[i]
  temp_data <- temp_object %>% 
    intensity_plot(variable_id = x,
                   color_by = 'class',
                   order_by = 'injection.order', 
                   interactive = FALSE)
  
  temp_plot <- temp_data$data %>% 
    arrange(batch, injection.order) %>% 
    mutate(x = seq(n())) %>%
    ggplot(aes(x = x, y = int, color = class)) +
    geom_point() +
    scale_x_discrete(labels = NULL) +
    ZZWtool::ZZWTheme() +
    theme(axis.text.y.left = element_text(hjust = 0.5)) +
    ggtitle(paste(istd_peak$name_istd[i])) +
    geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed') +
    xlab('Injection Order')
})

ggsave(
  filename = "./00_raw_data_processing/01_input_data_cleaning/plots/serrf_normalized_ISTD_plots.pdf", 
  plot = gridExtra::marrangeGrob(temp_plot, nrow=3, ncol=3), 
  width = 15, height = 9
)


# RSD calculation 
temp_qc_rsd <- temp_object %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class == 'poolQC') %>% 
  extract_expression_data() %>%
  apply(., 1, function(x){sd(x)/mean(x)*100})

temp_all_rsd <- temp_object %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class != 'poolQC') %>%
  extract_expression_data() %>%
  apply(., 1, function(x){sd(x)/mean(x)*100})

temp_object <- temp_object %>%
  activate_mass_dataset(what = 'variable_info') %>%
  mutate(rsd_qc = temp_qc_rsd,
         rsd_all_sample = temp_all_rsd)

# all feature
temp_object %>%
  extract_variable_info() %>%
  dplyr::summarise(mean_qc_rsd = mean(rsd_qc, na.rm = TRUE),
                   mean_all_sample_rsd = mean(rsd_all_sample, na.rm = TRUE),
                   qc_rsd_less_30 = sum(rsd_qc <= 30)/n(),
                   all_sample_rsd_less_30 = sum(rsd_all_sample <= 30)/n())

# ISTD  
temp_object %>%
  extract_variable_info() %>%
  dplyr::filter(variable_id %in% istd_peak$variable_id) %>%
  dplyr::summarise(mean_qc_rsd = mean(rsd_qc, na.rm = TRUE),
                   mean_all_sample_rsd = mean(rsd_all_sample, na.rm = TRUE),
                   qc_rsd_less_30 = sum(rsd_qc <= 30)/n(),
                   all_sample_rsd_less_30 = sum(rsd_all_sample <= 30)/n())

rm(list = ls())




################################################################################
# HILIC negative ---------------------------------------------------------------
rm(list = ls());gc()

  # 1. Create tidymass object --------------------------------------------------
setwd('~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_neg/')

create_ibd_object(peak_table = 'Peak-table.csv',
                  merged_worklist = '~/Project/00_IBD_project/Data/00_IBD_worklists & sample_id/merged_worklist/merged_worklist_hilic_neg.xlsx',
                  dir_path = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_neg/00_raw_data_processing/',
                  mode = 'hilic_neg')



  # 2. data cleaning -----------------------------------------------------------
load('./00_raw_data_processing/01_input_data_cleaning/object.RData')

object_hilic_neg <- object
rm(object);gc()

# plot global PCA for all batches & samples before cleaning
temp_data_trans <- object_hilic_neg@expression_data %>% 
  rotate_df()
pca_pos <- pca(temp_data_trans, center = TRUE, scale = TRUE) 


plotIndiv(pca_pos,
          group = object_hilic_neg@sample_info$batch,
          ellipse = FALSE,
          legend = TRUE,
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = object_hilic_neg@sample_info$sample_id,
          title = 'Before cleaning')

  #   2.a remove low-quality features ------------------------------------------

# object: replace 0 value as NA
object_hilic_neg@expression_data <- object_hilic_neg@expression_data %>% 
  mutate(across(where(is.numeric), ~na_if(., 0)))

# criteria of removal:
# MV number >= 80% pool QC
# MV Non-IBD >= 50%
# IBD (B1->B1, B1->B2, B1->B3) >= 50%, either group

pool_qc_id <- object_hilic_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "poolQC") %>%
  pull(sample_id)

non_IDB_id <- object_hilic_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == "Non-IBD") %>%
  pull(sample_id)

B1_id <- object_hilic_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B1') %>%
  pull(sample_id)

B2_id <- object_hilic_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B2') %>%
  pull(sample_id)

B3_id <- object_hilic_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(phenotype_group == 'B1->B3') %>%
  pull(sample_id)

object_hilic_neg <- object_hilic_neg %>%
  mutate_variable_na_freq(according_to_samples = pool_qc_id) %>% 
  mutate_variable_na_freq(according_to_samples = non_IDB_id) %>% 
  mutate_variable_na_freq(according_to_samples = B1_id) %>% 
  mutate_variable_na_freq(according_to_samples = B2_id) %>% 
  mutate_variable_na_freq(according_to_samples = B3_id)

# extract_variable_info(object_hilic_neg) %>% head()
object_hilic_neg <- 
  object_hilic_neg %>% 
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2 & (na_freq.1 < 0.5 | na_freq.2 < 0.5 | na_freq.3 < 0.5 | na_freq.4 < 0.5))

dir.create('./00_raw_data_processing/01_input_data_cleaning', showWarnings = FALSE, recursive = TRUE)
save(object_hilic_neg, file = './00_raw_data_processing/01_input_data_cleaning/01_object_hilic_neg_noise_removal_240829.RData')

object_hilic_neg@expression_data %>% dim()

# # show missing value plot again
# show_missing_values(object = object_hilic_neg, show_column_names = FALSE, percentage = TRUE)


  #   2.b missing value imputation ---------------------------------------------
get_mv_number(object_hilic_neg)

object_hilic_neg <- 
  impute_mv(object = object_hilic_neg, method = "knn")

# check whether introduce zero values by KNN. If it exist, remove them
temp <- which(object_hilic_neg@expression_data < 0.0001, arr.ind = TRUE)
if (length(temp) > 0) {
  peak_removal <- which(object_hilic_neg@expression_data < 0.0001, arr.ind = TRUE) %>% rownames() %>% unique()
  temp_idx <- match(peak_removal, rownames(object_hilic_neg@expression_data))
  object_hilic_neg@expression_data[temp_idx,] %>% 
    rotate_df() %>% 
    rownames_to_column() %>% 
    tidyr::pivot_longer(cols = -rowname) %>% 
    dplyr::group_by(name) %>% 
    summarise(sum(value < 0.0001)) 
  
  object_hilic_neg <- object_hilic_neg %>% 
    activate_mass_dataset(what = 'variable_info') %>% 
    filter(!(variable_id %in% peak_removal))
}


object_hilic_neg@sample_info <- object_hilic_neg@sample_info %>% 
  mutate(batch = as.character(batch))

get_mv_number(object_hilic_neg)
show_missing_values(object = object_hilic_neg, show_column_names = FALSE, percentage = TRUE)

save(object_hilic_neg, 
     file = './00_raw_data_processing/01_input_data_cleaning/02_object_hilic_neg_imputation.RData')


  #   2.c normalization --------------------------------------------------------
# evaluate before normalization
load('./00_raw_data_processing/01_input_data_cleaning/02_object_hilic_neg_imputation.RData')
object_hilic_neg <- object_hilic_neg %>% 
  activate_mass_dataset(what = 'sample_info') %>% 
  filter(class != 'nistQC') %>% 
  mutate(class = case_when(class == 'poolQC' ~ 'QC',
                           class == 'sample' ~ 'Subject')) %>% 
  mutate(injection.order = seq(n()))

# original PCA plot
temp_object <- object_hilic_neg
temp_data_trans <- temp_object@expression_data %>% 
  rotate_df()

pca_pos <- mixOmics::pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$class
pca_batch <- temp_object@sample_info$batch
pca_var_id <- temp_object@sample_info$sample_id

dir.create('./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots', showWarnings = FALSE, recursive = TRUE)
pdf(file = './00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/original_PCA.pdf', 
    width = 8, height = 6, onefile = TRUE)

plotIndiv(pca_pos, 
          group = pca_group, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = pca_var_id,
          title = 'SVR normalization before normalization')

dev.off()

# show ISTD intensity plot 
istd_peak <- match_istd(object = object_hilic_neg, mode = 'hilic_neg', lib_version = 'v2', mz_ppm = 10, rt_second = 30)
order_injection <- object_hilic_neg@sample_info$sample_id %>% stringr::str_detect('PoolQC06') %>% which()

# temp_plot <- lapply(istd_peak$variable_id, function(x){
#   object_hilic_neg %>% 
#     intensity_plot(variable_id = x,
#                    color_by = 'class',
#                    order_by = 'injection.order', 
#                    interactive = FALSE) + 
#     scale_x_discrete(labels = NULL) +
#     ZZWtool::ZZWTheme() +
#     theme(axis.text.y.left = element_text(hjust = 0.5)) +
#     ggtitle(paste(x, 'Before intra-batch normalization')) +
#     geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed')
# })
# 
# ggsave(
#   filename = "./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/ISTD_plots_before_intrabatch_normalization.pdf", 
#   plot = gridExtra::marrangeGrob(temp_plot, nrow=1, ncol=1), 
#   width = 15, height = 9
# )

temp_plot <- lapply(seq_along(istd_peak$variable_id), function(i){
  x <- istd_peak$variable_id[i]
  temp_data <- temp_object %>% 
    intensity_plot(variable_id = x,
                   color_by = 'class',
                   order_by = 'injection.order', 
                   interactive = FALSE)
  
  temp_plot <- temp_data$data %>% 
    arrange(batch, injection.order) %>% 
    mutate(x = seq(n())) %>%
    ggplot(aes(x = x, y = int, color = class)) +
    geom_point() +
    scale_x_discrete(labels = NULL) +
    ZZWtool::ZZWTheme() +
    theme(axis.text.y.left = element_text(hjust = 0.5)) +
    ggtitle(paste(istd_peak$name_istd[i])) +
    geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed') +
    xlab('Injection Order')
})

ggsave(
  filename = "./00_raw_data_processing/01_input_data_cleaning/00_normalization_plots/ISTD_plots_before_intrabatch_normalization_2.pdf", 
  plot = gridExtra::marrangeGrob(temp_plot, nrow=3, ncol=3), 
  width = 15, height = 9
)


# normalization via the SERRF 

load('./00_raw_data_processing/01_input_data_cleaning/02_object_hilic_neg_imputation.RData')

library(tidyverse)
library(tidymass)
library(sjmisc)
library(DoddLabQC)
library(DoddLabTool)
library(DoddLabSERRF)


# convert data for serrf
table_serrf <- convert_object_for_serrf(object = object_hilic_neg)

dir.create('./00_raw_data_processing/01_input_data_cleaning/01_serrf', 
           showWarnings = FALSE, recursive = TRUE)

saveRDS(table_serrf, 
        file = './00_raw_data_processing/01_input_data_cleaning/01_serrf/raw_data_240829.RDS')

# run serrf normalization
normalize_serrf(data_file = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_neg/00_raw_data_processing/01_input_data_cleaning/01_serrf/raw_data_240829.RDS',
                dir_path = '~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_neg/00_raw_data_processing/01_input_data_cleaning/01_serrf/')

# 
# merge SERRF result and update the tidymass object
setwd('~/Project/00_IBD_project/Data/20240828_IBD_B001_B042_analysis/hilic_neg/')

library(tidyverse)
library(tidymass)
library(sjmisc)
library(DoddLabQC)
library(DoddLabTool)
library(DoddLabSERRF)

load('./00_raw_data_processing/01_input_data_cleaning/02_object_hilic_neg_imputation.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/serrf_result/serrf_result.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/p.RData')
load('./00_raw_data_processing/01_input_data_cleaning/01_serrf/f.RData')

object_hilic_neg_serrf <- merge_serrf_result(object = object_hilic_neg, 
                                             serrf_result = serrf_result, 
                                             p = p,
                                             f = f)

save(object_hilic_neg_serrf,
     file = './00_raw_data_processing/01_input_data_cleaning/03_object_hilic_neg_serrf.RData')

rm(list = c('object_hilic_neg', 'serrf_result', 'p', 'f'));gc()


# evaluate the SERRF normalization result 
temp_object <- object_hilic_neg_serrf %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class %in% c('sample', 'poolQC'))


temp_data_trans <- temp_object@expression_data %>% 
  # activate_mass_dataset(what = 'expression_data') %>%
  # mutate(across(everything(), as.numeric)) %>%
  rotate_df()

pca_pos <- mixOmics::pca(temp_data_trans, center = TRUE, scale = TRUE) 
pca_group <- temp_object@sample_info$class
pca_group2 <- temp_object@sample_info$group
pca_batch <- temp_object@sample_info$batch
pca_var_id <- temp_object@sample_info$sample_id

dir.create('./00_raw_data_processing/01_input_data_cleaning/plots', showWarnings = FALSE, recursive = TRUE)
pdf(file = './00_raw_data_processing/01_input_data_cleaning/plots/serrf_normalized_PCA.pdf', 
    width = 8, height = 6, onefile = TRUE)

plotIndiv(pca_pos, 
          group = pca_group, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_batch, 
          ellipse = FALSE, 
          legend = TRUE, 
          legend.title = 'Batch number',
          # pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

plotIndiv(pca_pos, 
          group = pca_group2, 
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Class',
          pch = 16,
          ind.names = pca_var_id,
          title = 'SERRF normalization')

dev.off()



# show ISTD intensity plot 
istd_peak <- match_istd(object = temp_object, 
                        mode = 'hilic_neg', 
                        lib_version = 'v2', 
                        mz_ppm = 10, 
                        rt_second = 30)

order_injection <- temp_object@sample_info$sample_id %>% stringr::str_detect('PoolQC06') %>% which()

temp_plot <- lapply(seq_along(istd_peak$variable_id), function(i){
  x <- istd_peak$variable_id[i]
  temp_data <- temp_object %>% 
    intensity_plot(variable_id = x,
                   color_by = 'class',
                   order_by = 'injection.order', 
                   interactive = FALSE)
  
  temp_plot <- temp_data$data %>% 
    arrange(batch, injection.order) %>% 
    mutate(x = seq(n())) %>%
    ggplot(aes(x = x, y = int, color = class)) +
    geom_point() +
    scale_x_discrete(labels = NULL) +
    ZZWtool::ZZWTheme() +
    theme(axis.text.y.left = element_text(hjust = 0.5)) +
    ggtitle(paste(istd_peak$name_istd[i])) +
    geom_vline(xintercept = order_injection, colour = 'red', linetype = 'dashed') +
    xlab('Injection Order')
})

ggsave(
  filename = "./00_raw_data_processing/01_input_data_cleaning/plots/serrf_normalized_ISTD_plots.pdf", 
  plot = gridExtra::marrangeGrob(temp_plot, nrow=3, ncol=3), 
  width = 15, height = 9
)


# RSD calculation 
temp_qc_rsd <- temp_object %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class == 'poolQC') %>% 
  extract_expression_data() %>%
  apply(., 1, function(x){sd(x)/mean(x)*100})

temp_all_rsd <- temp_object %>% 
  activate_mass_dataset(what = 'sample_info') %>%
  filter(class != 'poolQC') %>%
  extract_expression_data() %>%
  apply(., 1, function(x){sd(x)/mean(x)*100})

temp_object <- temp_object %>%
  activate_mass_dataset(what = 'variable_info') %>%
  mutate(rsd_qc = temp_qc_rsd,
         rsd_all_sample = temp_all_rsd)

# all feature
temp_object %>%
  extract_variable_info() %>%
  dplyr::summarise(mean_qc_rsd = mean(rsd_qc, na.rm = TRUE),
                   mean_all_sample_rsd = mean(rsd_all_sample, na.rm = TRUE),
                   qc_rsd_less_30 = sum(rsd_qc <= 30)/n(),
                   all_sample_rsd_less_30 = sum(rsd_all_sample <= 30)/n())

# ISTD  
temp_object %>%
  extract_variable_info() %>%
  dplyr::filter(variable_id %in% istd_peak$variable_id) %>%
  dplyr::summarise(mean_qc_rsd = mean(rsd_qc, na.rm = TRUE),
                   mean_all_sample_rsd = mean(rsd_all_sample, na.rm = TRUE),
                   qc_rsd_less_30 = sum(rsd_qc <= 30)/n(),
                   all_sample_rsd_less_30 = sum(rsd_all_sample <= 30)/n())

rm(list = ls())

