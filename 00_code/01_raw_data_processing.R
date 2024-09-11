################################################################################
# Raw data processing ----------------------------------------------------------
# The script run in Stanford SCG cluster.

# Package version:
# DoddLabRawMS_0.1.16: https://github.com/DoddLab/DoddLabRawMS
# xcms: 3.20.0

# install required packages
# install.packages("devtools")
# devtools::install_github("DoddLab/DoddLabRawMS")
# BiocManager::install("xcms")

################################################################################
# C18 positive -----------------------------------------------------------------

# load required packages
library(dplyr)
library(tidyr)
library(xcms)
library(DoddLabRawMS)

# set working directory
# I prefer use absolute address. Note: change the folder
setwd("/labs/ddodd2/ZhiweiZhou/20231201_IBD_new_batch/c18_pos")

path <- '/labs/ddodd2/ZhiweiZhou/20231201_IBD_new_batch/c18_pos'

# raw data procssing
parameter_set <- initialize_raw_parameter_class(column = 'c18')
parameter_set@para_peak_detection$nSlaves <- 20
process_raw_data(parameter_set = parameter_set,
                 path = path)


################################################################################
# C18 negative -----------------------------------------------------------------

# load required packages
library(dplyr)
library(tidyr)
library(xcms)
library(DoddLabRawMS)

# set working directory
# I prefer use absolute address. Note: change the folder
setwd("/labs/ddodd2/ZhiweiZhou/20231201_IBD_new_batch/c18_neg")

path <- '/labs/ddodd2/ZhiweiZhou/20231201_IBD_new_batch/c18_neg'

# raw data procssing
parameter_set <- initialize_raw_parameter_class(column = 'c18')
parameter_set@para_peak_detection$nSlaves <- 20
process_raw_data(parameter_set = parameter_set,
                 path = path)


################################################################################
# HILIC positive ---------------------------------------------------------------

# load required packages
library(dplyr)
library(tidyr)
library(xcms)
library(DoddLabRawMS)

# set working directory
# I prefer use absolute address. Note: change the folder
setwd("/labs/ddodd2/ZhiweiZhou/20231201_IBD_new_batch/hilic_pos")

path <- '/labs/ddodd2/ZhiweiZhou/20231201_IBD_new_batch/hilic_pos'

# raw data procssing
parameter_set <- initialize_raw_parameter_class(column = 'hilic')
parameter_set@para_peak_detection$nSlaves <- 20
parameter_set@para_peak_detection$ppm <- 30
process_raw_data(parameter_set = parameter_set,
                 path = path)

################################################################################
# HILIC negative ---------------------------------------------------------------

# load required packages
library(dplyr)
library(tidyr)
library(xcms)
library(DoddLabRawMS)

# set working directory
# I prefer use absolute address. Note: change the folder
setwd("/labs/ddodd2/ZhiweiZhou/20231201_IBD_new_batch/hilic_neg")

path <- '/labs/ddodd2/ZhiweiZhou/20231201_IBD_new_batch/hilic_neg'

# raw data procssing
parameter_set <- initialize_raw_parameter_class(column = 'hilic')
parameter_set@para_peak_detection$nSlaves <- 20
parameter_set@para_peak_detection$ppm <- 30
process_raw_data(parameter_set = parameter_set,
                 path = path)

