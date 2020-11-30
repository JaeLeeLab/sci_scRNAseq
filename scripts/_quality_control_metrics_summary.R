
####### Count-based quality control - Final metrics summary #######


# This script is to summarize all quality control metric values for cells that
# passed QC. I also include supplementary read-level information.


# Data Import -------------------------------------------------------------


# Load libraries and set directories
require('ggplot2')
require('Matrix')
require('dplyr')
data_in <- './data/QC_filtered_feature_bc_matrix/'
# data_out <- './data/QC_filtered_feature_bc_matrix/'
results_out <- './results/quality_control_final_summary/'
# ref_out <- './ref/'
# dir.create(path = data_out)
dir.create(path = results_out)


