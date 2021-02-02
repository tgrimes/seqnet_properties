# These scripts are based on version 1.1.0 of the SeqNet package.
library(devtools)
install_version("SeqNet", version = "1.1.0", repos = "http://cran.us.r-project.org")
library(SeqNet)

# Additional packages that are used in these scripts:
library(tidyverse)
library(igraph)
library(viridis)

# Create a directory to store the output, if none exists.
if(!dir.exists("output")) dir.create("output")

# Analyze the properties of the tissue-specific networks.
source(file.path("src", "summarize_humanbase.R"))
init_human_base()

# Analyze the properties of the organism networks.
dir_inetbio_summary <- file.path("output", "inetbio_summary")
file_list <- dir(dir_inetbio_summary)
df_inetbio <- NULL
for(file in file_list) {
  df_inetbio <- rbind(df_inetbio,
                      readRDS(file.path(dir_inetbio_summary, file)))
}
# Remove the tissue-specfic drosophila networks.
df_inetbio <- df_inetbio[!grepl("TS", df_inetbio$tissue), ]

# Simulate networks from the ER, WS, BA, and SeqNet generators.
n_sim <- 50 # Number of networks to generate for each simulation setting.
source(file.path("src", "simulation_from_ER.R"))
source(file.path("src", "simulation_from_WS.R"))
source(file.path("src", "simulation_from_BA.R"))
source(file.path("src", "simulation_from_SeqNet.R"))

# Create figures to summarize the statistical properties of each generators
# with respect to each topological measure. These are provided in the
# supplementary materials.
source(file.path("src", "create_supplementary_figures"))
