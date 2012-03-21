
##********************************************************##
##    LPA report                                          ##
##********************************************************##

# Input data
# gene_list_files <- c(
#   "../data/lpa_results_09_11_2011/LPA1_gene_list.txt",
#   "../data/lpa_results_09_11_2011/LPA2_gene_list.txt",
#   "../data/lpa_results_09_11_2011/LPA3_gene_list.txt",
#   "../data/lpa_results_09_11_2011/LPA4_gene_list.txt",
#   "../data/lpa_results_09_11_2011/LPA5_gene_list.txt"
#   )  
gene_list_files <- c(
  "../data/lpa_results_01_03_2012/lpa1_genes_01_03_2012.txt",
  "../data/lpa_results_01_03_2012/lpa2_genes_01_03_2012.txt",
  "../data/lpa_results_01_03_2012/lpa3_genes_01_03_2012.txt",
  "../data/lpa_results_01_03_2012/lpa4_genes_01_03_2012.txt",
  "../data/lpa_results_01_03_2012/lpa5_genes_01_03_2012.txt"
  )
group_names <- c("LPA1","LPA2","LPA3","LPA4","LPA5")

# Output data

  # algorithm params
numinterm <- 1
zero_intermediates <- TRUE

  # output 
network_outdir <- "../executions/lpa_results"


#source("interactome_analysis.R")