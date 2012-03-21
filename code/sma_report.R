
##********************************************************##
##    LPA report                                          ##
##********************************************************##

gene_list_files <- c(
  "../data/sma/known_genes.txt",
  "../data/sma/sma_disease_recessive_selected_genes_maf01.txt",
  "../data/sma/sma_activator_dominant_selected_genes_maf01.txt"
)
group_names <- c("Known_genes","DISEASE_RECESSIVE","ACTIVATOR_DOMINANT")

# Output data

  # algorithm params
numinterm <- 1
zero_intermediates <- TRUE

  # output 
network_outdir <- "../executions/sma"


if(!file.exists(network_outdir)){
  dir.create(network_outdir)
}

#source("interactome_analysis.R")