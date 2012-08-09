
home_path <- getwd()

source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("postprocess.r")

#
# INPUT DATA
#

gene_list_files <- c(
  "../data/lpa_results_09_03_2012_interactome/lpa1_selected_genes.txt",
  "../data/lpa_results_09_03_2012_interactome/lpa2_selected_genes.txt",
  "../data/lpa_results_09_03_2012_interactome/lpa3_selected_genes.txt",
  "../data/lpa_results_09_03_2012_interactome/lpa4_selected_genes.txt",
  "../data/lpa_results_09_03_2012_interactome/lpa5_selected_genes.txt"
)

group_names <- c("LPA1","LPA2","LPA3","LPA4","LPA5")


#
# PREPARE DATA
#

raw_gene_list2 <- load_family_genes(gene_list_files=gene_list_files,group_names=group_names,hasHeader=T)

interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))

global_score_methods <- c("no_zero_mean","no_zero_max")

#
# COMPUTE SCORE
#

multi_score <- compute_multi_score(raw_gene_list,interactomes,default_score,global_score_methods=global_score_methods,verbose=T)