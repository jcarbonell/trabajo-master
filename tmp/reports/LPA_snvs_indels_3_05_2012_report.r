
home_path <- getwd()

source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("postprocess.r")

#
# INPUT DATA
#

gene_list_files <- c(
  "../data/LPA_snvs_indels_3_05_2012/lpa1_selected_gene_list.txt",
  "../data/LPA_snvs_indels_3_05_2012/lpa2_selected_gene_list.txt",
  "../data/LPA_snvs_indels_3_05_2012/lpa3_selected_gene_list.txt",
  "../data/LPA_snvs_indels_3_05_2012/lpa4_selected_gene_list.txt",
  "../data/LPA_snvs_indels_3_05_2012/lpa5_selected_gene_list.txt"
  )

group_names <- c("LPA1","LPA2","LPA3","LPA4","LPA5")


#
# PREPARE DATA
#

raw_gene_list <- load_family_genes(gene_list_files=gene_list_files,group_names=group_names)

interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))

global_score_methods <- c("no_zero_mean","no_zero_max")

#
# COMPUTE SCORE
#

multi_score <- compute_multi_score(raw_gene_list,interactomes,default_score,global_score_methods=global_score_methods,verbose=T)