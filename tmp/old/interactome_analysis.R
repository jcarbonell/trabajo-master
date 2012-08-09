
source("igraph_utils.R")
source("nh_simulation_lib.R")
source("shortestpath-mariangela.R")
source("nh_lib.R")
source("score.r")
source("utils.r")


home_path <- getwd()


# FAMILY GENES LOADING

raw_gene_list <- load_family_genes(group_names,gene_list_files)
genes_by_family <- lapply(raw_gene_list,length)
all_raw_genes <- unique(unlist(raw_gene_list))


# INTERACTOMES LOADING

interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
ninteractomes <- length(interactomes)

# SUB-NETWORK INITIALIZATION

#subnets <- get_subnets(all_raw_genes,interactomes,5)


# distance_matrices <- list()
# for(i in 1:length(interactomes)){
#   cat("Computing distance matrix" ,names(interactomes)[i],"...\n")
#   distance_matrices[[ names(interactomes)[i] ]] <- compute_distance_matrix(interactomes[[i]])
#   cat("      cleaning memory...\n")
#   gc()
# }


# SCORE COMPUTATION


scores_frames <- list()
for(i in 1:length(interactomes)){
  cat("Computing score with",names(interactomes)[i],"interactome\n")
  scores_frames[[ names(interactomes)[i] ]] <- compute_score(raw_gene_list,interactomes[[i]],default_score2)
}


# SUMMARY

scored_genes <- unique(unlist(lapply(scores_frames,function(score_frame) return(rownames(score_frame$score_table)))))

score_matrix <- matrix(rep(0,length(scored_genes)*ninteractomes),ncol=ninteractomes)
rownames(score_matrix) <- scored_genes
colnames(score_matrix) <- names(interactomes)
for(i in 1:ninteractomes){
  score_matrix[,i] <- scores_frames[[i]]$score_table[scored_genes,"score"]
}
global_score <- rowSums(score_matrix)
score_rank <- order(global_score,decreasing=T)
score_matrix <- cbind(global_score,score_matrix)
score_matrix <- score_matrix[score_rank,]

score_matrix_frame <- as.data.frame(score_matrix)
score_matrix_frame$fams <- scores_frames[[i]]$fam_hits[rownames(score_matrix)]

score_matrix_file <- paste(network_outdir,"/score_matrix.txt",sep="")
write(paste("gene\t",paste(colnames(score_matrix_frame),collapse="\t"),sep=""),file=score_matrix_file)
write.table(score_matrix_frame,file=score_matrix_file,sep="\t",quote=F,col.names=F,append=T)




