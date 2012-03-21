
simulate_family_genes <- function(nfams,numgenes,distance_matrix=NULL,interactome=NULL,fam_name_prefix="fam_") {

  # family names
  fam_names <- paste(fam_name_prefix,c(1:nfams),sep="")
  
  # select a disease gene per family
  if(is.null(distance_matrix)){
#     cat("computing distance matrix...\n")
    distance_matrix <- compute_distance_matrix(interactome)
  } else {
#     cat("distance matrix provided...\n")
  }
  
  correct <- F
  while(correct==F){
    row <- as.integer(runif(1,min=1,max=nrow(distance_matrix)))
    pivot <- colnames(distance_matrix)[row]
    neighborhood <- names(which(distance_matrix[row,]<3 & distance_matrix[row,]>0))
    if(length(neighborhood)>nfams){
      fam_disease_genes <- c(pivot,sample(neighborhood,nfams-1))
      names(fam_disease_genes) <- fam_names
      correct <- T
    } 
  }
    
  # add other random genes
  interactors <- rownames(distance_matrix)
  raw_gene_list <- list()
  random_gene_list <- list()
  for(i in 1:nfams){
    random_list <- sample(interactors,numgenes)
    raw_gene_list[[fam_names[i]]] <- c(as.character(fam_disease_genes[i]),random_list)  
    random_gene_list[[fam_names[i]]] <- random_list 
  }

  return(list(
    raw_gene_list=raw_gene_list,
    random_gene_list=random_gene_list,
    pivot=pivot,
    fam_disease_genes=fam_disease_genes
  ))
  
}


evaluate_prioritization <- function(scores_frame,fam_simulation,paint=T){
  
  fam_genes_scores <- scores_frame$score_table[fam_simulation$fam_disease_genes, "score" ]
  names(fam_genes_scores) <- fam_simulation$fam_disease_genes
  
  fam_genes_positions <- sapply(fam_simulation$fam_disease_genes, function(x) which(scores_frame$score_table$gene==x))
  fam_genes_relative_positions <- 1-(fam_genes_positions/nrow(scores_frame$score_table))
  mean_relative_position <- mean(fam_genes_relative_positions)
  
  fam_genes_positions_with_inter <- sapply(fam_simulation$fam_disease_genes,function(x) which(scores_frame$score_table_with_intermediates$gene==x))
  fam_genes_relative_positions_with_inter <- 1-(fam_genes_positions_with_inter/nrow(scores_frame$score_table_with_intermediates))
  mean_relative_position_with_inter <- mean(fam_genes_relative_positions_with_inter)
    
  if(paint){
    par(mfrow=c(2,1))
    paint_evaluation(scores_frame$score_table$score,fam_genes_scores)
    paint_evaluation(scores_frame$score_table_with_intermediates$score,fam_genes_scores,main="Score distribution (with intermediates)")
  }
  
  return(list(
    fam_genes_scores=fam_genes_scores,
    fam_genes_positions=fam_genes_positions,
    fam_genes_relative_positions=fam_genes_relative_positions,
    fam_genes_positions_with_inter=fam_genes_positions_with_inter,
    fam_genes_relative_positions_with_inter=fam_genes_relative_positions_with_inter,
    mean_relative_position=mean_relative_position,
    mean_relative_position_with_inter=mean_relative_position_with_inter
  ))
  
}

paint_evaluation <- function(all_scores,fam_genes_scores,main="Score distribution",bins=20){
  
  h <- hist(all_scores,bins,probability=T,main=main,xlab="score")
  lines(density(all_scores),col="blue")
  for(i in fam_genes_scores){
    lines(c(i,i),c(0,max(h$density*0.9)),col="red",lwd=2)  
  }
  
}