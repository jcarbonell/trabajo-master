

default_score <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
  
  scores_by_distance <- numeric(length(weights))
  
  distance <- 0
  distance_prior <- 1  
  
  for(j in 1:length(weights)){
                  
      if(distance==0){
        selected_distance_nodes <- node
        distance_prior <- 1
      } else {
        node_distances <- distance_matrix[node,]
        selected_distance_nodes <- names(node_distances[which(node_distances==distance)])
        distance_prior <- distance_priors[distance]
      }
        
      selected_distance_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
      selected_distance_fams <- setdiff(selected_distance_fams,inter_family)
            
      if(is.null(selected_distance_fams)){
          scores_by_distance[j] <- 0
      } else {                
          scores_by_distance[j] <- (length(selected_distance_fams)*weights[j])#/distance_prior                   
      }
  
      distance <- distance+1
      
  }

  return(list(
    scores_by_distance=scores_by_distance,
    score=sum(scores_by_distance)
  ))
  
}

default_score2 <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
  
  scores_by_distance <- numeric(length(weights))
  
  distance <- 0
  distance_prior <- 1  
  
  for(j in 1:length(weights)){
    
    if(distance==0){
      selected_distance_nodes <- node
      distance_prior <- 1
    } else {
      node_distances <- distance_matrix[node,]
      selected_distance_nodes <- names(node_distances[which(node_distances==distance)])
      distance_prior <- distance_priors[distance]
    }
    
    node_fam <- unique(unlist(super_list_fams[node]))
    
    selected_distance_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
    selected_distance_fams <- setdiff(selected_distance_fams,inter_family)
    selected_distance_fams <- setdiff(selected_distance_fams,node_fam)
        
    if(is.null(selected_distance_fams)){
      scores_by_distance[j] <- 0
    } else {                
      scores_by_distance[j] <- (length(selected_distance_fams)*weights[j])#/distance_prior                   
    }
    
    distance <- distance+1
    
  }
  
  return(list(
    scores_by_distance=scores_by_distance,
    score=sum(scores_by_distance)
    ))
  
}

default_score3 <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
  
  scores_by_distance <- numeric(length(weights))
  
  distance <- 0
  distance_prior <- 1  
  
  for(j in 1:length(weights)){
    
    if(distance==0){
      selected_distance_nodes <- node
      distance_prior <- 1
    } else {
      node_distances <- distance_matrix[node,]
      selected_distance_nodes <- names(node_distances[which(node_distances==distance)])
      distance_prior <- distance_priors[distance]
    }
    
    node_fam <- unique(unlist(super_list_fams[node]))
    
    selected_distance_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
    selected_distance_fams <- setdiff(selected_distance_fams,inter_family)
    selected_distance_fams <- setdiff(selected_distance_fams,node_fam)
        
    if(is.null(selected_distance_fams)){
      scores_by_distance[j] <- 0
    } else {                
      scores_by_distance[j] <- (length(selected_distance_fams)*weights[j])/distance_prior                   
    }
    
    distance <- distance+1
    
  }
  
  return(list(
    scores_by_distance=scores_by_distance,
    score=sum(scores_by_distance)
    ))
  
}


compute_score <- function(raw_gene_list,interactome,score_function,radio=5,numinter=1,inter_family="INTER"){
  
  # interactome pruning
  all_raw_genes <- unique(unlist(raw_gene_list))
  subnet <- get.all.shortest.paths.Josete(interactome,all_raw_genes,5)
  subnet_summary <- describe_interactome(subnet,plot=F)
  distance_priors <- subnet_summary["Median",]/sum(subnet_summary["Median",])
#   cat("     (",nrow(subnet),"interactions)\n")
    
  # init interactors
  interactors <- unique(c(subnet[,1],subnet[,3]))
  ninteractors <- length(interactors)
  
  # init super list
  super_list <- unique(c(all_raw_genes,interactors))
  nsuper <- length(super_list)
  
  # compute fam hits
  family_names <- names(raw_gene_list)
  nfams <- length(unique(family_names))
  genes_by_family <- sapply(raw_gene_list,length)
    # obtain fams by gene
  super_list_fams <- list()
  get_gene_fams <- function(gene){
    gene_fams <- names(which(lapply(raw_gene_list,function(fam) is.element(gene,fam))==T))
    if(length(gene_fams)==0){
      gene_fams <- inter_family
    }
    return(gene_fams)
  }
  super_list_fams <- lapply(super_list,get_gene_fams)
  names(super_list_fams) <- super_list

  # compute node distances
  distance_matrix <- compute_distance_matrix(subnet)
    
  # compute weights
  steps <- radio+1
  distances <- 0:radio
  raw_weights <- 1-pnorm(radio:(radio*2),mean=radio,sd=2)
  weights <- raw_weights/sum(raw_weights)

  # init working variables
  node_scores <- numeric(length(super_list))
  names(node_scores) <- super_list
  scores_by_distance <- matrix(rep(0,nsuper*steps),ncol=steps)
  rownames(scores_by_distance) <- super_list
  colnames(scores_by_distance) <- distances
  uninteractors <- c()
  
  # compute score by genes
  for(w in 1:nsuper){
    
    node <- super_list[w]
   
    i <- which(interactors==node)

    # is it an interactor??
    if(length(i)==0){
        
      if(is.null(super_list_fams[[node]])){
        
        node_scores[w] <- 0
        
      } else {
                
        uninteractors <- c(uninteractors,node)
        
        node_score <- score_function(node,distance_matrix,weights[1],distance_priors,super_list_fams,inter_family=inter_family)
        node_scores[w] <- node_score$score
        scores_by_distance[w,1] <- node_score$scores_by_distance
        
      }
          
    } else {
       
      node_score <- score_function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family=inter_family)
      node_scores[w] <- node_score$score
      scores_by_distance[w,] <- node_score$scores_by_distance
          
    }
    
  }


  # normalization
  node_scores <- node_scores/max(node_scores)
  
  # fam_hits  
  get_fam_hit <- function(fams){return(paste(fams,collapse=","))}
  fam_hits <- sapply(super_list_fams,get_fam_hit)
  
  # sorting index
  sorting_index <- rev(order(node_scores))

  # result table
  result_table <- data.frame(
    gene=super_list[sorting_index],
    score=node_scores[sorting_index],  
    fams=fam_hits[sorting_index]
  )
  
  # remove intermediates
  no_inter <- which(result_table$fams!=inter_family)
  
  # return
  scores <- list(
#     sorting_index=sorting_index,
#     super_list=super_list,
#     node_scores=node_scores,
#     fam_hits=fam_hits,
    distance_priors=distance_priors,
    
    uninteractor_ratio=length(uninteractors)/length(all_raw_genes),
    uninteractors=uninteractors,
    
    score_table=result_table[no_inter,],
    score_table_with_intermediates=result_table,
    
    scores_by_distance=scores_by_distance,
    
    fam_hits=fam_hits
    
  )
  
  return(scores)

}



