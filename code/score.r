

default_score <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
  
  scores_by_distance <- numeric(length(weights))
  
  distance <- 0
  distance_prior <- 1  
  
#   for(j in 1:length(weights)){
  
  per_distance <- function(j){
    
      distance = j-1
      
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
          scores_by_distance <- 0
      } else {                
          scores_by_distance <- (length(selected_distance_fams)*weights[j])#/distance_prior                   
      }
      
      return(scores_by_distance)
        
  }

  scores_by_distance <- sapply(1:length(weights),per_distance)
  
  return(list(
    scores_by_distance=scores_by_distance,
    score=sum(scores_by_distance)
  ))
  
}

default_score1b <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
  
  scores_by_distance <- numeric(length(weights))
  
  distance <- 0
  distance_prior <- 1  
  
#   for(j in 1:length(weights)){
  
  per_distance <- function(j){
    
      distance = j-1
                  
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
      
      if(distance>0){
        node_fam <- unique(unlist(super_list_fams[node]))
        selected_distance_fams <- setdiff(selected_distance_fams,node_fam)
      }
  
      if(is.null(selected_distance_fams)){
          scores_by_distance <- 0
      } else {                
          scores_by_distance <- (length(selected_distance_fams)*weights[j])#/distance_prior                   
      }
  
      return(scores_by_distance)

  }

  scores_by_distance <- sapply(1:length(weights),per_distance)
  
  score <- sum(scores_by_distance) - weights[1]
        
  return(list(
    scores_by_distance=scores_by_distance,
    score=score
    ))
  
}

default_score2 <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
  
  scores_by_distance <- numeric(length(weights))
  
  distance <- 0
  distance_prior <- 1  
  
#   for(j in 1:length(weights)){
  
  per_distance <- function(j){
    
      distance = j-1
    
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
      scores_by_distance <- 0
    } else {                
      scores_by_distance <- (length(selected_distance_fams)*weights[j])#/distance_prior                   
    }
    
    return(scores_by_distance)
    
  }
  
  scores_by_distance <- sapply(1:length(weights),per_distance)
  
  return(list(
    scores_by_distance=scores_by_distance,
    score=sum(scores_by_distance)
    ))
  
}

default_score2b <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
  
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
  
  score <- sum(scores_by_distance) - 1*weights[j]
        
  return(list(
    scores_by_distance=scores_by_distance,
    score=score
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

default_score4 <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
  
  scores_by_distance <- numeric(length(weights))
  
  distance <- 0
  distance_prior <- 1  
  
  for(j in 1:length(weights)){
    
    if(distance==0){
      selected_distance_nodes <- node
      distance_prior <- distance_priors[1]
    } else {
      node_distances <- distance_matrix[node,]
      selected_distance_nodes <- names(node_distances[which(node_distances==distance)])
      distance_prior <- distance_priors[distance]
    }
    
    node_fam <- unique(unlist(super_list_fams[node]))
    
    selected_distance_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
    selected_distance_fams <- setdiff(selected_distance_fams,inter_family)
    if(distance>0){
      selected_distance_fams <- setdiff(selected_distance_fams,node_fam)
    }
    
    if(is.null(selected_distance_fams)){
      scores_by_distance[j] <- 0
    } else {                
      scores_by_distance[j] <- (length(selected_distance_fams)*weights[j])/(distance_prior)
    }
    
    
    distance <- distance+1
    
  }
  
  score <- sum(scores_by_distance) - 1*weights[j]
        
  return(list(
    scores_by_distance=scores_by_distance,
    score=score
    ))
  
}

compute_score <- function(raw_gene_list,interactome,score_function,radio=5,numinter=1,inter_family="INTER",verbose=T){
  
  # interactome pruning
  all_raw_genes <- unique(unlist(raw_gene_list))
  subnet <- get.all.shortest.paths.Josete(interactome,all_raw_genes,5)
  subnet_summary <- describe_interactome(subnet,plot=F)
  distance_priors <- subnet_summary["Median",]/sum(subnet_summary["Median",])
    
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
#   for(w in 1:nsuper){    
#     node <- super_list[w]
   
  per_node <- function(node){  
    
    i <- which(interactors==node)

    # is it an interactor??
    if(length(i)==0){
        
      if(is.null(super_list_fams[[node]])){
        
        node_score <- 0
        scores_by_distance <- rep(0,length(weights))
        
      } else {
                   
        comp <- score_function(node,distance_matrix,weights[1],distance_priors,super_list_fams,inter_family=inter_family)
        node_score <- comp$score
        scores_by_distance <- c(comp$scores_by_distance,rep(0,length(weights)-1))
        
      }
      
      neighborhood <- c()
      neighbour_distance <- c()
          
    } else {
       
      comp <- score_function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family=inter_family)
      node_score <- comp$score
      scores_by_distance <- comp$scores_by_distance
          
      neighborhood <- colnames(distance_matrix)[which(distance_matrix[node,]>0 & distance_matrix[node,]<=radio)]
      neighbour_distance <- distance_matrix[node,neighborhood]
      
    }
    
    return(list(
      node_score=node_score,
      scores_by_distance=scores_by_distance,
      neighborhood=neighborhood,
      neighbour_distance=neighbour_distance
    ))
    
  }
  
#   }

  # process
  processed <- lapply(super_list,per_node)
  
  # recover score
  node_scores <- sapply(processed,function(x)x$node_score)
  names(node_scores) <- super_list
  scores_by_distance <- do.call("rbind",lapply(processed,function(x)x$scores_by_distance))
  rownames(scores_by_distance) <- super_list
  colnames(scores_by_distance) <- distances  
  neighborhoods <- lapply(processed,function(x)x$neighborhood)
  names(neighborhoods) <- super_list
  neighbour_distance <- lapply(processed,function(x)x$neighbour_distance)
  names(neighbour_distance) <- super_list
  
  # group score
  gene_group_score <- function(gene){
    neighborhood <- neighborhoods[[gene]]
    neighbour_distance <- neighbour_distance[[gene]]
    neighborhood_scores <- node_scores[neighborhood]
    neighborhood_weights <- weights[neighbour_distance]
    group_score <- median(neighborhood_scores*neighborhood_weights)
    group_score
  }
  group_scores <- sapply(super_list,gene_group_score)
  
  # uninteractors
  uninteractors <- super_list[]
 
  # normalization
  node_scores <- (node_scores-min(node_scores))/(max(node_scores)-min(node_scores))
  group_scores <- (group_scores-min(group_scores))/(max(group_scores)-min(group_scores)) 
  
  # fam_hits  
  get_fam_hit <- function(fams){return(paste(fams,collapse=","))}
  fam_hits <- sapply(super_list_fams,get_fam_hit)
  
  # sorting index
  sorting_index <- rev(order(node_scores))

  # result table
  result_table <- data.frame(
    gene=super_list[sorting_index],
    score=node_scores[sorting_index],
    group_score=group_scores[sorting_index],
    fams=fam_hits[sorting_index],
    stringsAsFactors=F
  )
  
  # remove intermediates
  no_inter <- which(result_table$fams!=inter_family)
  
  # return
  scores <- list(
    # score
    score_table=result_table[no_inter,],
    score_table_with_intermediates=result_table,    
    scores_by_distance=scores_by_distance,    
    fam_hits=fam_hits,
    # neighborhood
    neighborhoods=neighborhoods,
    neighbour_distance=neighbour_distance,
    # weights
    distance_priors=distance_priors,
    weights=weights,
    # uninteractors
    uninteractors=uninteractors,
    uninteractor_ratio=length(uninteractors)/length(all_raw_genes)
  )
  
  return(scores)

}


compute_multi_score <- function(raw_gene_list,interactomes,score_function,radio=5,numinter=1,inter_family="INTER",verbose=T,global_score_methods=NULL){
  
  # run and evaluate
  scores_frames <- list()
  
  for(i in 1:length(interactomes)){
    
    if(verbose==T){
      cat("Computing score with",names(interactomes)[i],"interactome\n")
    }
    scores_frames[[ names(interactomes)[i] ]] <- compute_score(raw_gene_list,interactomes[[i]],score_function,verbose=verbose,radio=radio)
    
  }

  global_score_table <- compute_global_score_table(scores_frames,inter=F,method=global_score_methods)
  
  ###################################
#   cat("preparing group score\n")
#   
#   get_neighborhood <- function(score_frame,gene){
#     score_frame$neighborhoods[gene]
#   }  
#   get_global_neighborhood <- function(gene){
#     unique(unlist(lapply(scores_frames,get_neighborhood,gene=gene)))
#   }
#   all_raw_genes <- rownames(global_score_table)
#   neighborhoods <- lapply(all_raw_genes,get_global_neighborhood)
#   global_score_table$group_score <- compute_neighborhood_score(global_score_table,neighborhoods)
  
#   per_interactome <- function(interactome){
#     get_subnet(all_raw_genes,interactome,radio)  
#   }
#   subnets <- lapply(interactomes,per_interactome)
#   pseudo_interactome <- unique(do.call("rbind",subnets))
#   group_score <- compute_group_score(global_score_table,pseudo_interactome)
#   global_score_table$group_score <- group_score
  
  global_group_score_table <- compute_global_score_table(scores_frames,inter=F,method=global_score_methods,field="group_score")
  group_global_scores <- global_group_score_table[,global_score_methods]
  colnames(group_global_scores) <- paste("group_",global_score_methods,sep="")
  global_score_table <- cbind(global_score_table,group_global_scores)
  
  to_mix <- cbind(global_score_table[,global_score_methods[1]],global_group_score_table[,global_score_methods[1]])
  global_score_table$mix <- apply(to_mix,1,max)
    
  cluster_score <- function(gene,clusters,single_score){
    gene_cluster <- clusters$clusters[gene]
    neighborhood <- names(which(clusters$clusters==gene_cluster))
    median(single_score[neighborhood])
  }
  clusters <- compute_clusters(rownames(global_score_table),interactomes)
  single_score <- global_score_table[,global_score_methods[1]]
  names(single_score) <- rownames(global_score_table)
  global_score_table$gene_cluster_score <- sapply(rownames(global_score_table),cluster_score,clusters=clusters,single_score=single_score)
  
  ###################################
  
  
  global_score_frame <- as.data.frame(global_score_table)
  global_score_frame$fams <- get_fams_from_score_frame(scores_frames[[1]],rownames(global_score_table))
  
  global_score_table_with_inter <- compute_global_score_table(scores_frames,inter=T,method=global_score_methods)  
  
  ###################################
  
  global_group_score_table_with_inter <- compute_global_score_table(scores_frames,inter=T,method=global_score_methods,field="group_score")
  group_global_scores_with_inter <- global_group_score_table_with_inter[,global_score_methods]
  colnames(group_global_scores_with_inter) <- paste("group_",global_score_methods,sep="")
  global_score_table_with_inter <- cbind(global_score_table_with_inter,group_global_scores_with_inter)
  
  to_mix_with_inter <- cbind(global_score_table_with_inter[,global_score_methods[1]],global_group_score_table_with_inter[,global_score_methods[1]])
  global_score_table_with_inter$mix <- apply(to_mix_with_inter,1,max)
  
#   all_raw_genes <- rownames(global_score_table_with_inter)
#   per_interactome <- function(interactome){
#     v1 <- which(interactome$V1 %in% all_raw_genes)
#     v3 <- which(interactome$V3 %in% all_raw_genes)
#     interactome[union(v1,v3),]
#   }
#   subnets <- lapply(interactomes,per_interactome)
#   pseudo_interactome <- unique(do.call("rbind",subnets))
#   clusters_with_inter <- compute_clusters(rownames(global_score_table_with_inter),pseudo_interactome=pseudo_interactome)
#   global_score_table_with_inter$gene_cluster_score <- sapply(rownames(global_score_table_with_inter),cluster_score,clusters=clusters_with_inter,single_score=global_score_table_with_inter[,global_score_methods[1]])
  global_score_table_with_inter$gene_cluster_score <- rep(0,nrow(global_score_table_with_inter))
  global_score_table_with_inter[rownames(global_score_table),global_score_table$gene_cluster_score]
  
  ###################################
  
  global_score_frame_with_inter <- as.data.frame(global_score_table_with_inter)
  global_score_frame_with_inter$fams <- get_fams_from_score_frame(scores_frames[[1]],rownames(global_score_table_with_inter))
  
  return(list(
    scores_frame_list = scores_frames,
    global_score_table = global_score_table,
    global_score_frame = global_score_frame,
    global_score_table_with_inter = global_score_table_with_inter,
    global_score_frame_with_inter = global_score_frame_with_inter
  ))
  
}

compute_global_score_table <- function(score_frames,inter=F,method="max",field="score"){
  
  # define scored genes
  if(inter){
    scored_genes <- unique(unlist(lapply(score_frames,function(score_frame) return(rownames(score_frame$score_table_with_intermediates)))))  
  } else {
    scored_genes <- unique(unlist(lapply(score_frames,function(score_frame) return(rownames(score_frame$score_table)))))
  }
  
  # all scores in a matrix
  ninteractomes <- length(score_frames)
  score_matrix <- matrix(rep(0,length(scored_genes)*ninteractomes),ncol=ninteractomes)
  rownames(score_matrix) <- scored_genes
  colnames(score_matrix) <- names(interactomes)
  for(i in 1:ninteractomes){
    if(inter){
      score_matrix[,i] <- score_frames[[i]]$score_table_with_intermediates[scored_genes,field]
    } else {
      score_matrix[,i] <- score_frames[[i]]$score_table[scored_genes,field]
    }
  }
  
  # all scores in a matrix
#   group_score_matrix <- matrix(rep(0,length(scored_genes)*ninteractomes),ncol=ninteractomes)
#   rownames(group_score_matrix) <- scored_genes
#   colnames(group_score_matrix) <- paste("group_",names(interactomes),sep="")
#   for(i in 1:ninteractomes){
# #     print(head(score_frames[[i]]$score_table[scored_genes,]))
#     if(inter){
#       group_score_matrix[,i] <- score_frames[[i]]$score_table_with_intermediates[scored_genes,"group_score"]
#     } else {
#       group_score_matrix[,i] <- score_frames[[i]]$score_table[scored_genes,"group_score"]
#     }
#   }
#   
  # compute global score
  global_score <- compute_global_score(score_matrix,method=method)
  if(ncol(global_score)>1){
    score_rank <- order(global_score[,1],decreasing=T)  
  } else {
    score_rank <- order(global_score,decreasing=T)
  }

  score_matrix <- cbind(global_score,score_matrix)
  score_matrix <- score_matrix[score_rank,]
  
  # mount data frame
  score_matrix_frame <- as.data.frame(score_matrix)
  score_matrix_frame[is.na(score_matrix_frame)] <- 0
  
  return(score_matrix_frame)
  
}

compute_global_score <- function(score_matrix,method=NULL){
  
  if(is.null(method)){
    method <- "sum"
  }
  
  nmethods <- length(method)
  
  global_score <- matrix(rep(0,nrow(score_matrix)*nmethods),ncol=nmethods)
  colnames(global_score) <- method
  rownames(global_score) <- rownames(score_matrix)
  
  for(m in 1:nmethods){
    method_function <- get(method[m])
    global_score[,m] <- apply(score_matrix,1,method_function)
  }
 
  return(global_score)
}

compute_group_score <- function(global_score_table,pseudo_interactome){
  cat("computing group score\n")
  cat("...with ",nrow(global_score_table),"genes\n")
  cat("...with ",nrow(pseudo_interactome),"interactions\n")
  raw_genes <- rownames(global_score_table)
  print(head(raw_genes))
  cat("...computing distance matrix\n")
  super_dm <- compute_distance_matrix(pseudo_interactome)
  clean_genes <- raw_genes[raw_genes %in% rownames(super_dm)]
  dm <- super_dm[clean_genes,]# get_gene_distance_matrix(raw_genes,pseudo_interactome) # interactomes[["binding"]])

  get_neighbourhood <- function(row){
    names(which(row>0 & row<=1))
  }
  neighbors <- apply(dm,1,get_neighbourhood)

  single_scores <- global_score_table[,1]
  names(single_scores) <- raw_genes

  group_score <- function(genes){
    if(length(genes)==0){
      return(0)
    } else {
      return(mean(single_scores[genes]))
    }
  }
  raw_group_scores <- sapply(neighbors,group_score)
 
  group_scores <- raw_group_scores[raw_genes]
  names(group_scores) <- raw_genes
  group_scores[is.na(group_scores)] <- 0
  
  return(group_scores)
  
}

compute_neighborhood_score <- function(global_score_table,neighborhoods,column=1){
    
  raw_genes <- rownames(global_score_table)
  
  single_scores <- global_score_table[,column]
  names(single_scores) <- raw_genes

  neighborhood_score <- function(genes){
    if(length(genes)==0){
      return(0)
    } else {
      return(mean(single_scores[genes]))
    }
  }
  raw_neighborhood_scores <- sapply(neighborhoods,neighborhood_score)
 
  neighborhood_scores <- raw_neighborhood_scores[raw_genes]
  names(neighborhood_scores) <- raw_genes
  neighborhood_scores[is.na(neighborhood_scores)] <- 0
  
  return(neighborhood_scores)

  
  
}


no_zero_mean <- function(scores){
  
  no_zero <- which(scores!=0)
  if(length(no_zero)==0){
    score <- 0
  } else {
    score <- mean(scores[no_zero])
  }
  return(score)
  
}

no_zero_min <- function(scores){
  
  no_zero <- which(scores!=0)
  if(length(no_zero)==0){
    score <- 0
  } else {
    score <- min(scores[no_zero])
  }
  return(score)
  
}

no_zero_max <- function(scores){
  
  no_zero <- which(scores!=0)
  if(length(no_zero)==0){
    score <- 0
  } else {
    score <- max(scores[no_zero])
  }
  return(score)
  
}

get_fams_from_score_frame <- function(scores_frame,genes,inter_family="INTER"){
  
  fams <- scores_frame$fam_hits[genes]  
  fams[is.na(fams)] <- inter_family
  names(fams) <- genes
  return(fams)
  
}

