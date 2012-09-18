
final_score <- function(node,distance_matrix,super_list_fams,family_weights=NULL,inter_family="INTER"){
  
  score <- 0
  
  if(node %in% colnames(distance_matrix)){
  
    node_distances <- distance_matrix[node,]
    
    fam_list <- unlist(unique(super_list_fams))    
    fam_list <- setdiff(fam_list,inter_family)

    if(is.null(family_weights)){
      family_weights <- rep(1,length(fam_list))
      names(family_weights) <- fam_list
    }
              
    for(i in 1:length(fam_list)){
              
      all_fam_genes <- names(which(lapply(super_list_fams,function(x) fam_list[i] %in% x)==T))
      fam_gene_distances <- node_distances[all_fam_genes]
           
      # direct intersection term
      fix_term <- 0
      if(node %in% all_fam_genes){
        fix_term <- family_weights[fam_list[i]]
       } 
     
      # interaction term
      interaction_term <- 0
      if(!(node %in% all_fam_genes)){
        selected_gene_distances <- fam_gene_distances[which(fam_gene_distances>=quantile(fam_gene_distances,.85,na.rm=T))]
        interaction_term <- mean(selected_gene_distances,na.rm=T)
      }
      
      score <- score + fix_term + interaction_term
      
      #cat(node,":",fam_list[i],"     ",fix_term," + ", interaction_term,"\n")
      
    }

  } else {
    score <- NA
  }
  
  return(list(
    score=score
  ))
  
}








# 
# default_score3c <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
# 
#   node_fam <- unique(unlist(super_list_fams[node]))
#   
#   score <- 0
#   
#   if(node %in% colnames(distance_matrix)){
#   
#     node_distances <- distance_matrix[node,]
#   
#     interactors <- node_distances#[which(node_distances>0 & node_distances<=5)]
#   
#     interactor_names <- names(interactors)
#     
#     fam_list <- unlist(unique(super_list_fams))
#     fam_list <- setdiff(fam_list,node_fam)
#     fam_list <- setdiff(fam_list,inter_family)
#               
#     for(i in 1:length(fam_list)){
#       all_fam_genes <- names(which(lapply(super_list_fams,function(x) fam_list[i] %in% x)==T))
#       fam_distances <- node_distances[all_fam_genes]
#       
# #       score <- score +min(fam_distances,na.rm=T)
#     #score <- score + mean(fam_distances,na.rm=T)
#        #score <- score + sum(fam_distances[which(fam_distances>=quantile(fam_distances,.5,na.rm=T))],na.rm=T)
#     score <- score + sum(fam_distances[which(fam_distances>=quantile(fam_distances,.95,na.rm=T))],na.rm=T)
#     }
# 
#   } else {
#     score <- NA
#   }
#   
#   return(list(
#     scores_by_distance=rep(0,6),
#     score=score
#     ))
#   
# }
# 
# default_score3d <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
#   
#   weight <- function(distance){
#     1-pnorm(distance*2,mean=5,sd=1.1)
#   }
#   
#   node_fam <- unique(unlist(super_list_fams[node]))
#   
#   score <- 0
#   
#   if(node %in% colnames(distance_matrix)){
#   
#     node_distances <- distance_matrix[node,]
#   
#     interactors <- node_distances[which(node_distances>0 & node_distances<=5)]
#   
#     interactor_names <- names(interactors)
#     
#     fam_list <- unlist(unique(super_list_fams))
#     fam_list <- setdiff(fam_list,node_fam)
#     fam_list <- setdiff(fam_list,inter_family)
#       
#     all_fam_genes <- c()
#     for(i in 1:length(fam_list)){
#       fam_genes <- names(which(lapply(super_list_fams,function(x) fam_list[i] %in% x)==T))
#       all_fam_genes <- c(all_fam_genes,fam_genes)
#     }
#     all_fam_genes <- unique(all_fam_genes)
#     
#     fam_distances <- node_distances[all_fam_genes]      
# #     score <- score + 5-min(fam_distances,na.rm=T)
# #     score <- score + sum(fam_distances,na.rm=T)
#     score <- score + sum(fam_distances[which(fam_distances>=quantile(fam_distances,.95,na.rm=T))],na.rm=T)
# 
#     
#   } else {
#     score <- NA
#   }
#   
#   return(list(
#     scores_by_distance=rep(0,6),
#     score=score
#     ))
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# score_kernel <- function(node,distance_matrix,weights,distance_priors,super_list_fams,ACCOUNT_OWN_FAMILY=F,SCALE_BY_FAMILY_SIZE=F,SCALE_BY_DISTANCE_PRIORS=F,inter_family="INTER"){
#   #
#   # Computes gene score
#   #
#   # parameters:
#   
#   #  node: gene to score
#   #  distance_matrix: interactome distance matrix
#   #  weights: Weights of scores by distance
#   #  distance_priors: glogal density distribution of gene-to-gene distances
#   #  super_list_fams: Family ids of selected and intermediates gene
#   #  
#   #  ACCOUNT_OWN_FAMILY: Define whether the own node family must be filtered 
#   #  SCALE_BY_FAMILY_SIZE: Scale every family contribution by family size
#   #  SCALE_BY_DISTANCE_PRIORS: Scale every distance contribution by distance priors
#   #
#   # return:
#   # 
#   #  gene score
#   #
#   #
#   
#   family_freqs <- table(as.character(super_list_fams))
#   
#   scores_by_distance <- numeric(length(weights))
#   
#   distance <- 0
#   distance_prior <- 1
#   
#   for(j in 1:length(weights)){
#     
#     if(distance==0){
#       selected_distance_nodes <- node
#       distance_prior <- distance_priors[1]
#     } else {
#       node_distances <- distance_matrix[node,]
#       selected_distance_nodes <- names(node_distances[which(node_distances==distance)])
#       distance_prior <- distance_priors[distance]
#     }
#     
#     node_fam <- unique(unlist(super_list_fams[node]))
#     
#     selected_distance_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
#     selected_distance_fams <- setdiff(selected_distance_fams,inter_family)
#     if(distance>0 & ACCOUNT_OWN_FAMILY==F){
#       selected_distance_fams <- setdiff(selected_distance_fams,node_fam)
#     }
#     
#     if(is.null(selected_distance_fams)){
#       scores_by_distance[j] <- 0
#     } else {
#       selected_distance_fams_freqs <- table(unlist(super_list_fams[selected_distance_nodes]))[selected_distance_fams]
#       selected_distance_fams_sizes <- family_freqs[selected_distance_fams]
#       scores_by_distance[j] <- sum(selected_distance_fams_freqs/selected_distance_fams_sizes) * weights[j]
#     }    
#     
#     distance <- distance+1
#     
#   }
#   
#   score <- sum(scores_by_distance) - 1*weights[j]
#   
#   return(list(
#     scores_by_distance=scores_by_distance,
#     score=score
#     ))
#   
# }
# 
# 
# default_score <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
#   
#   scores_by_distance <- numeric(length(weights))
#   
#   distance <- 0
#   
#   per_distance <- function(j){
#     
#       distance = j-1
#       
#       if(distance==0){
#         selected_distance_nodes <- node
#       } else {
#         node_distances <- distance_matrix[node,]
#         selected_distance_nodes <- names(node_distances[which(node_distances==distance)])
#       }
#         
#       selected_distance_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
#       selected_distance_fams <- setdiff(selected_distance_fams,inter_family)
#   
#       if(is.null(selected_distance_fams)){
#           scores_by_distance <- 0
#       } else {
#           scores_by_distance <- (length(selected_distance_fams)*weights[j])                   
#       }
#       
#       return(scores_by_distance)
#         
#   }
# 
#   scores_by_distance <- sapply(1:length(weights),per_distance)
#   
#   return(list(
#     scores_by_distance=scores_by_distance,
#     score=sum(scores_by_distance)
#   ))
#   
# }
# 
# default_score2 <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
#   
#   scores_by_distance <- numeric(length(weights))
#   
#   distance <- 0
#   
#   per_distance <- function(j){
#     
#       distance = j-1
#                   
#       if(distance==0){
#         selected_distance_nodes <- node
#       } else {
#         node_distances <- distance_matrix[node,]
#         selected_distance_nodes <- names(node_distances[which(node_distances==distance)])
#       }
#         
#       selected_distance_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
#       selected_distance_fams <- setdiff(selected_distance_fams,inter_family)      
#       if(distance>0){
#         node_fam <- unique(unlist(super_list_fams[node]))
#         selected_distance_fams <- setdiff(selected_distance_fams,node_fam)
#       }
#   
#       if(is.null(selected_distance_fams)){
#           scores_by_distance <- 0
#       } else {                
#           scores_by_distance <- (length(selected_distance_fams)*weights[j])
#       }
#   
#       return(scores_by_distance)
# 
#   }
# 
#   scores_by_distance <- sapply(1:length(weights),per_distance)
#   
#   score <- sum(scores_by_distance) - weights[1]
#         
#   return(list(
#     scores_by_distance=scores_by_distance,
#     score=score
#     ))
#   
# }
# 
# default_score3 <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
#   
#   scores_by_distance <- numeric(length(weights))
#   
#   distance <- 0
#   distance_prior <- 1
#   
#   per_distance <- function(j){
#     
#     distance = j-1
#     
#     if(distance==0){
#       selected_distance_nodes <- node
#       distance_prior <- 1
#     } else {
#       node_distances <- distance_matrix[node,]
#       selected_distance_nodes <- names(node_distances[which(node_distances==distance)])
#       distance_prior <- distance_priors[distance]
#     }
#     
#     node_fam <- unique(unlist(super_list_fams[node]))
#     
#     selected_distance_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
#     selected_distance_fams <- setdiff(selected_distance_fams,inter_family)
#     if(distance>0){
#       node_fam <- unique(unlist(super_list_fams[node]))
#       selected_distance_fams <- setdiff(selected_distance_fams,node_fam)
#     }
#         
#     if(is.null(selected_distance_fams)){
#       scores_by_distance <- 0
#     } else {                
# #       scores_by_distance <- (length(selected_distance_fams)*weights[j])/distance_prior
#       scores_by_distance <- (length(selected_distance_fams)*weights[j])
# 
#     }
#     
#     return(scores_by_distance)
#     
#   }
#   
#   scores_by_distance <- sapply(1:length(weights),per_distance)
#   
#   return(list(
#     scores_by_distance=scores_by_distance,
#     score=sum(scores_by_distance)
#     ))
#   
# }
# 
# default_score3b <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
#   
#   weight <- function(distance){
#     1-pnorm(distance*2,mean=5,sd=1.1)
#   }
#   
#   node_fam <- unique(unlist(super_list_fams[node]))
#   
#   score <- length(node_fam)*weight(0)
#   
#   if(node %in% colnames(distance_matrix)){
#   
#     node_distances <- distance_matrix[node,]
#   
#     interactors <- node_distances[which(node_distances>0 & node_distances<=5)]
#     interactor_names <- names(interactors)
#       
#     for(i in 1:length(interactors)){
#      
#       selected_distance_fams <- unique(unlist(super_list_fams[interactor_names[i]]))
#       selected_distance_fams <- setdiff(selected_distance_fams,node_fam)
#       selected_distance_fams <- setdiff(selected_distance_fams,inter_family)
#       
#       if(!is.null(selected_distance_fams)){        
#         score <- score + as.numeric(length(selected_distance_fams)*weight(interactors[i]))
#       }
#         
#     }
# 
#   }
#   
#   return(list(
#     scores_by_distance=rep(0,6),
#     score=score
#     ))
#   
# }
# 
# 
# 
# 
# default_score4 <- function(node,distance_matrix,weights,distance_priors,super_list_fams,inter_family="INTER"){
#   
#   family_freqs <- table(as.character(super_list_fams))
#   
#   scores_by_distance <- numeric(length(weights))
#   
#   distance <- 0
#   distance_prior <- 1
#   
#   for(j in 1:length(weights)){
#     
#     if(distance==0){
#       selected_distance_nodes <- node
#       distance_prior <- distance_priors[1]
#     } else {
#       node_distances <- distance_matrix[node,]
#       selected_distance_nodes <- names(node_distances[which(node_distances==distance)])
#       distance_prior <- distance_priors[distance]
#     }
#     
#     node_fam <- unique(unlist(super_list_fams[node]))
#     
#     selected_distance_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
#     selected_distance_fams <- setdiff(selected_distance_fams,inter_family)
#     if(distance>0){
#       selected_distance_fams <- setdiff(selected_distance_fams,node_fam)
#     }
#     
#     if(is.null(selected_distance_fams)){
#       scores_by_distance[j] <- 0
#     } else {
#       selected_distance_fams_freqs <- table(unlist(super_list_fams[selected_distance_nodes]))[selected_distance_fams]
#       selected_distance_fams_sizes <- family_freqs[selected_distance_fams]
#       scores_by_distance[j] <- sum(selected_distance_fams_freqs/selected_distance_fams_sizes) * weights[j]
#     }    
#     
#     distance <- distance+1
#     
#   }
#   
#   score <- sum(scores_by_distance) - 1*weights[j]
#         
#   return(list(
#     scores_by_distance=scores_by_distance,
#     score=score
#     ))
#   
# }
# 
# 
