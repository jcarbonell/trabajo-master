
 
#################################################################################################################################################
# UTILS
#################################################################################################################################################


nh_get_neiborhood <- function(node,sif,radio=5){
  
  network <- graph.data.frame(matrix(sif[,c(1,3)],ncol=2))
  
  all_distances <- shortest.paths(network,mode="all")
  
  node_names <- V(network)$name
  
  neighbors <- node_names[which(all_distances[which(node_names==node),]<=radio)]
 
  return(nh_filter_interactions(neighbors,sif))
  
}

nh_filter_interactions <- function(selected_nodes,sif){
  
  left <- which(is.element(sif[,1],selected_nodes))
  right <- which(is.element(sif[,3],selected_nodes))
  
  interactions <- unique(c(left,right))
  
  return(sif[interactions,])
  
}

nh_get_standard_layout <- function(sif){
  
   network <- graph.data.frame(sif[,c(1,3)])
   
   layout <- layout.graphopt(network)
   rownames(layout) <- V(network)$name
   
   return(layout)
   
}

nh_get_standard_layout <- function(sif){
  
   network <- graph.data.frame(sif[,c(1,3)])
   
   layout <- layout.graphopt(network,niterInteger=1000)
   #layout <- layout.fruchterman.reingold(network)
   rownames(layout) <- V(network)$name
   
   return(layout)
   
}

nh_remove_inter_inter_interactions <- function(subnet,family_hits){
  
  left_fams <- family_hits[subnet[,1]]
  left <- which(left_fams!="inter")
  
  right_fams <- family_hits[subnet[,3]]
  right <- which(right_fams!="inter")
  
  subnet[union(left,right),]
  
}

nh_get_nodes <- function(subnet){
  return(unique(c(subnet[,1],subnet[,3])))
}

nh_get_uninteractors <- function(gene_list,interactome){
  
  left <- which(!is.element(gene_list,interactome[,1]))
  right <- which(!is.element(gene_list,interactome[,3]))
  
  gene_list[intersect(left,right)]
  
}

nh_add_pseudo_interactions <- function(subnet,uninteractors){
  
  ui <- data.frame(cbind(uninteractors,"pp",uninteractors))
  names(ui) <- names(subnet)
  rbind(subnet,ui)
  
}

#################################################################################################################################################
# SCORE COMPUTATIONS
#################################################################################################################################################


nh_score <- function(raw_gene_list,sif,radio=5){
    
  # init interactors
  interactors <- unique(c(sif[,1],sif[,3]))
  ninteractors <- length(interactors)
  
  # init super list
  all_raw_genes <- unique(unlist(raw_gene_list))
  super_list <- unique(c(all_raw_genes,interactors))
  nsuper <- length(super_list)
  
  # compute fam hits
  family_names <- names(raw_gene_list)
  nfams <- length(unique(family_names))
  genes_by_family <- sapply(raw_gene_list,length)
  
  super_list_hits <- character(nsuper)
  super_list_nfams <- numeric(nsuper)
  super_list_fams <- list()
    
  for(i in 1:nsuper){
        
    thefams <- names(which(lapply(raw_gene_list,function(fam) is.element(super_list[i],fam))==T))

    if(length(thefams)==0){
      super_list_hits[i] <- "inter"
      super_list_nfams[i] <- 0
      super_list_fams[super_list[i]] <- c() 
    } else {
      super_list_hits[i] <- paste(thefams,collapse="_")
      super_list_nfams[i] <- length(thefams)
      super_list_fams[super_list[i]] <- list(thefams)
    }  
      
  }
  
  names(super_list_hits) <- super_list
  names(super_list_nfams) <- super_list
  
  
  # COMPUTE SCORE ######################################################
  
    # init igraph network
  network <- graph.data.frame(sif[,c(1,3)])
  V(network)$label <- V(network)$name
  node_names <- V(network)$name
  
  # compute node distances
  all_distances <- shortest.paths(network,mode="all")
    
  # compute weights
  steps <- radio+1
  distances <- 0:radio
  raw_weights <- 1-pnorm(radio:(radio*2),mean=radio,sd=2)
  weights <- raw_weights/sum(raw_weights)

  # init working variables
  nfc_node_scores <- numeric(length(super_list))
  names(nfc_node_scores) <- super_list
  nfc_scores_and_distances <- matrix(rep(0,nsuper*steps),ncol=steps)
  rownames(nfc_scores_and_distances) <- super_list
  colnames(nfc_scores_and_distances) <- distances
  
  pnh_node_scores <- numeric(length(super_list))
  names(pnh_node_scores) <- super_list
  pnh_scores_and_distances <- matrix(rep(0,nsuper*steps),ncol=steps)
  rownames(pnh_scores_and_distances) <- super_list
  colnames(pnh_scores_and_distances) <- distances
  

  for(w in 1:nsuper){
    
    i <- which(interactors==super_list[w])[1]

    if(length(i)==0){
      
      node <- super_list[w]
      
      if(is.null(super_list_fams[[node]])){
        node_fams <- c()
        nfc_score <- 0
        pnh_score <- 0
      } else {        
        node_fams <- unique(unlist(super_list_fams[[node]]))        
        nfc_score <- length(node_fams)*weights[1]
        pnh_score <- length(node_fams)*weights[1]      
      }
          
    } else {
      
      node <- interactors[i]
     
      #cat("processing interactor",node,"...\n")
      
      nfc_score <- 0
      pnh_score <- 0
      
      for(j in 1:length(distances)){
        
        distance <- distances[j]
        
        #cat("    distance",distance," -> ")
        
        selected_distance_nodes <- NA
        
        if(distance==0){
        
          if(is.null(super_list_fams[[node]])){
            node_fams <- c()
            nfc_score <- 0
            pnh_score <- 0
          } else {
            node_fams <- unique(unlist(super_list_fams[[node]]))          
            nfc_score <- nfc_score + length(node_fams)*weights[j]
            pnh_score <- pnh_score + length(node_fams)*weights[j]
          }
          
            
        } else {
                  
          selected_distance_nodes <- node_names[which(all_distances[which(V(network)$name==node),]==distance)]
                    
          if(length(selected_distance_nodes)>0){
                                
            # way one
            node_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
            nfc_score <- nfc_score + length(node_fams)*weights[j]
          
            # way two
            allfams <- unlist(super_list_fams[selected_distance_nodes])
            pfams <- numeric(nfams)
            for(k in 1:nfams){
              from_fam <- length(which(allfams==family_names[k]))
              single_p <- (genes_by_family[k]/nsuper)*(length(selected_distance_nodes)/nsuper)
              multi_p <- single_p  ^ from_fam
              pfams[k] <- multi_p          
            }
            
            node_fams <- unique(allfams)
            pnh_score <- pnh_score + sum(1/pfams)*weights[j]
            
          }
        
        }
        
        nfc_scores_and_distances[i,j] <- nfc_score
        pnh_scores_and_distances[i,j] <- pnh_score
          
        #cat(length(selected_distance_nodes),":",length(node_fams)," = ",nfc_score," ", pnh_score,"\n")
                
      }  
              
  #         
  #     if(zero_intermediates & is.null(super_list_fams[[node]])){        
  #       nfc_node_scores[w] <- 0
  #       pnh_node_scores[w] <- 0
  #     } else {
        
  #     }  
    
    }
        nfc_node_scores[w] <- nfc_score
        pnh_node_scores[w] <- pnh_score
    
  }

  scores <- list(
    nfc=nfc_node_scores,
    pnh=pnh_node_scores,
    hits=super_list_hits,
    fams=super_list_fams,
    nfams=super_list_nfams,
    all_distances=all_distances
    )
  
  return(scores)

}

nh_multi_score <- function(raw_gene_list,sif,radio=5,inter_family="INTER"){
    
  # init interactors
  interactors <- unique(c(sif[,1],sif[,3]))
  ninteractors <- length(interactors)
  
  # init super list
  all_raw_genes <- unique(unlist(raw_gene_list))
  super_list <- unique(c(all_raw_genes,interactors))
  nsuper <- length(super_list)
  
  # compute fam hits
  family_names <- names(raw_gene_list)
  nfams <- length(unique(family_names))
  genes_by_family <- sapply(raw_gene_list,length)
  
  super_list_hits <- character(nsuper)
  super_list_nfams <- numeric(nsuper)
  super_list_fams <- list()
    
  # obtain gene fams
  get_gene_fams <- function(gene){
    gene_fams <- names(which(lapply(raw_gene_list,function(fam) is.element(gene,fam))==T))
    if(length(gene_fams)==0){
      gene_fams <- inter_family
    }
    return(gene_fams)
  }
  super_list_fams <- lapply(super_list,get_gene_fams)
  names(super_list_fams) <- super_list
  
  
  # COMPUTE SCORE ######################################################
  
  # init igraph network
  network <- graph.data.frame(sif[,c(1,3)])
  V(network)$label <- V(network)$name
  node_names <- V(network)$name
  
  # compute node distances
  all_distances <- shortest.paths(network,mode="all")
    
  # compute weights
  steps <- radio+1
  distances <- 0:radio
  raw_weights <- 1-pnorm(radio:(radio*2),mean=radio,sd=2)
  weights <- raw_weights/sum(raw_weights)

  # init working variables
  nfc_node_scores <- numeric(length(super_list))
  names(nfc_node_scores) <- super_list
  nfc_scores_and_distances <- matrix(rep(0,nsuper*steps),ncol=steps)
  rownames(nfc_scores_and_distances) <- super_list
  colnames(nfc_scores_and_distances) <- distances
  
  pnh_node_scores <- numeric(length(super_list))
  names(pnh_node_scores) <- super_list
  pnh_scores_and_distances <- matrix(rep(0,nsuper*steps),ncol=steps)
  rownames(pnh_scores_and_distances) <- super_list
  colnames(pnh_scores_and_distances) <- distances
  

  for(w in 1:nsuper){
    
    i <- which(interactors==super_list[w])[1]

    if(length(i)==0){
      
      node <- super_list[w]
      
      if(is.null(super_list_fams[[node]])){
        node_fams <- c()
        nfc_score <- 0
        pnh_score <- 0
      } else {        
        node_fams <- unique(unlist(super_list_fams[[node]]))        
        nfc_score <- length(node_fams)*weights[1]
        pnh_score <- length(node_fams)*weights[1]      
      }
          
    } else {
      
      node <- interactors[i]
     
      #cat("processing interactor",node,"...\n")
      
      nfc_score <- 0
      pnh_score <- 0
      
      for(j in 1:length(distances)){
        
        distance <- distances[j]
        
        #cat("    distance",distance," -> ")
        
        selected_distance_nodes <- NA
        
        if(distance==0){
        
          if(is.null(super_list_fams[[node]])){
            node_fams <- c()
            nfc_score <- 0
            pnh_score <- 0
          } else {
            node_fams <- unique(unlist(super_list_fams[[node]]))          
            nfc_score <- nfc_score + length(node_fams)*weights[j]
            pnh_score <- pnh_score + length(node_fams)*weights[j]
          }
          
            
        } else {
                  
          selected_distance_nodes <- node_names[which(all_distances[which(V(network)$name==node),]==distance)]
                    
          if(length(selected_distance_nodes)>0){
                                
            # way one
            node_fams <- unique(unlist(super_list_fams[selected_distance_nodes]))
            nfc_score <- nfc_score + length(node_fams)*weights[j]
          
            # way two
            allfams <- unlist(super_list_fams[selected_distance_nodes])
            pfams <- numeric(nfams)
            for(k in 1:nfams){
              from_fam <- length(which(allfams==family_names[k]))
              single_p <- (genes_by_family[k]/nsuper)*(length(selected_distance_nodes)/nsuper)
              multi_p <- single_p  ^ from_fam
              pfams[k] <- multi_p          
            }
            
            node_fams <- unique(allfams)
            pnh_score <- pnh_score + sum(1/pfams)*weights[j]
            
          }
        
        }
        
        nfc_scores_and_distances[i,j] <- nfc_score
        pnh_scores_and_distances[i,j] <- pnh_score
          
        #cat(length(selected_distance_nodes),":",length(node_fams)," = ",nfc_score," ", pnh_score,"\n")
                
      }  
              
  #         
  #     if(zero_intermediates & is.null(super_list_fams[[node]])){        
  #       nfc_node_scores[w] <- 0
  #       pnh_node_scores[w] <- 0
  #     } else {
        
  #     }  
    
    }
        nfc_node_scores[w] <- nfc_score
        pnh_node_scores[w] <- pnh_score
    
  }

  scores <- list(
    nfc=nfc_node_scores,
    pnh=pnh_node_scores,
    hits=super_list_hits,
    fams=super_list_fams,
    nfams=super_list_nfams,
    all_distances=all_distances
    )
  
  return(scores)

}
  
#################################################################################################################################################
# CLASSIFY
#################################################################################################################################################

  
nh_kmeans_classify <- function(score){
 
  classif <- kmeans(as.matrix(score),quantile(score,c(0.1,0.99)))
  
  dispersion <- sd(score[names(classif$cluster[classif$cluster==2])])
  
  breakpoint <- classif$centers[2] - dispersion
  
  selected=score[score>=breakpoint]
  unselected=score[score<breakpoint]
  
  list(score=score, breakpoint=breakpoint, dispersion=dispersion,selected=selected, selected_genes=names(selected), unselected=unselected, unselected_genes=names(unselected))
  
}

nh_quantile_classify <- function(score,quantile){
 
  breakpoint <- quantile(score,quantile)
  
  selected=score[score>=breakpoint]
  unselected=score[score<breakpoint]
  
  list(score=score, breakpoint=breakpoint,selected=selected, selected_genes=names(selected), unselected=unselected, unselected_genes=names(unselected))
  
}
  
nh_paint_classification <- function(classif, title=NULL){
  
  all <- c(classif$selected,classif$unselected)
  
  if(is.null(title)){
    title <- "Score classification"
  }
  h <- hist(all,50,probability=T,main=title,xlab="score")
  
  selected_ratio <- length(classif$selected)/length(classif$unselected)

  if(length(classif$unselected)>0){
    du <- density(classif$unselected)
    lines(du$x,du$y*(1-selected_ratio),col="blue")    
  }
  if(length(classif$selected)>0){
    ds <- density(classif$selected)
    lines(ds$x,ds$y*selected_ratio,col="red")
  }    
  legend((max(h$mids)-min(h$mids))*0.75,max(h$density),legend=c("selected","unselected"),col=c("red","blue"),lwd=1)
  
}
  
  
  
