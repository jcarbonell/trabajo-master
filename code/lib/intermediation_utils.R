
compute_intermediation <- function(origin,network,vertex_indexes=NULL,all_genes=NULL){
  
  # vertex indexes
  if(is.null(vertex_indexes)){
    vertex_indexes <- get_vertex_indexes(network)
  }
  
  # all genes
  if(is.null(all_genes)){
    all_genes <- names(vertex_indexes)
  }
  
  # undescribed origin
  if(is.na(vertex_indexes[origin])){
    
    # empty intermediation
    intermediation <- rep(0,length(all_genes))
    names(intermediation) <- all_genes
    
  } else {
    
    # all origin shortest path from
    all_sp <- get.shortest.paths(network,from=vertex_indexes[origin])
    
    # gene intermediation
    raw_intermediation <- table(unlist(all_sp))    
    intermediation <- raw_intermediation[as.character(vertex_indexes[all_genes])]
    intermediation[is.na(intermediation)] <- 0
    names(intermediation) <- all_genes
    intermediation <- intermediation/length(all_sp)
      
  }
  
  return(intermediation)
  
}

process_intermediation <- function(interactome,outfile=NULL){
  
  cat(">>>>>>>",outfile,"\n")
  
  interactome <- interactome[which(interactome[,1]!=interactome[,2]),]
  network <- graph.data.frame(interactome[,c(1,3)])  
  network <- simplify(network)
  
  all_genes <- unique(c(interactome[,1],interactome[,3]))
  vertex_indexes <- get_vertex_indexes(network)
  
  intermediation <- do.call("rbind",lapply(all_genes,compute_intermediation,network=network,vertex_indexes=vertex_indexes,all_genes=all_genes))
  rownames(intermediation) <- all_genes
  
  if(!is.null(outfile)){
    write.table(format(digits=5,intermediation),file=outfile,sep="\t",quote=F)  
  } 
  
  return(intermediation)
  
}

# vertex info
get_vertex_indexes <- function(network){ 
  min_index <- 0
  if(length(get.vertex.attribute(network,"name",index=0))==0) {
    min_index <- 1
  }
  max_index <- min_index + length(V(network)) - 1
  vertex_indexes <- seq(min_index,max_index)
  names(vertex_indexes) <- get.vertex.attribute(network,"name",index=vertex_indexes)
  return(vertex_indexes)
}
