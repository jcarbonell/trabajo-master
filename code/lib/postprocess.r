
compute_clusters <- function(all_raw_genes,interactomes=NULL,pseudo_interactome=NULL,maxdist=0,radio=5){
  
  if(is.null(pseudo_interactome)){
    per_interactome <- function(interactome){
      get_subnet(all_raw_genes,interactome,radio)  
    }
    subnets <- lapply(interactomes,per_interactome)
    pseudo_interactome <- unique(do.call("rbind",subnets))
    subnet <- get_subnet(all_raw_genes,pseudo_interactome,maxdist,verbose=F)
    network <- graph.data.frame(subnet[,c(1,3)])
  } else {
    network <- graph.data.frame(pseudo_interactome[,c(1,3)])
  }
    
  subnet_clusters <- clusters(network)
  
  names(subnet_clusters$membership) <- V(network)$name
  arg_clusters <- subnet_clusters$membership[all_raw_genes]
  names(arg_clusters) <- all_raw_genes
  arg_clusters[is.na(arg_clusters)] <- -1
  
  nout <- length(which(arg_clusters==-1))
  nin <- length(all_raw_genes) - nout 
  
#   cat("Found",nin,"genes in",subnet_clusters$no,"clusters (and",nout,"genes ungroupped)\n")
  
  return(list(
    pseudo_interactome=pseudo_interactome,
    subnet_clusters=subnet_clusters,
    clusters=arg_clusters,
    nin=nin,
    nout=nout
  ))  
  
}

