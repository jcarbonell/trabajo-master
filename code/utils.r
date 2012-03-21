
load_family_genes <- function(group_names,gene_list_files){
  
  ngroups <-length(gene_list_files)
  raw_gene_list <- list()
  for(i in 1:ngroups){
    if(file.info(gene_list_files[i])$size==0){
      raw_gene_list[[group_names[i]]] <- data.frame()
    } else {
      raw_gene_list[[group_names[i]]] <- read.table(gene_list_files[i],stringsAsFactors=F,header=F)[,1]
    }
  }

  return(raw_gene_list)
  
}


load_interactomes <- function(interactomes_path){
  
  interactomes <- list()
   
  interactomes[["binding"]] <- load_interactome(paste(interactomes_path,"/binding_hgnc.sif",sep=""))
  interactomes[["functional"]] <- load_interactome(paste(interactomes_path,"/functional_hgnc.sif",sep=""))
  interactomes[["ptmod"]] <- load_interactome(paste(interactomes_path,"/ptmod_hgnc.sif",sep=""))
  interactomes[["regulation"]] <- load_interactome(paste(interactomes_path,"/regulation_hgnc.sif",sep=""))
  interactomes[["textmining"]] <- load_interactome(paste(interactomes_path,"/textmining_hgnc.sif",sep=""))
  #interactomes[["fake"]] <- load_interactome(paste(interactomes_path,"/fake_interactome.sif",sep=""))

  return(interactomes)
  
}

load_interactome <- function(interactome_file){
    interactome <- read.table(file=interactome_file, header = FALSE, sep = "\t", quote = "",stringsAsFactors = FALSE, as.is=TRUE)
    return(interactome)
}

get_subnets <- function(genes,interactomes,ninter){
  
  subnets <- list()
  
  for(i in 1:length(interactomes)){
    subnets[[ names(interactomes)[i] ]] <- get.all.shortest.paths.Luz(interactomes[[i]],genes,ninter)  
    gc()
  }
  
  return(subnets)
  
}

get_subnet <- function(genes,interactome,ninter){  
  return(get.all.shortest.paths.Josete(interactome,genes,ninter))
}

compute_distance_matrix <- function(sif){
  
  # init igraph network
  network <- graph.data.frame(sif[,c(1,3)])
  V(network)$label <- V(network)$name
  
  # compute node distances
  distance_matrix <- shortest.paths(network,mode="all")
  colnames(distance_matrix) <- V(network)$name
  rownames(distance_matrix) <- V(network)$name
  
  return(distance_matrix)
  
}

describe_interactome <- function(interactome,plot=T,max_dist=15){
    
  distance_matrix <- compute_distance_matrix(interactome)
  
  per_gene <- function(gene_distances,max_dist=15){
    neighbors <- numeric(max_dist)    
    for(i in 1:15){    
       neighbors[i] <- length(which(gene_distances==i))
    }
    return(neighbors)
  }
  
  neighbors <- t(apply(distance_matrix,1,per_gene))
  dinter_summary <- apply(neighbors,2,summary)
  
  if(plot){
    boxplot(neighbors,main="Number of genes by distance")
    lines(spline(dinter_summary["Median",],n=100),col="red")
  }
  
  return(dinter_summary)
  
}

convert_interactome_ids <- function(interactome_file,xref_file,outfile,id_col=1){
  
  interactome <- load_interactome(interactome_file)
  xref <- read.table(file=xref_file,header=T, sep = "\t", quote = "",stringsAsFactors = FALSE, as.is=TRUE)
  rownames(xref) <- xref[,1]
  converted_interactome <- cbind(xref[interactome[,1],2],interactome[,2],xref[interactome[,3],2])
  write.table(converted_interactome,file=outfile,sep="\t",row.names=F,col.names=F,quote=F)
  
}


get.all.shortest.paths.Josete <- function(sif, list, numinterm,verbose=F){

  # clean  files
  list <- clean.list.exits(sif, list)
  sif <- clean.sif(sif)
  
  # get igraph
  my.igraph <- graph.data.frame(sif[,c(1,3)], directed=FALSE, vertices=NULL) # get an igraph  
  destiny <- (which(V(my.igraph)$name %in% list)-1)
  
  # get shortest paths
#   per_node <- function(i){
#     return(get.all.shortest.paths(my.igraph,from=(which(V(my.igraph)$name==i)-1),to=destiny,mode = "all"))
#   }
#   out2 <- lapply(list,per_node)
  

  out <- NULL
  for(i in list){
    out <- c(out, get.shortest.paths(my.igraph,from=(which(V(my.igraph)$name==i)-1),to=destiny,mode = "all"))
  }
  
  if(length(out)>1){
    
    paths <- extract.paths.from.list(out, my.igraph, numinterm)

    if(verbose){
      print(paste(numinterm, " intermediates allowed"))
    }
    # get subnet from nodes in shortest path
    subnet <- get.subnet.from.nodes.in.paths(paths, sif, my.igraph)

    if(verbose){
      get.mydata.report(subnet, list)
    }
  } else {
    subnet <- unlist(out)
  }
  return(subnet)
}




