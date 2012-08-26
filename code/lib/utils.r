
load_family_genes <- function(group_names,gene_list_files,hasHeader=F){
  
  ngroups <-length(gene_list_files)
  raw_gene_list <- list()
  for(i in 1:ngroups){
    if(file.info(gene_list_files[i])$size==0){
      raw_gene_list[[group_names[i]]] <- data.frame()
    } else {
      raw_gene_list[[group_names[i]]] <- read.delim(gene_list_files[i],stringsAsFactors=F,header=hasHeader,sep="\t")[,1]
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

get_meta_interactome <- function(interactomes){  
  return(unique(do.call("rbind",interactomes)))
}

get_subnets <- function(genes,interactomes,ninter){
  
  subnets <- list()
  
  for(i in 1:length(interactomes)){
    subnets[[ names(interactomes)[i] ]] <- get.all.shortest.paths.Luz(interactomes[[i]],genes,ninter)  
    gc()
  }
  
  return(subnets)
  
}

get_subnet <- function(genes,interactome,ninter,verbose=F){  
  return(get.all.shortest.paths.Josete(interactome,genes,ninter,verbose=verbose))
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

describe_gene_series <- function(genes_series,interactome){
 per_serie <- function(genes,interactome){
  get.all.shortest.paths.Josete(interactome,genes,numinterm=5)  
  }
  all_subnets <- do.call("rbind",lapply(genes_series,per_serie,interactome=interactome))
  describe_interactome(all_subnets)
}

convert_interactome_ids <- function(interactome_file,xref_file,outfile,id_col=1){
  
  interactome <- load_interactome(interactome_file)
  xref <- read.table(file=xref_file,header=T, sep = "\t", quote = "",stringsAsFactors = FALSE, as.is=TRUE)
  rownames(xref) <- xref[,1]
  converted_interactome <- cbind(xref[interactome[,1],2],interactome[,2],xref[interactome[,3],2])
  write.table(converted_interactome,file=outfile,sep="\t",row.names=F,col.names=F,quote=F)
  
}

get.igraph <- function(interactome=NULL,interactome.file=NULL){
  if(is.null(interactome)){
    my.graph <- read.table(interactome.file, quote=NULL, header=FALSE, sep="\t", stringsAsFactors=FALSE)[,c(1,3)]  
    my.graph <- my.graph[which(my.graph[,1]!=my.graph[,2]),]
  } else {
    my.graph <- interactome[which(interactome[,1]!=interactome[,3]),][,c(1,3)]
  }
  my.igraph <- simplify(graph.data.frame(my.graph, directed=FALSE))
  return(my.igraph)
}

get_interactome_params <- function(interactome){
    
  my.graph <- interactome[,c(1,3)]
  my.graph <- my.graph[which(my.graph[,1]!=my.graph[,2]),]
  my.igraph <- simplify(graph.data.frame(my.graph, directed=FALSE))
  
  gene_degree <- degree(my.igraph)
  gene_betweenness <- betweenness(my.igraph)
  gene_closeness <- closeness(my.igraph)
  gene_burt <- constraint(my.igraph)
  
  params <- data.frame(degree=gene_degree,betweenness=gene_betweenness,closeness=gene_closeness,burt=gene_burt)
  rownames(params) <- get.vertex.attribute(my.igraph,"name")
  
  return(params)
  
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
  if(length(list)>1){
    for(i in list){
      out <- c(out, get.shortest.paths(my.igraph,from=(which(V(my.igraph)$name==i)-1),to=destiny,mode = "all"))
    }
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

add_hgnc_symbol <- function(ens2hgnc_file,data_file,ens_col,outfile=NULL){
  
  if(is.null(outfile)){
    outfile <- paste(data_file,".hgnc",sep="") 
  }
  
  raw <- read.delim(ens2hgnc_file,sep="\t",header=T,stringsAsFactors=F)
  ens2hgnc <- raw$Associated.Gene.Name
  names(ens2hgnc) <- raw$Ensembl.Gene.ID
  
  data <- read.delim(data_file,sep="\t",header=F,stringsAsFactors=F)
  
  data$hgnc <-ens2hgnc[data[,1]]
  
  write.table(data,file=outfile,quote=F,row.names=F,col.names=F,sep="\t")
}


load_known_disease_genes_series <- function(annot_file,nfams){
  
  size <- length(nfams)
  
  dannot <- read.delim(annot_file,sep="\t",header=T,stringsAsFactors=F)
  dsize <- by(dannot,dannot$disease,nrow,simplify=T)
  
  fam_disease_genes <- list()
  selected_diseases <- character(size)

  for(i in 1:size){
    available_indexes <- which(dsize>=nfams[i])
    if(length(available_indexes)>0){
      selected_index <- sample(available_indexes,1)
      selected_disease <- names(dsize)[selected_index]
      selected_diseases[i] <- selected_disease
      fam_disease_genes[[i]] <- sample(dannot$gene[which(dannot$disease==selected_disease)],nfams[i])
    } else {
      fam_disease_genes[[i]] <- NULL
    }
  } 
  
  return(list(
    fam_disease_genes=fam_disease_genes,
    selected_diseases=selected_diseases
  ))
  
}

get_gene_distance_matrix <- function(genes,interactome){
  
  ngenes <- length(genes)
  subnet <- get.all.shortest.paths.Josete(interactome,genes,numinterm=5)
  subnet_genes <- c(subnet$V1,subnet$V3)
  
  distance_matrix <- matrix(rep(Inf,ngenes*ngenes),ncol=ngenes)
  colnames(distance_matrix) <- genes
  rownames(distance_matrix) <- genes
  
  if(!is.null(subnet)){
    raw_distance_matrix <- compute_distance_matrix(subnet)
    for(i in 1:ngenes){
      for(j in 1:ngenes){
        if((genes[i] %in% subnet_genes) & (genes[j] %in% subnet_genes)){        
          distance <- raw_distance_matrix[genes[i],genes[j]]
          distance_matrix[genes[i],genes[j]] <- distance
          distance_matrix[genes[j],genes[i]] <- distance
        }
      }
    }

  }
  
  return(distance_matrix)
}


get_interactome_genes <- function(interactome){
    return(unique(c(interactome[,1],interactome[,3])))
}
