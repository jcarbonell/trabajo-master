#####################################################################################################################################
##                                                                                                                                 ##
## Luz Garcia Alosno                                                                                                               ##
## 21/07/2011                                                                                                                      ##
##                                                                                                                                 ##
## SCRIPT to extract SHORTEST PATHS from a list of nodes using a sif file                                                          ##
## input: 1. a sif file, 2. a node list, · number of intermediates allowed                                                         ##
## output: sif file wich represents your paths                                                                                     ##
##                                                                                                                                 ##
## use library(igraph)                                                                                                             ##
##                                                                                                                                 ##
######################################################################################################################################


library(igraph)          

# my functions
clean.list.exits <- function(sif, list,verbose=F){
  if(verbose){
    print("Removing duplicates and nodes that don't match with your sif file from your list")
    print(paste("- initial list length: ", length(list)))
  }
  interactors <- unique(c(sif[,1], sif[,3]))
  list <-unique( list[ list[] %in% interactors])
  if(verbose){
    print(paste("- final list length: ", length(list)))
  }
  return(list)
}

clean.sif <- function(sif){
  sif <- sif[!is.na(sif[,1]),]
  sif <- sif[!is.na(sif[,3]),]
  return(sif)
}

extract.paths.from.list <- function(out, my.igraph, numinterm){
  #filter paths by length
  path.lengths <- cbind(c(1:length(out)),  sapply(out, length))  
  c.path.lengths <- path.lengths[path.lengths[,2] <= (numinterm + 2),]
  out.filtered <- out[c.path.lengths[,1]]

  # extract paths
  size <- sort(unique(path.lengths[,2]), decreasing=TRUE)[1]
  paths <-  NULL
  for (num in 1:size){
    a <-sapply(out.filtered, "[", num)
    paths <- cbind(paths, a)
  }
  return(paths)
}


get.subnet.from.nodes.in.paths <- function(paths, sif, my.igraph){
  # transform paths in non-redundant sif(two columns) with all interactions between the nodes in shortest paths
  nodes.paths <- NULL
  for( k in 1:length(paths[1,]) ){
    nodes.paths <- c(nodes.paths, paths[,k])
  }
  nodes.paths <- nodes.paths[!is.na(nodes.paths)] 
  subnet <- sif[sif[,1] %in% V(my.igraph)[nodes.paths]$name,]
  subnet <- subnet[subnet[,3] %in% V(my.igraph)[nodes.paths]$name ,] 
  return(subnet)
}


get.mydata.report <- function(subnet, list){
    print(paste("your list contains:", length(list), " nodes"))
    print(paste("your subnet contains:", dim(subnet)[1], " interactions"))
    interactors <- unique(c(subnet[,1],subnet[,3]))
    print(paste("your subnet contains:", length(interactors), " interactors"))
    interactors.list <- interactors[interactors[] %in% list]
    print(paste("your subnet contains:", length(interactors.list), " interactors from list"))
}


get.all.shortest.paths.Luz <- function(sif, list, numinterm,verbose=F){
  # clean  files
  list <- clean.list.exits(sif, list)
  sif <- clean.sif(sif)
    
  # get igraph
  my.igraph <- graph.data.frame(sif[,c(1,3)], directed=FALSE, vertices=NULL) # get an igraph
  if(verbose){
    print("RESULTS")
  }
  
  # get shortest paths
  out <- NULL
  for(i in list){
    out <- c(out, get.all.shortest.paths(my.igraph, from= (which(V(my.igraph)$name %in% i)-1),
                                  to=(which(V(my.igraph)$name %in% list)-1), mode = "all"))
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
