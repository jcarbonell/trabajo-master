

source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/evaluation.r")
source("lib/postprocess.r")
source("lib/random_walk.r")
source("lib/intermediation_utils.R")

home_path <- getwd()

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))



bin <- get.igraph(interactomes$binding)
bin_ix <- get_vertex_indexes(bin)

fun <- get.igraph(interactomes$functional)
fun_ix <- get_vertex_indexes(fun)









fake <- load_interactome("../interactomes/fake_interactome.sif")

#
sif <- interactomes$binding
network <- graph.data.frame(sif[,c(1,3)])
fake_net <- graph.data.frame(fake[,c(1,3)])

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

# 
g <- "ATG5"
h <- "DAGLB"

all_genes <- unique(c(sif[,1],sif[,3]))
sub_inter <- get_subnet(genes=sample(all_genes,500),sif,ninter=5)
subnet <- graph.data.frame(sub_inter[,c(1,3)])
vertex_indexes <- get_vertex_indexes(subnet)

cont <- 0
for(g in names(vertex_indexes)){ 
  pp <- get.all.shortest.paths(subnet,from=vertex_indexes[g])
  tt <- table(unlist(pp))
  tt_clean <- tt[as.character(vertex_indexes)]
  names(tt_clean) <- vertex_indexes
  tt_clean[is.na(tt_clean)] <- 0
  tt_norm <- 1 - tt_clean/length(pp)
  names(tt_norm) <- names(vertex_indexes)
  cont <- cont +1
  print(cont)
}
  
# system.time(pp <- lapply(vertex_indexes,function(x) get.all.shortest.paths(network,from=x)))

p_dists <- numeric(length(vertex_indexes))
names(p_dists) <- vertex_indexes
p_count <- numeric(length(vertex_indexes))
names(p_count) <- vertex_indexes
for(i in 1:length(pp)){
    path <- pp[[i]]
    if(length(path)>1){
      dists <- 2:length(path)
      path <- path[dists]
      p_dists[path] <- p_dists[path] + dists
      p_count[path] <- p_count + 1
    }
}

p_distances <- p_dists / p_count^2



n <- length(vertex_indexes)

pp <- get.all.shortest.paths(network)

p_dists <- matrix(0,nrow=n,ncol=n)
colnames(p_dists) <- vertex_indexes
rownames(p_dists) <- vertex_indexes

p_count <- matrix(0,nrow=n,ncol=n)
colnames(p_count) <- vertex_indexes
rownames(p_count) <- vertex_indexes

p_count <- numeric(length(vertex_indexes))
names(p_count) <- vertex_indexes
for(i in 1:length(pp)){
    path <- pp[[i]]
    if(length(path)>1){
      dists <- 2:length(path)
      path <- path[dists]
      p_dists[path] <- p_dists[path] + dists
      p_count[path] <- p_count + 1
    }
}

p_distances <- p_dists / p_count^2



pp <- get.all.shortest.paths(network,from=vertex_indexes[g])
sp_to_edges <- function(sp){
  l <- length(sp)
  if(l<3){
    if(l==2){
      return(sp)
    } else {
      return(c())
    }
  } else {
    paths <- c()
    for(i in 2:length(sp)){
      paths <- rbind(paths,c(sp[i-1],sp[i]))
    }
    return(paths)
  }
}
links <- do.call("rbind",lapply(pp,sp_to_edges))
named_links <- cbind(get.vertex.attribute(network,"name",index=links[,1]),get.vertex.attribute(network,"name",index=links[,2]))
subnet <- graph.data.frame(named_links)
subnet_indexes <- get_vertex_indexes(subnet)
bet <- betweenness(subnet,v=subnet_indexes)
names(bet) <- names(subnet_indexes)

