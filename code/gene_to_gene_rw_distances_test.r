
source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/evaluation.r")
source("lib/postprocess.r")
source("lib/random_walk.r")

library(igraph)

home_path <- getwd()

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))

# ptmod
process <- function(interactome, name){
  cat(name,">>>>>>>>>>>>>>>>>>>>>>>>>\n")
  adm <- get_interactome_adjacency_matrix(interactome)
  get_gene_to_gene_random_walk_distance(adm,paste(name,"_g2g_rwd.txt",sep=""))
  rm(adm)
  gc()  
}

# process interactomes
process(interactomes$ptmod,"ptmod")
process(interactomes$binding,"binding")
process(interactomes$functional,"functional")
process(interactomes$regulation,"regulation")
process(interactomes$textmining,"textmining")

