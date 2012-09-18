
source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/evaluation.r")
source("lib/postprocess.r")
source("lib/random_walk.r")
source("lib/intermediation_utils.R")

library(igraph)

home_path <- getwd()

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))


# process interactomes
process_intermediation(interactomes$ptmod,"ptmod_intermediation2.txt")
#process_intermediation(interactomes$binding,"binding_intermediation.txt")
process_intermediation(interactomes$functional,"functional_intermediation2.txt")
#process_intermediation(interactomes$regulation,"regulation_intermediation.txt")
#process_intermediation(interactomes$textmining,"textmining_intermediation.txt")

