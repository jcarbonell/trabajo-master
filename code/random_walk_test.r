
source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/evaluation.r")
source("lib/postprocess.r")
source("lib/random_walk.r")

library(igraph)

home_path <- getwd()

# input params
nfams <- 3
ngenes <- 200
nsims <- 10
global_score_methods <- c("no_zero_mean","no_zero_max")
with_known_genes <- T

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactomes <- interactomes[c("binding","ptmod","functional")]

# simulate family set series
if(with_known_genes){
  all_genes <- read.table("../misc/all_hgnc_symbols.txt",header=F,stringsAsFactors=F)$V1
  family_set_series <- simulate_known_family_set_series(paste(home_path,"/../misc/omim_clean_annot.txt",sep=""),all_genes,nfams=nfams,ngenes=ngenes,nsims=nsims)
} else {
  meta_interactome  <- get_meta_interactome(interactomes)
  meta_interactome_matrix <- compute_distance_matrix(head(meta_interactome,50000))
  family_set_series <- simulate_family_set_series(nfams=nfams,ngenes=ngenes,nsims=nsims,distance_matrix=meta_interactome_matrix)  
}

# run and evaluate
system.time(
  sre <- run_and_evaluate(family_set_series,interactomes,score_function_names="default_score3",test.inter=F,global_score_methods=global_score_methods)
)
system.time(
  sre_rw <- run_and_evaluate_random_walk(family_set_series,interactomes,global_score_methods=global_score_methods)
)

# paint evaluation
par(mfrow=c(2,1))
paint_global_score(sre_rw$score_runs[[1]],"rw")
paint_global_score(sre$score_runs[[1]],"jnf")


paint_score_density(sre_rw$score_runs[[1]],score_name="binding")




all_binding_genes <- get_interactome_genes(interactomes$binding)
binding_adm <- get_interactome_adjacency_matrix(interactomes$binding)
system.time(
rw <- RandomWalk.with.restart(all_binding_genes,"BRCA2",binding_adm)
)



weight <- function(distance){
  1-pnorm(distance*2,mean=5,sd=1.1)
}


