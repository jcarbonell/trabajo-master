
source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("simulation.r")
source("evaluation.r")

home_path <- getwd()

# input params
nfams <- 3
ngenes <- 20
nsims <- 10
# score_function_names <- c("default_score","default_score2","default_score3")
score_function_names <- c("default_score2","default_score4")
global_score_methods <- c("no_zero_mean","no_zero_min","no_zero_max","sum","max")

# prepare variables
nscores <- length(score_function_names)

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))

# simulate family set series
meta_interactome  <- get_meta_interactome(interactomes)
all_genes <- unique(unique(meta_interactome$V1),unique(meta_interactome$V3))

family_set_series <- simulate_known_family_set_series(paste(home_path,"/../misc/omim_clean_annot.txt",sep=""),all_genes,nfams=nfams,ngenes=ngenes,nsims=nsims)

# run and evaluate
sre <- run_and_evaluate(family_set_series,interactomes,score_function_names=score_function_names,global_score_methods=global_score_methods)

# paint evaluation
paint_global_scores(sre$score_runs)
paint_score_density(sre$score_runs[[1]],label=score_function_names[1])
paint_score_density(sre$score_runs[[2]],label=score_function_names[2])



known_sim <- load_known_disease_genes_series("../misc/omim_clean_annot.txt",rep(3,10))

par(mfrow=c(2,1))
dms <- lapply(known_sim$fam_disease_genes,get_gene_distance_matrix,interactome=interactomes[[""]])
get_distance_vector <- function(dm){
    distance_vector <- as.numeric(dm)
    distance_vector[distance_vector>0]
}
all_distances <- unlist(lapply(dms,get_distance_vector))
hist(all_distances)




