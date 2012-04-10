
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
global_score_methods <- c("no_zero_mean","no_zero_min","sum","median","min","max")

# prepare variables
nscores <- length(score_function_names)

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))

# simulate family set series
meta_interactome  <- get_meta_interactome(interactomes)
meta_interactome_matrix <- compute_distance_matrix(head(meta_interactome,50000))

family_set_series <- simulate_family_set_series(nfams=nfams,ngenes=ngenes,nsims=nsims,distance_matrix=meta_interactome_matrix)

# run and evaluate
sre <- run_and_evaluate(family_set_series,interactomes,score_function_names=score_function_names,global_score_methods=global_score_methods)

# paint evaluation
paint_global_scores(sre$score_runs)


paint_score_density(sre$score_runs[[1]],label=score_function_names[1])
paint_score_density(sre$score_runs[[2]],label=score_function_names[2])


annot_file <- "../misc/omim_clean_annot.txt"
nfams <- c(3,5,3,5)






