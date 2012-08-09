
source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/evaluation.r")
source("lib/postprocess.r")

home_path <- getwd()

# input params
nfams <- 3
ngenes <- 50
nsims <- 10

score_function_names <- c("default_score","default_score4")
global_score_methods <- c("no_zero_mean","no_zero_max")

# prepare variables
nscores <- length(score_function_names)

# load interactomes
raw_interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactomes <- list(binding=raw_interactomes$binding)

# simulate family set series
meta_interactome  <- get_meta_interactome(interactomes)

# known genes simulation
# all_genes <- unique(unique(meta_interactome$V1),unique(meta_interactome$V3))
# family_set_series <- simulate_known_family_set_series(paste(home_path,"/../misc/omim_clean_annot.txt",sep=""),all_genes,nfams=nfams,ngenes=ngenes,nsims=nsims)

# pure simulated genes
meta_interactome_matrix <- compute_distance_matrix(meta_interactome[sample(1:nrow(meta_interactome),50000),])
family_set_series <- simulate_family_set_series(nfams=nfams,ngenes=ngenes,nsims=nsims,distance_matrix=meta_interactome_matrix)

# run and evaluate
sre <- run_and_evaluate(family_set_series,interactomes,score_function_names=score_function_names,global_score_methods=global_score_methods)

# paint evaluation
paint_global_scores(sre$score_runs)
paint_score_density(sre$score_runs[[1]],label=score_function_names[1])
paint_score_density(sre$score_runs[[2]],label=score_function_names[2])

