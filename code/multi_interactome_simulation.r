
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

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
meta_interactome  <- get_meta_interactome(interactomes)

# precompute distance matrices
meta_interactome_matrix <- compute_distance_matrix(head(meta_interactome,50000))

# simulate family genes
family_set_series <- simulate_family_set_series(nfams=nfams,ngenes=ngenes,nsims=nsims,distance_matrix=meta_interactome_matrix)

# run and evaluate
sre <- run_and_evaluate(family_set_series,interactomes,score_function_name="default_score2",global_score_method=c("no_zero_mean","sum","median","min","max"))

# paint evaluation
all_quantiles <- do.call("rbind",lapply(sre$global_evaluation_list,function(x) x$genes_quantile))

boxplot(all_quantiles)