
source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/evaluation.r")
#source("lib/modeling.r")
source("lib/postprocess.r")

home_path <- getwd()

# input params

nfams <- 3
ngenes <- 20
nsims <- 10
score_function_names <- c("default_score","default_score4")
global_score_methods <- c("no_zero_mean","no_zero_max")

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactome <- interactomes[["binding"]]
distance_matrix <- compute_distance_matrix(interactome)
interactomes <- list(binding=interactome)

# simulate family genes
all_genes <- read.delim("../misc/all_hgnc_symbols.txt",stringsAsFactors=F,header=F)$V1
family_set_series <- simulate_iterative_family_set_series(nfams,ngenes,nsims,all_genes=all_genes,distance_matrix=distance_matrix)

# run and evaluate
sre <- run_and_evaluate(family_set_series,interactomes,score_function_names=score_function_names,global_score_methods=global_score_methods)

# paint evaluation
paint_global_scores(sre$score_runs)
paint_score_density(sre$score_runs[[1]],label=score_function_names[1])
paint_score_density(sre$score_runs[[2]],label=score_function_names[2])

# run and evaluate
sims <- list()
for(i in 1:nscores){
  cat("\n\n\n*** SIMULATION ", i,"***\n\n")
  sims[[score_functions[i]]] <- run_and_evaluate(family_genes_set=family_genes_set,interactome=interactome,score_function_name=score_functions[i],distance_matrix=distance_matrix)  
}

# paint simulation trend
par(mfrow=c(1,nscores))
for(k in 1:nscores){
  paint_simulation_trend(sims[[k]]$evaluation_list,label=score_functions[k])  
}

# paint simulation trend comparison
paint_simulation_trend_comparison(sims,score_functions)

# score distribution
par(mfrow=c(1,nscores))
for(k in 1:nscores){
  paint_score_groups(sims[[k]]$scores_frame_list,sims[[k]]$evaluation_list,label=score_functions[k])  
}

# mean relative position
paint_mrp_comparison(sims,score_functions)



