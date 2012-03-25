
source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("simulation.r")
source("evaluation.r")

home_path <- getwd()

# input params

nfams <- 3
ngenes <- 20
nsims <- 30
score_functions <- c("default_score","default_score2","default_score3")

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactome <- interactomes[["binding"]]
distance_matrix <- compute_distance_matrix(interactome)

# do simulations
family_genes_set <- simulate_family_genes_set(nfams,ngenes,nsims,distance_matrix=distance_matrix)
sim1 <- simulate_and_evaluate(family_genes_set=family_genes_set,interactome=interactome,score_function_name=score_functions[1],distance_matrix=distance_matrix)
sim2 <- simulate_and_evaluate(family_genes_set=family_genes_set,interactome=interactome,score_function_name=score_functions[2],distance_matrix=distance_matrix)
sim3 <- simulate_and_evaluate(family_genes_set=family_genes_set,interactome=interactome,score_function_name=score_functions[3],distance_matrix=distance_matrix)
sims <- list(sim1,sim2,sim3)
nscores <- length(sims)

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





