
source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("simulation.r")
source("evaluation.r")

home_path <- getwd()

# input params

nfams <- 3
ngenes <- 20
nsims <- 5
score_functions <- c("default_score","default_score2","default_score3")


# prepare variables
nscores <- length(score_functions)


# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactome <- interactomes[["binding"]]
distance_matrix <- compute_distance_matrix(interactome)

# do simulations
family_genes_set <- simulate_family_genes_set(nfams,ngenes,nsims,distance_matrix=distance_matrix)
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



# require(multicore)
# p1 <- parallel(simulate_and_evaluate(nfams=3,ngenes=20,nsims=5,interactome=interactome,score_function_name=score_functions[1],distance_matrix=distance_matrix))
# p2 <- parallel(simulate_and_evaluate(nfams=3,ngenes=20,nsims=5,interactome=interactome,score_function_name=score_functions[2],distance_matrix=distance_matrix))
# p3 <- parallel(simulate_and_evaluate(nfams=3,ngenes=20,nsims=5,interactome=interactome,score_function_name=score_functions[3],distance_matrix=distance_matrix))
# sims <- collect(list(p1,p2,p3))
# names(sims) <- score_functions
# 







