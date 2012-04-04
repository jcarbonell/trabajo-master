
source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("simulation.r")
source("evaluation.r")
source("modeling.r")

home_path <- getwd()

# input params

nfams <- 3
ngenes <- 20
nsims <- 10
score_functions <- c("default_score","default_score2","default_score3")


# prepare variables
nscores <- length(score_functions)


# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactome <- interactomes[["binding"]]
distance_matrix <- compute_distance_matrix(interactome)

# simulate family genes
all_genes <- read.delim("../misc/all_hgnc_symbols.txt",stringsAsFactors=F,header=F)$V1
family_genes_set <- simulate_iterative_family_genes(nfams,ngenes,nsims,all_genes=all_genes,distance_matrix=distance_matrix)

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



