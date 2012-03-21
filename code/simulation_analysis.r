
source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("simulation.r")

home_path <- getwd()

nfams <- 3
ngenes <- 20

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactome <- interactomes[["binding"]]
distance_matrix <- compute_distance_matrix(interactome)

# simulate family genes
fam_simulation <- simulate_family_genes(nfams,ngenes,distance_matrix=distance_matrix)

# score computation
scores_frame <- compute_score(fam_simulation$raw_gene_list,interactome,default_score2)

# evaluation
evaluation <- evaluate_prioritization(scores_frame,fam_simulation)



