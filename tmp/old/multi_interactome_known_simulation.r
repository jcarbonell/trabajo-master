
source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("simulation.r")
source("evaluation.r")
source("postprocess.r")

home_path <- getwd()

# input params
nfams <- 3
ngenes <- 50
nsims <- 10

score_function_names <- c("default_score1b","default_score4")
global_score_methods <- c("no_zero_mean","no_zero_max")

# prepare variables
nscores <- length(score_function_names)

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactomes$textmining <- NULL

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
paint_score_density(sre$score_runs[[3]],label=score_function_names[3])
paint_score_density(sre$score_runs[[4]],label=score_function_names[4])

paint_score_density(sre$score_runs[[1]],label=score_function_names[1],"binding")
paint_score_density(sre$score_runs[[2]],label=score_function_names[2],"binding")
paint_score_density(sre$score_runs[[3]],label=score_function_names[3],"gene_cluster_score")
paint_score_density(sre$score_runs[[4]],label=score_function_names[4],"gene_cluster_score")

###############
###############
###############

family_set <- family_set_series[[1]]
interactome <- interactomes[["binding"]]
radio <- 2

all_raw_genes <- unique(unlist(family_set$raw_gene_list))

subnet <- get_subnet(all_raw_genes,interactomes[["binding"]],1)
dm_with_inter <- compute_distance_matrix(subnet)
dm <- subnet_dm_with_inter[all_raw_genes,all_raw_genes]
# subnet_summary <- describe_interactome(subnet)

get_neighbourhood <- function(row){
  names(which(row>0 & row<=1))
}

  
  
###############
###############
###############
global_score_table <- sre$score_runs$default_score$multi_score_list[[1]]$global_score_table
raw_genes <- rownames(global_score_table)
cl <- compute_clusters(raw_genes,interactomes)
genes <- names(cl$clusters)[which(cl$clusters>-1)]
net <- makeNetwork(cl$pseudo_interactome$V1,cl$pseudo_interactome$V3)
subbionet <- subNetwork(genes,net)
subbionet <- rmSelfLoops(subbionet)
plot3dModule(subbionet)
pval <- global_score_table[genes,1]
names(pval) <- genes
fb <- fitBumModel(pval, plot = T)
fit_scores <- scoreNodes(subbionet, fb, fdr = 0.001)

module <- runFastHeinz(subbionet,pval)
plotModule(module,diff.expr=pval)


  
###############
###############
###############

global_score_table <- sre$score_runs$default_score$multi_score_list[[1]]$global_score_table



which(names(single_scores) %in% family_set$fam_disease_genes)
which(names(group_scores_sorted) %in% family_set$fam_disease_genes)

###############
###############
###############
'known_sim <- load_known_disease_genes_series("../misc/omim_clean_annot.txt",rep(5,100))

par(mfrow=c(2,1))
dms <- lapply(known_sim$fam_disease_genes,get_gene_distance_matrix,interactome=meta_interactome)
get_distance_vector <- function(dm){
    distance_vector <- as.numeric(dm)
    distance_vector[distance_vector>0]
}
all_distances <- unlist(lapply(dms,get_distance_vector))
hist(all_distances)


interactomes <- list()
interactomes[["meta"]] <- meta_interactome

