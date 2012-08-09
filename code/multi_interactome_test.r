
source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/evaluation.r")
source("lib/postprocess.r")

home_path <- getwd()

# input params
nfams <- 3
ngenes <- 20
nsims <- 10
score_function_names <- c("default_score2","default_score3","default_score4")
global_score_methods <- c("no_zero_mean","no_zero_max")
with_known_genes <- T

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
#interactomes <- list(textmining=interactomes$textmining)


# simulate family set series
if(with_known_genes){
  all_genes <- read.table("../misc/all_hgnc_symbols.txt",header=F,stringsAsFactors=F)$V1
  family_set_series <- simulate_known_family_set_series(paste(home_path,"/../misc/omim_clean_annot.txt",sep=""),all_genes,nfams=nfams,ngenes=ngenes,nsims=nsims)
} else {
  meta_interactome  <- get_meta_interactome(interactomes)
  meta_interactome_matrix <- compute_distance_matrix(head(meta_interactome,50000))
  family_set_series <- simulate_family_set_series(nfams=nfams,ngenes=ngenes,nsims=nsims,distance_matrix=meta_interactome_matrix)  
}

# run and evaluate
sre <- run_and_evaluate(family_set_series,interactomes,score_function_names=score_function_names,global_score_methods=global_score_methods,test.inter=F)

# paint evaluation
paint_global_scores(sre$score_runs)
paint_score_density(sre$score_runs[[1]],label=score_function_names[1])
paint_score_density(sre$score_runs[[2]],label=score_function_names[2])
paint_score_density(sre$score_runs[[3]],label=score_function_names[3])

