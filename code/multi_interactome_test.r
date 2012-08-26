
source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/evaluation.r")
source("lib/postprocess.r")
source("lib/random_walk.r")

home_path <- getwd()

# input params
nfams <- 3
ngenes <- 50
nsims <- 10
# score_function_names <- c("default_score3c")
score_function_names <- c("default_score3c","default_score3")
global_score_methods <- c("no_zero_mean","no_zero_max")
with_known_genes <- T

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
#interactomes <- list(textmining=interactomes$textmining)
# interactomes <- interactomes[c("binding","ptmod","functional")]
# interactomes <- interactomes[c("functional","ptmod")]

# simulate family set series
if(with_known_genes){
  all_genes <- read.table("../misc/all_hgnc_symbols.txt",header=F,stringsAsFactors=F)$V1
  family_set_series <- simulate_known_family_set_series(paste(home_path,"/../misc/omim_clean_annot.txt",sep=""),all_genes,nfams=nfams,ngenes=ngenes,nsims=nsims)
} else {
  meta_interactome  <- get_meta_interactome(interactomes)
  meta_interactome_matrix <- compute_distance_matrix(head(meta_interactome,50000))
  family_set_series <- simulate_family_set_series(nfams=nfams,ngenes=ngenes,nsims=nsims,distance_matrix=meta_interactome_matrix)  
}

source("lib/evaluation.r")
source("lib/score.r")
# run and evaluate
sre <- run_and_evaluate(family_set_series,interactomes,score_function_names=score_function_names,global_score_methods=global_score_methods,test.inter=F)
# paint evaluation
paint_global_scores(sre$score_runs)


source("lib/evaluation.r"); paint_global_score(sre$score_runs[[1]])

paint_score_density(sre$score_runs[[1]],label=score_function_names[1])
paint_score_density(sre$score_runs[[2]],label=score_function_names[2])
paint_score_density(sre$score_runs[[3]],label=score_function_names[3])

source("lib/random_walk.r")
system.time(
  sre_rw <- run_and_evaluate_random_walk(family_set_series,interactomes,global_score_methods=global_score_methods)
)
paint_global_scores(sre_rw$score_runs)

X11()
paint_global_score(sre_rw$score_runs[[1]],"rw")
paint_score_density(sre_rw$score_runs[[1]],label=score_function_names[1])



for( i in 1:nsims){
 a <- sre$score_runs[[1]]$multi_score_list[[i]]$scores_frame_list$functional$score_table$score
 a <- sre2$score_runs[[1]]$multi_score_list[[i]]$scores_frame_list$functional$score_table$score
 print(sum(a-b,na.rm=T))
}

nsim <- 1
global_score_table1 <- sre$score_runs[[1]]$multi_score_list[[nsim]]$global_score_table
global_score_table2 <- sre2$score_runs[[1]]$multi_score_list[[nsim]]$global_score_table
setdiff(rownames(global_score_table1),rownames(global_score_table2))

# system.time(
#   adm <- read.table("../interactomes/binding_hgnc.am",stringsAsFactor=F, header=T)
# )
# 
# 
# system.time(
# adm <- read.table("../interactomes/test.txt",
#                   header=F,sep="\t",
#                   colClasses=rep("double",10179),
#                   comment.char=""
#               )
# )
# 
# system.time(
# adm <- scan("../interactomes/test.txt",what=double(),skip=1,comment.char = "",multi.line=F)
# )

