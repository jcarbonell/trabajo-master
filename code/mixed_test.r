
source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/score_functions.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/evaluation.r")
source("lib/postprocess.r")
source("lib/random_walk.r")

home_path <- getwd()

# input params
nfams <- 5
ndisease_genes <- 4
ngenes <- 150
nsims <- 10
score_function_names <- c("final_score")
global_score_methods <- c("no_zero_mean","no_zero_max")

# load interactomes
interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))

# simulate family set series
all_genes <- read.table("../misc/all_hgnc_symbols.txt",header=F,stringsAsFactors=F)$V1
family_set_series <- simulate_known_family_set_series(paste(home_path,"/../misc/omim_clean_annot.txt",sep=""),all_genes,nfams=nfams,ngenes=ngenes,nsims=nsims,ndisease_genes=ndisease_genes)

# load probability matrices (and synchronize with interactomes)
valid_interactomes <- c("binding","ptmod")
load("../interactomes/rwd/intermediation_matrices.rdata")
intermediation_matrices$binding <- as.matrix(intermediation_matrices$binding)
intermediation_matrices$ptmod <- as.matrix(intermediation_matrices$ptmod)
intermediation_matrices <- intermediation_matrices[valid_interactomes]
load("../interactomes/rwd/all_rwd.rdata")
probability_matrices <- probability_matrices[valid_interactomes]
interactomes <- interactomes[valid_interactomes]
gc()

# shortest path
source("lib/score.r")
source("lib/score_functions.r")
sre_nop <- run_and_evaluate(family_set_series,interactomes,score_function_names=score_function_names,global_score_methods=global_score_methods,test.inter=F)
paint_global_score(sre_nop$score_runs[[1]],label="Shortest path")
paint_score_density(sre_nop$score_runs[[1]],score_name="binding")

# intermediation
source("lib/score.r")
sre_inter <- run_and_evaluate(family_set_series,interactomes,probability_matrices=intermediation_matrices,score_function_names=score_function_names,global_score_methods=global_score_methods,test.inter=F)
paint_global_scores(sre_inter$score_runs)
paint_score_density(sre_inter$score_runs[[1]],score_name="binding")

# with random walk probability
source("lib/score.r")
sre <- run_and_evaluate(family_set_series,interactomes,probability_matrices=probability_matrices,score_function_names=score_function_names,global_score_methods=global_score_methods,test.inter=F)
paint_global_scores(sre$score_runs)

# pure random walk
sre_rw <- run_and_evaluate_random_walk(family_set_series,interactomes,global_score_methods=global_score_methods)
paint_global_scores(sre_rw$score_runs)

# fine analysis
paint_global_score(sre_nop$score_runs[[1]])
paint_global_score(sre_nop$score_runs[[2]])


paint_global_score(sre$score_runs[[1]])
paint_global_score(sre$score_runs[[2]])
paint_global_score(sre_rw$score_runs[[1]])

paint_score_density(sre$score_runs[[1]],label=score_function_names[1],score_name="binding")
paint_score_density(sre_nop$score_runs[[1]],label=score_function_names[1])




adm <- get_interactome_adjacency_matrix(interactomes$binding)

fam_index <- 5
fam_genes <- family_set_series[[fam_index]]$raw_gene_list[[1]]
rest_of_genes <- unique(c(family_set_series[[fam_index]]$raw_gene_list[[2]],family_set_series[[fam_index]]$raw_gene_list[[3]]))
rw <- RandomWalk.with.restart(fam_genes,rest_of_genes,adm)

hist(rw$final.prob,50)
rw$final.prob[fam_genes[1]]

nas <- sum(is.na(rw$final.prob))
nas
no_nas <- length(rw$final.prob)-nas
no_nas 
sum(rw$final.prob<=rw$final.prob[fam_genes[1]],na.rm=T)/no_nas

rest_quantiles <- c()
for(i in 1:nsims){

  for(j in 1:nfams){
    
    fam_gene <- family_set_series[[i]]$fam_disease_genes[j]
    rest_of_genes <- family_set_series[[i]]$fam_disease_genes[-j]
    
    if(fam_gene %in% rownames( probability_matrices$binding)){
      node_dis <- probability_matrices$binding[fam_gene,]
      rest_dis <- node_dis[rest_of_genes]
      rest_quantiles <- c(rest_quantiles, sapply(rest_dis, function(x) sum(node_dis<=x)/length(node_dis)))
      
    }
  }

}

boxplot(rest_quantiles,outline=F,ylim=c(0,1))
points(rep(1,length(rest_quantiles))+runif(length(rest_quantiles),min=-0.03,max=0.03),rest_quantiles,pch=20,cex=0.4)



family_set <- family_set_series[[9]]
disease_genes <- family_set$fam_disease_gene
all_raw_genes <- unique(unlist(family_set$raw_gene_list))
super_list <- all_raw_genes
super_list_fams <- list()
get_gene_fams <- function(gene){
  gene_fams <- names(which(lapply(family_set$raw_gene_list,function(fam) is.element(gene,fam))==T))
  if(length(gene_fams)==0){
    gene_fams <- inter_family
  }
  return(gene_fams)
}
super_list_fams <- lapply(super_list,get_gene_fams)
names(super_list_fams) <- super_list


yeah <- function(node){
  default_score3c(node,probability_matrices$binding,weights=rep(1,10),distance_priors=rep(1,10),super_list_fams=super_list_fams)
}
yeah2 <- function(node){
  default_score3d(node,probability_matrices$binding,weights=rep(1,10),distance_priors=rep(1,10),super_list_fams=super_list_fams)
}

subnet <- get.all.shortest.paths.Josete(interactome,all_raw_genes,5)
subnet.igraph <- get.igraph(subnet)
subnet.wc <- walktrap.community(subnet.igraph)
names(subnet.wc$modularity) <- subnet.wc$labels
names(subnet.wc$membership) <- subnet.wc$labels

cluster.score <- function(node,scores,membership){  
  group <- membership[node]
  if(is.na(group)){
    return(NA)
  } else {    
    neighbours <- setdiff(names(membership)[which(membership==group)],node)
    neighbours <- neighbours[neighbours %in% names(scores)]
    if(length(neighbours)==0){
      return(scores[node])
    } else {
      neighbours_scores <- scores[neighbours]
      return(mean(neighbours_scores[which(neighbours_scores>=quantile(neighbours_scores,0.75,na.rm=T))],na.rm=T))
    }
  }  
}
cluster.score(disease_genes[3],scoresc,subnet.wc$membership)

source("lib/score.r")
scoresc <- c()
scoresd <- c()
for(node in all_raw_genes){
    aa <- yeah(node)
    scoresc <- c(scoresc,aa$score)
    bb <- yeah2(node)
    scoresd <- c(scoresd,bb$score)
}
names(scoresc) <- all_raw_genes
names(scoresd) <- all_raw_genes

# scoresc <- scoresc/ip[all_raw_genes,"burt"]
# scoresc <- scoresc*ip[all_raw_genes,"degree"]
# scoresc <- scoresc/subnet.wc$modularity[all_raw_genes]
scoresc_cs <- sapply(all_raw_genes,cluster.score,scores=scoresc,membership=wc$membership)
scoresc <- scoresc*0.5 + scoresc_cs*0.5
# scoresc <- scoresc/scoresc_cs
# scoresd <- scoresd*ip[all_raw_genes,"degree"]



sorted_scoresc <- sort(scoresc,decreasing=T,na.last=T)
sorted_scoresd <- sort(scoresd,decreasing=T,na.last=T)
thec <- sapply(disease_genes,function(x) sum(sorted_scoresc<=sorted_scoresc[x],na.rm=T)/sum(!is.na(sorted_scoresc)))
thed <- sapply(disease_genes,function(x) sum(sorted_scoresd<=sorted_scoresd[x],na.rm=T)/sum(!is.na(sorted_scoresd)))
disease_genes
rbind(thec,thed)




