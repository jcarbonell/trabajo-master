
simulate_family_set_series <- function(nfams=NULL,ngenes=NULL,nsims=NULL,distance_matrix=NULL,interactome=NULL,fam_name_prefix="fam_",verbose=T) {
  
  # default
  if(is.null(nsims)){
    nsims <- 5
  }
  if(verbose){
    cat("Simulating",nsims,"simulations...\n")
  }
  
  # simulate set
  family_genes_set <- list()
  for(k in 1:nsims){
    family_genes_set[[k]] <- simulate_family_set(nfams,ngenes,distance_matrix=distance_matrix,interactome=interactome)
  }
  return(family_genes_set)
}

simulate_iterative_family_set_series <- function(nfams=NULL,ngenes=NULL,nsims=NULL,all_genes=NULL,distance_matrix=NULL,interactome=NULL,fam_name_prefix="fam_",verbose=T) {
  
  # default
  if(is.null(nsims)){
    nsims <- 5
  }
  if(verbose){
    cat("Simulating",nsims,"simulations...\n")
  }
  
  # prepare init family genes
  init_family_genes <- simulate_family_set(nfams,ngenes,distance_matrix=distance_matrix)
  
  # init inter genes
  init_genes <- unique(unlist(init_family_genes$raw_gene_list))
  prepared_genes <- setdiff(all_genes,init_genes)
  
  # simulate
  family_genes_set <- list()
  for(i in 1:nsims){
    inters <- sample(prepared_genes,nfams)
    new_fam <- init_family_genes
    for(k in 1:nfams){
      new_fam$raw_gene_list[[k]] <- c(new_fam$raw_gene_list[[k]],inters[k])
      new_fam$random_gene_list[[k]] <- c(new_fam$random_gene_list[[k]],inters[k])
    }
    new_fam$ngenes <- new_fam$ngenes + 1
    family_genes_set[[i]] <- new_fam
  }
  
  return(family_genes_set)
}

simulate_family_set <- function(nfams=NULL,ngenes=NULL,distance_matrix=NULL,interactome=NULL,fam_name_prefix="fam_") {

  # defaults
  if(is.null(nfams)){
    nfams <- 3  
  }
  if(is.null(ngenes)){
    ngenes <- 20
  }
  
  # family names
  fam_names <- paste(fam_name_prefix,c(1:nfams),sep="")
  
  # select a disease gene per family
  if(is.null(distance_matrix)){
#     cat("computing distance matrix...\n")
    distance_matrix <- compute_distance_matrix(interactome)
  } else {
#     cat("distance matrix provided...\n")
  }
  
  correct <- F
  while(correct==F){
    row <- as.integer(runif(1,min=1,max=nrow(distance_matrix)))
    pivot <- colnames(distance_matrix)[row]
    neighborhood <- names(which(distance_matrix[row,]<3 & distance_matrix[row,]>0))
    if(length(neighborhood)>nfams){
      fam_disease_genes <- c(pivot,sample(neighborhood,nfams-1))
      names(fam_disease_genes) <- fam_names
      correct <- T
    } 
  }
    
  # add other random genes
  interactors <- rownames(distance_matrix)
  raw_gene_list <- list()
  random_gene_list <- list()
  for(i in 1:nfams){
    random_list <- sample(interactors,ngenes-1)
    raw_gene_list[[fam_names[i]]] <- c(as.character(fam_disease_genes[i]),random_list)  
    random_gene_list[[fam_names[i]]] <- random_list 
  }

  return(list(
    raw_gene_list=raw_gene_list,
    random_gene_list=random_gene_list,
    pivot=pivot,
    fam_disease_genes=fam_disease_genes,
    nfams=nfams,
    ngenes=ngenes
  ))
  
}


simulate_run_and_evaluate <- function(nfams=3,ngenes=20,nsims=5,interactomes,distance_matrix=NULL,score_function_name="default_score",verbose=T,global_score_method=NULL){
  
  # prepare distance matrix
  if(is.null(distance_matrix)){
    if(verbose){cat(">>constructing distance_matrix for family simulation..\n")}
    distance_matrix <- compute_distance_matrix(interactomes[["binding"]])
  }
  
  # simulation
  cat("simulating",nfams,"families",nsims,"times\n")
  family_genes_set <- simulate_family_set_series(nfams,ngenes,nsims,distance_matrix=distance_matrix)
  
  sim <- run_and_evaluate(family_genes_set,interactomes,distance_matrix,score_function_name=score_function_name,verbose=verbose)
  
}

run_and_evaluate <- function(family_set_series,interactomes,distance_matrix=NULL,score_function_name="default_score",verbose=T,global_score_method=NULL){
  
  start <- proc.time()
  
  message <- function(...){
    if(verbose){
      cat(">>",...,"\n",sep="")
    }
  }
  
  # ninteractomes
  ninteractomes <- length(interactomes)
  
  # nsims
  nsims <- length(family_set_series)
  
  # init function
  score_function <- get(score_function_name)
  
  # family set simulation params
  nfams_vector <- sapply(family_set_series,function(x) x$nfams)
  ngenes_vector <- sapply(family_set_series,function(x) x$ngenes)
  nfams_range <- range(nfams_vector)
  ngenes_range <- range(ngenes_vector)
  
  if(verbose){
    cat("Run and evaluate",date(),"\n")
    message("score function: ",score_function_name)
    message("params: steps=",nsims,", nfams=[",paste(nfams_range,collapse=","),"],ngenes=[",paste(ngenes_range,collapse=","),"]")
  }
    
  # init working lists
  multi_score_list <- list()
  global_evaluation_list <- list()
  global_evaluation_with_inter_list <- list()
  
#   evaluation_list <- list()
  
  for(k in 1:nsims){
  
    message("running step",k)
  
    # run
    multi_score_list[[k]] <- compute_multi_score(family_set_series[[k]]$raw_gene_list,interactomes,score_function,global_score_method=global_score_method)
    
    # evaluate
    global_evaluation_list[[k]] <- evaluate_global_prioritization(multi_score_list[[k]]$global_score_table,family_set_series[[k]])
    global_evaluation_with_inter_list[[k]] <- evaluate_global_prioritization(multi_score_list[[k]]$global_score_table_with_inter,family_set_series[[k]])
  }

  end <- proc.time() - start
  
  message("finished in ",end["elapsed"]," seconds\n")
  
  return(list(
      score_function_name=score_function_name,
      sim_params=list(        
        family_set_series=family_set_series,
        nfams_vector=nfams_vector,
        ngenes_vector=ngenes_vector,
        nfams_range=nfams_range,
        ngenes_range=ngenes_range,
        nsims=length(family_set_series)
      ),
      multi_score_list=multi_score_list,
      global_evaluation_list=global_evaluation_list,
      global_evaluation_with_inter_list=global_evaluation_with_inter_list,
      info=list(
        date=Sys.Date(),
        time=end["elapsed"]
      )
  ))
  
}


