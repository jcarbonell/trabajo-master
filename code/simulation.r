
simulate_family_genes_set <- function(nfams,ngenes,nsims,distance_matrix=NULL,interactome=NULL,fam_name_prefix="fam_") {
  family_genes_set <- list()
  for(k in 1:nsims){
    family_genes_set[[k]] <- simulate_family_genes(nfams,ngenes,distance_matrix=distance_matrix,interactome=interactome)
  }
  return(family_genes_set)
}

simulate_family_genes <- function(nfams,ngenes,distance_matrix=NULL,interactome=NULL,fam_name_prefix="fam_") {

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
    random_list <- sample(interactors,ngenes)
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


simulate_and_evaluate <- function(family_genes_set=NULL,nfams=3,ngenes=20,nsims=5,interactome,distance_matrix=NULL,score_function_name="default_score",verbose=T){
  
  start <- proc.time()
  
  message <- function(...){
    if(verbose){
      cat(">>",...,"\n",sep="")
    }
  }
  
  # nsims
  if(!is.null(family_genes_set)){
    nsims <- length(family_genes_set)
  }
  
  score_function <- get(score_function_name) 
  
  if(verbose){
    cat("Simulation",date(),"\n")
    message("score function: ",score_function_name)
    message("params: steps=",nsims,", nfams=",nfams,",ngenes=",ngenes)
  }
    
  # prepare distance matrix
  if(is.null(distance_matrix)){
    if(verbose){cat(">>constructing distance_matrix..\n")}
    distance_matrix <- compute_distance_matrix(interactome)
  }
  
  # simulation
  if(is.null(family_genes_set)){
    message("simulating ",nsims," families")
    family_genes_set <- simulate_family_genes_set(nfams,ngenes,nsims,distance_matrix=distance_matrix)
  }

  # init working lists
  scores_frame_list <- list()
  evaluation_list <- list()
  
  for(k in 1:nsims){
  
    message("simulation step",k)
    
    # score computation
    #cat("  computing score...\n")
    scores_frame_list[[k]] <- compute_score(family_genes_set[[k]]$raw_gene_list,interactome,score_function)
    
    # evaluation
    #cat("  evaluating prioritization...\n")
    evaluation_list[[k]] <- evaluate_prioritization(scores_frame_list[[k]],family_genes_set[[k]],paint=F)
        
  }

  end <- proc.time() - start
  
  message("finished in ",end["elapsed"]," seconds\n")
  
  return(list(
      score_function_name=score_function_name,
      sim_params=list(        
        family_genes_set=family_genes_set,
        nfams=nfams,
        ngenes=ngenes,
        nsims=length(family_genes_set)
      ),
      scores_frame_list=scores_frame_list,
      evaluation_list=evaluation_list,
      info=list(
        date=Sys.Date(),
        time=end["elapsed"]
      )
  ))
  
}
