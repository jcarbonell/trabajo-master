
get.column.normalized.adjency.matrix <- function(my.igraph){
  igraph.adjM <-  get.adjacency(my.igraph)
  column.normalized.igraph.adjM <- matrix(NA, ncol=dim(igraph.adjM)[2], nrow=dim(igraph.adjM)[1])
  for (i in 1:dim(igraph.adjM)[1]){
    the_sum <- sum(igraph.adjM[,i])
    #if(the_sum>0){
      column.normalized.igraph.adjM[,i] <- igraph.adjM[,i] / the_sum
    #}
  }
  colnames(column.normalized.igraph.adjM) <- colnames(igraph.adjM)
  rownames(column.normalized.igraph.adjM) <- rownames(igraph.adjM)
  return(column.normalized.igraph.adjM)
}


get_interactome_adjacency_matrix <- function(interactome){
  
  # clean interactome
  my.graph <- interactome[,c(1,3)]
  my.graph <- my.graph[which(my.graph[,1]!=my.graph[,2]),]
  my.igraph <- simplify(graph.data.frame(my.graph, directed=FALSE))
  
  # get adjacency matrix
  igraph.adjM <-  get.adjacency(my.igraph)
  
  # normalize by column
  column.normalized.igraph.adjM <- matrix(NA, ncol=dim(igraph.adjM)[2], nrow=dim(igraph.adjM)[1])
  for (i in 1:dim(igraph.adjM)[1]){
    the_sum <- sum(igraph.adjM[,i])
    if(the_sum>0){
    column.normalized.igraph.adjM[,i] <- igraph.adjM[,i] / the_sum
    }
  }
  colnames(column.normalized.igraph.adjM) <- colnames(igraph.adjM)
  rownames(column.normalized.igraph.adjM) <- rownames(igraph.adjM)
  
  return(column.normalized.igraph.adjM)
  
}


RandomWalk.with.restart <- function(final.nodes.name, start.nodes.name, igraph.norm.adjM, prob.restart=0.3,thresh=1e-08,max.iter=50){
  
  all.nodes.names <- colnames(igraph.norm.adjM)
  prob.t.0 <- rep(0, length(all.nodes.names))
  prob.t.0[ which(all.nodes.names %in% start.nodes.name) ] <- 1 / length(start.nodes.name)
  prob.t <- prob.t.0
  a <- 0; b <- 1; count <- 1
  
  change <- c()
  
  last_change <- max(abs(b-a), na.rm=T)
  
  while ( last_change > thresh &&  count < max.iter ){
    
    a <- prob.t    
    prob.t <- (1 - prob.restart) * ( igraph.norm.adjM %*% prob.t ) + prob.restart * prob.t.0
    b <- prob.t
    
    count <- count + 1
    last_change <- max(abs(b-a), na.rm=T)
    change <- c(change,last_change)
    
  }
  
  all.nodes.names_pos <- 1:length(all.nodes.names)
  names(all.nodes.names_pos) <- all.nodes.names
  final.nodes.name_pos <- all.nodes.names_pos[final.nodes.name]   
    
  final.prob <- as.numeric(prob.t[final.nodes.name_pos])
  names(final.prob) <- final.nodes.name
  
#   final.prob[is.na(final.prob)] <- 0
  
  #final.prob <- as.numeric(prob.t[which(all.nodes.names %in% final.nodes.name)])
  
  return(list(
      final.prob = final.prob,
      count = count,
      change = change,
      thresh = thresh,
      max.iter = max.iter
    )
      
  )
         
  
}


compute_random_walk_score <- function(raw_gene_list,interactome,inter_family="INTER",verbose=T){
    
  ad_matrix <- get_interactome_adjacency_matrix(interactome)  
  
  #net <- get.igraph(interactome)
  #ad_matrix <- read.delim(gzfile(ad_matrix_file),sep="\t",header=T,row.names=T)
  
  # profile families
  all_raw_genes <- unique(unlist(raw_gene_list)) 
  query_genes <- unique(all_raw_genes)
  all_genes <- all_raw_genes #intersect(all_raw_genes,rownames(ad_matrix))
  uninteractors <- setdiff(all_raw_genes,rownames(ad_matrix))
  fam_hits <- compute_fam_hits(raw_gene_list)
  
  scores <- numeric(length(all_genes))
  names(scores) <- all_genes
  
  for(i in 1:length(raw_gene_list)) {
    
    fam_genes <- raw_gene_list[[i]]
    rest_of_genes <- setdiff(all_genes,fam_genes)
    
    rw <- RandomWalk.with.restart(fam_genes,rest_of_genes,ad_matrix,prob.restart=0.0,thresh=1e-05,max.iter=10)    
    fam_scores <- rw$final.prob
    
    scores[fam_genes] <- scores[fam_genes] + fam_scores
    
  }
  
  # free mem
  rm(ad_matrix)
  gc()
  
   # normalization
  node_scores <- (scores-min(scores,na.rm=T))/(max(scores,na.rm=T)-min(scores,na.rm=T))
#   group_scores <- (group_scores-min(group_scores))/(max(group_scores)-min(group_scores)) 
  
  # fam_hits  
  get_fam_hit <- function(fams){return(paste(fams,collapse=","))}
  collapsed_fam_hits <- sapply(fam_hits,get_fam_hit)
  
  # page ranke weight
  my.graph <- interactome[,c(1,3)]
  my.graph <- my.graph[which(my.graph[,1]!=my.graph[,2]),]
  my.igraph <- simplify(graph.data.frame(my.graph, directed=FALSE))
  pr <- page.rank(my.igraph)
  page_ranks <- pr$vector
  names(page_ranks) <- get.vertex.attribute(my.igraph,"name")
  node_page_ranks <- page_ranks[query_genes]
  names(node_page_ranks) <- query_genes
  node_page_ranks[is.na(node_page_ranks)] <- min(node_page_ranks,na.rm=T)
  node_page_ranks <- (node_page_ranks-min(node_page_ranks))/(max(node_page_ranks)-min(node_page_ranks))
#   node_page_ranks <- node_page_ranks +1 
  
#   node_scores <- node_scores*(node_page_ranks)
#   node_scores <- node_page_ranks
  
  # sorting index
  sorting_index <- rev(order(node_scores))
  
  # result table
  result_table <- data.frame(
    gene=names(node_scores[sorting_index]),
    score=node_scores[sorting_index],
#     group_score=group_scores[sorting_index],
    fams=collapsed_fam_hits[sorting_index],
    stringsAsFactors=F
  )
  
  # return
  scores <- list(
    
    # score
    score_table=result_table,
    
    fam_hits=fam_hits,
        
    # uninteractors
    uninteractors=uninteractors,
    uninteractor_ratio=length(uninteractors)/length(all_raw_genes)
    
  )
  
  return(scores)
  
  
}


compute_random_walk_multi_score <- function(raw_gene_list,interactomes,inter_family="INTER",verbose=T,global_score_methods=NULL){
  
  # run and evaluate
  scores_frames <- list()
  
  for(i in 1:length(interactomes)){
    
    if(verbose==T){
      cat("Computing score with",names(interactomes)[i],"interactome\n")
    }
    
    scores_frames[[ names(interactomes)[i] ]] <- compute_random_walk_score(raw_gene_list,interactomes[[i]],verbose=verbose)
    
  }
  
  global_score_table <- compute_global_score_table(scores_frames,inter=F,method=global_score_methods)  
  global_score_frame <- as.data.frame(global_score_table)
  global_score_frame$fams <- get_fams_from_score_frame(scores_frames[[1]],rownames(global_score_table))
 
  return(list(
    scores_frame_list = scores_frames,
    global_score_table = global_score_table,
    global_score_frame = global_score_frame
  ))
  
}


compute_fam_hits <- function(raw_gene_list,genes=NULL){
  
  # gene list
  if(is.null(genes)){
    genes <- all_raw_genes <- unique(unlist(raw_gene_list))  
  }
  
  # compute fam hits
  family_names <- names(raw_gene_list)
  nfams <- length(unique(family_names))  
  
  # obtain fams by gene
  super_list_fams <- list()
  get_gene_fams <- function(gene){
    gene_fams <- names(which(lapply(raw_gene_list,function(fam) is.element(gene,fam))==T))
    if(length(gene_fams)==0){
      gene_fams <- inter_family
    }
    return(gene_fams)
  }
  super_list_fams <- lapply(genes,get_gene_fams)
  names(super_list_fams) <- genes

  return(super_list_fams)
  
}

get_gene_to_gene_random_walk_distance <- function(ad_matrix,outfile,genes=NULL){
  
  if(is.null(genes)){
    genes <- rownames(ad_matrix)
  }  
  ngenes <- length(genes)
    
  write(paste(genes,collapse="\t"),file=outfile)
  
  total <- 0
  step <- ngenes/100
  percent <- 0
  cont <- 0
  
  cat("Processed 0 percent...\n")
  
  for(i in 1:ngenes){
    
    rw <- RandomWalk.with.restart(genes,genes[i],ad_matrix)
    write(paste(genes[i],paste(format(rw$final.prob,digits=4,scientific=F),collapse="\t"),sep="\t"),file=outfile,append=T)
    
    cont = cont+1
    if(cont>=step){
      total <- total + cont
      percent <- format((total/ngenes)*100,digits=2)        
      cat("Processed",percent,"percent...\n")
      cont <- 0
    }
    
  }
  
  if(percent<100){
    cat("Processed 100 percent...\n")
  }
      
}
compute_gene_to_gene_random_walk_distance <- function(ad_matrix,genes=NULL,verbose=F){
  
  if(is.null(genes)){
    genes <- rownames(ad_matrix)
  }  
  ngenes <- length(genes)
  
  distance_matrix <- matrix(0,ncol=ngenes,nrow=ngenes)
  colnames(distance_matrix) <- genes
  rownames(distance_matrix) <- genes
  
  total <- 0
  step <- ngenes/100
  percent <- 0
  cont <- 0
  
  if(verbose) {cat("Processed 0 percent...\n")}
  
  for(i in 1:ngenes){
    
    rw <- RandomWalk.with.restart(genes,genes[i],ad_matrix)
    distance_matrix[i,] <- rw$final.prob
  
    cont = cont+1
    if(cont>=step){
      total <- total + cont
      percent <- format((total/ngenes)*100,digits=2)        
      if(verbose) {cat("Processed",percent,"percent...\n")}
      cont <- 0
    }
    
  }

  if(verbose) {
    if(percent<100){
      cat("Processed 100 percent...\n")
    }
  }

  distance_matrix
}


load_gene_to_gene_random_walk_distance <- function(distance_file){  
  as.matrix.data.frame(read.delim(distance_file,sep="\t",stringsAsFactors=F,header=T))
}


