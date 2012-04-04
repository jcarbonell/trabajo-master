
compare_scores <- function(score_functions,nfams=NULL,ngenes=NULL,nsims=NULL,interactome=NULL,distance_matrix=NULL,paint=T,verbose=T){
    
  # prepare variables
  print(score_functions)
  nscores <- length(score_functions)

  # load interactomes
  if(is.null(distance_matrix)){
    if(verbose){
      cat("Computing distance matrix...\n")
    }
    distance_matrix <- compute_distance_matrix(interactome)
  }

  # do simulations
  family_genes_set <- simulate_family_genes_set(nfams,ngenes,nsims,distance_matrix=distance_matrix)
  sims <- list()
  for(i in 1:nscores){
    cat("\n\n\n*** SIMULATION ", i,"***\n\n")
    sims[[score_functions[i]]] <- simulate_and_evaluate(family_genes_set=family_genes_set,interactome=interactome,score_function_name=score_functions[i],distance_matrix=distance_matrix)  
  }
  
  # paint simulation trend comparison
  if(paint){
    paint_simulation_trend_comparison(sims,score_functions)
  }
  
  return(sims)
  
}

evaluate_prioritization <- function(scores_frame,fam_simulation,paint=T){
  
  fam_genes_scores <- scores_frame$score_table[fam_simulation$fam_disease_genes, "score" ]
  names(fam_genes_scores) <- fam_simulation$fam_disease_genes
  
  fam_genes_positions <- sapply(fam_simulation$fam_disease_genes, function(x) which(scores_frame$score_table$gene==x))
  fam_genes_relative_positions <- 1-(fam_genes_positions/nrow(scores_frame$score_table))
  mean_relative_position <- mean(fam_genes_relative_positions)
  
  fam_genes_positions_with_inter <- sapply(fam_simulation$fam_disease_genes,function(x) which(scores_frame$score_table_with_intermediates$gene==x))
  fam_genes_relative_positions_with_inter <- 1-(fam_genes_positions_with_inter/nrow(scores_frame$score_table_with_intermediates))
  mean_relative_position_with_inter <- mean(fam_genes_relative_positions_with_inter)
    
  if(paint){
    par(mfrow=c(2,1))
    paint_scores(scores_frame$score_table$score,fam_genes_scores)
    paint_scores(scores_frame$score_table_with_intermediates$score,fam_genes_scores,main="Score distribution (with intermediates)")
  }
  
  return(list(
    fam_genes_scores=fam_genes_scores,
    fam_genes_positions=fam_genes_positions,
    fam_genes_relative_positions=fam_genes_relative_positions,
    fam_genes_positions_with_inter=fam_genes_positions_with_inter,
    fam_genes_relative_positions_with_inter=fam_genes_relative_positions_with_inter,
    mean_relative_position=mean_relative_position,
    mean_relative_position_with_inter=mean_relative_position_with_inter
  ))
  
}

evaluate_global_prioritization <- function(global_score_table,family_set){

  nfams <- family_set$nfams
  ninteractomes <- ncol(global_score_table)

  genes_position <- matrix(rep(0,ninteractomes*nfams),nrow=nfams)

  rownames(genes_position) <- family_set$fam_disease_genes
  colnames(genes_position) <- colnames(global_score_table)
  
  genes_relative_position <- matrix(rep(0,ninteractomes*nfams),nrow=nfams)
  
  for(i in 1:ninteractomes){
    col <- colnames(global_score_table)[i]
    
    score <- as.numeric(global_score_table[,col])
    names(score) <- rownames(global_score_table)
    
    sorted_score <- sort(score,decreasing=T)
    get_pos <- function(gene){
      which(names(sorted_score)==gene)
    }
    genes_position[,i] <- sapply(family_set$fam_disease_genes,get_pos)
    
  }
  
  genes_quantile <- 1 - (genes_position / nrow(global_score_table))
  
  return(list(
      genes_position=genes_position,
      genes_quantile=genes_quantile
  ))
  
}

paint_scores <- function(all_scores,fam_genes_scores,main="Score distribution",bins=20){
  
  h <- hist(all_scores,bins,probability=T,main=main,xlab="score")
  lines(density(all_scores),col="blue")
  for(i in fam_genes_scores){
    lines(c(i,i),c(0,max(h$density*0.9)),col="red",lwd=2)  
  }
  
}

paint_simulation_trend <- function(evaluation_list,label="Relative position across simulations"){
    
  mean_relative_position <- sapply(evaluation_list,function(x) x$mean_relative_position)
  mean_relative_position_with_inter <- sapply(evaluation_list,function(x) x$mean_relative_position_with_inter)
  plot(mean_relative_position,col="red",ylim=c(0,1),pch=20,xlab="step",ylab="relative position",main=label)
  points(mean_relative_position_with_inter,col="blue",ylim=c(0,1),pch=20)
  legend(1,0.2,legend=c("input genes","input genes with intermediates"),col=c("red","blue"),cex=0.6,lwd=1)
  mean_mean_relative_position <- mean(mean_relative_position)
  sd_mean_relative_position <- sd(mean_relative_position)
  lines(c(1,nsims),rep(mean_mean_relative_position,2),col="red")
  lines(c(1,nsims),rep(mean_mean_relative_position-sd_mean_relative_position,2),col="red",lty=2)
  lines(c(1,nsims),rep(mean_mean_relative_position+sd_mean_relative_position,2),col="red",lty=2)
  mean_mean_relative_position_inter <- mean(mean_relative_position_with_inter)
  sd_mean_relative_position_inter <- sd(mean_relative_position_with_inter)
  lines(c(1,nsims),rep(mean_mean_relative_position_inter,2),col="blue")
  lines(c(1,nsims),rep(mean_mean_relative_position_inter-sd_mean_relative_position_inter,2),col="blue",lty=2)
  lines(c(1,nsims),rep(mean_mean_relative_position_inter+sd_mean_relative_position_inter,2),col="blue",lty=2)

}

paint_simulation_trend_comparison <- function(sims,sim_labels){
  
  nsteps <- sims[[1]]$sim_params$nsims
  scol <- 2
  thecols <- paste(rgb(col2rgb(seq(1,scol*length(sims),scol)),maxColorValue=255),"88",sep="")
  thefillcols <- paste(rgb(col2rgb(seq(1,scol*length(sims),scol)),maxColorValue=255),"11",sep="")
  for(i in 1:length(sims)){
    mean_relative_position <- sapply(sims[[i]]$evaluation_list,function(x) x$mean_relative_position)
    mean_mean_relative_position <- mean(mean_relative_position)
    sd_mean_relative_position <- sd(mean_relative_position)
    print(mean_mean_relative_position)
    if(i==1){
      plot(mean_relative_position,col=thecols[i],ylim=c(0,1),pch=20,xlab="step",ylab="relative position")
    } else {
      points(mean_relative_position,col=thecols[i],ylim=c(0,1),pch=20)  
    }
    lines(c(1,nsims),rep(mean_mean_relative_position,2),col=thecols[i])
  #   lines(c(1,nsims),rep(mean_mean_relative_position-sd_mean_relative_position,2),col=thecols[i],lty=2)
  #   lines(c(1,nsims),rep(mean_mean_relative_position+sd_mean_relative_position,2),col=thecols[i],lty=2)
  #   rect(1,mean_mean_relative_position-sd_mean_relative_position,nsteps,mean_mean_relative_position+sd_mean_relative_position,col=thefillcols[i],border=F)
    
  }
  legend(1,0.3,cex=0.5,lty=1,legend=sim_labels,col=thecols)

}


paint_score_groups <- function(scores_frame_list,evaluation_list,label="Scores"){
  
  fam_genes_score <- as.vector(sapply(evaluation_list,function(x) x$fam_genes_score))
  all_scores <- unlist(sapply(scores_frame_list, function(x) x$score_table$score))
  all_scores_with_inter <- unlist(sapply(scores_frame_list, function(x) x$score_table_with_inter$score))
  boxplot(list(all_with_inter=all_scores_with_inter,all=all_scores,fam=fam_genes_score),main=label,ylabel="score")
  
}

paint_mrp_comparison <- function(sims,sim_labels){
  
  get_sim_mrp <- function(sim){
    list(
      mean_relative_positions=unlist(lapply(sim$evaluation_list,function(x) x$mean_relative_position)),
      mean_relative_positions_with_inter=unlist(lapply(sim$evaluation_list,function(x) x$mean_relative_position_with_inter))
    )
  }

  sims_mrp <- lapply(sims,get_sim_mrp)
  mrps <- lapply(sims_mrp,function(x) x$mean_relative_positions)
  names(mrps) <- sim_labels
  mrps_wi <- lapply(sims_mrp,function(x) x$mean_relative_positions_with_inter)
  names(mrps_wi) <- sim_labels
  par(mfrow=c(1,2))
  boxplot(mrps,main="mean relative_position",ylim=c(0,1),cex.axis=0.7,cex.main=0.7,boxwex=0.5)
  boxplot(mrps_wi,main="mean relative_position_with_inter",ylim=c(0,1),cex.axis=0.7,cex.main=0.7,boxwex=0.5)  
  
}


