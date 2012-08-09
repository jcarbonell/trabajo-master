
compare_scores <- function(score_functions,nfams=NULL,ngenes=NULL,nsims=NULL,interactome=NULL,distance_matrix=NULL,paint=T,verbose=T){
    
  # prepare variables
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

  # target genes quantiles
  genes_position <- matrix(rep(0,ninteractomes*nfams),nrow=nfams)
  rownames(genes_position) <- family_set$fam_disease_genes
  colnames(genes_position) <- colnames(global_score_table)
      
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
  
  # split global_score_table by (INTER,ND_FAM_GENES,DISEASE)
  genes <- rownames(global_score_table)
  gene_fams <- sapply(genes,get_gene_fams,family_set=family_set)
  target_indexes <- which(genes %in% family_set$fam_disease_genes)
  inter_indexes <- which(gene_fams =="INTER")
  random_indexes <- setdiff(1:length(genes),target_indexes)
  random_indexes <- setdiff(random_indexes,inter_indexes)
  target_global_score_table <- global_score_table[target_indexes,]
  inter_global_score_table <- global_score_table[inter_indexes,]
  random_global_score_table <- global_score_table[random_indexes,]
  group_freqs <- c(nrow(inter_global_score_table),nrow(random_global_score_table),nrow(target_global_score_table))
  names(group_freqs) <- c("INTER","RANDOM","TARGET")
  
  return(list(
      genes_position=genes_position,
      genes_quantile=genes_quantile,
      target_global_score_table=target_global_score_table,
      inter_global_score_table=inter_global_score_table,
      random_global_score_table=random_global_score_table,
      group_freqs=group_freqs
  ))
  
}

get_gene_fams <- function(gene,family_set,inter_family="INTER"){
  gene_pos_in_fams <- sapply(family_set$raw_gene_list,function(x) which(x==gene))
  fams <- names(which(lapply(gene_pos_in_fams,length)>0))
  if(length(fams)==0){
    fams <- inter_family
  }
  return(fams)
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


paint_global_scores <- function(score_runs){
  nscores <- length(score_runs)
  print(names(score_runs))
  par(mfrow=c(nscores,1))
  for(i in 1:nscores){
    paint_global_score(score_runs[[i]],names(score_runs)[i])
  }
}

paint_global_score <- function(score_run,label="Score"){
  
  all_quantiles <- do.call("rbind",lapply(score_run$global_evaluation_list,function(x) x$genes_quantile))
  boxplot(all_quantiles,main=label,cex.axis=0.7,names=F,ylim=c(0,1))
  labels = colnames(all_quantiles)
  text(1:length(labels), par("usr")[3] - 0.10, cex = 0.7,srt = 45, adj = 1, labels, xpd = TRUE)
  lines(c(0,length(labels)+1),c(0.5,0.5),lty=2)
  
  score_summary <- apply(all_quantiles,2,summary)
  print(score_summary)
  return(score_summary)
  
}

paint_score_density <- function(score_run=NULL,target_scores=NULL,random_scores=NULL,inter_scores=NULL,score_name=NULL,label="Score density",norm_densities=F,inter=F){
  
  thecols <- c("grey","cornflowerblue","darkolivegreen4","brown2")
  names(thecols) <- c("all","inter","random","target")
  
  if(is.null(score_name)){
    score_name <- colnames(score_run$multi_score_list[[1]]$global_score_table)[1]
    print(score_name)
  }
  
  if(is.null(score_run)){
    all_scores <- c(target_scores,random_scores,inter_scores)
  } else {
    if(inter==T){
      all_scores <- do.call("rbind",lapply(score_run$multi_score_list,function(x) x$global_score_table_with_inter))
      inter_scores <- do.call("rbind",lapply(score_run$global_evaluation_with_inter_list,function(x) x$inter_global_score_table))
      random_scores <- do.call("rbind",lapply(score_run$global_evaluation_with_inter_list,function(x) x$random_global_score_table))
      target_scores <- do.call("rbind",lapply(score_run$global_evaluation_with_inter_list,function(x) x$target_global_score_table))
    } else {
      all_scores <- do.call("rbind",lapply(score_run$multi_score_list,function(x) x$global_score_table))
      inter_scores <- do.call("rbind",lapply(score_run$global_evaluation_list,function(x) x$inter_global_score_table))
      random_scores <- do.call("rbind",lapply(score_run$global_evaluation_list,function(x) x$random_global_score_table))
      target_scores <- do.call("rbind",lapply(score_run$global_evaluation_list,function(x) x$target_global_score_table))
      print(summary(random_scores))
    }
  }
  
  all_density <- density(all_scores[,score_name])
  if(nrow(inter_scores)>0){
    inter_density <- density(inter_scores[,score_name])
  }
  random_density <- density(random_scores[,score_name])
  target_density <- density(target_scores[,score_name])
  
  if(norm_densities){
    if(nrow(inter_scores)>0){
      inter_ratio <- nrow(inter_scores)/nrow(all_scores)
    }
    random_ratio <- nrow(random_scores)/nrow(all_scores)
    target_ratio <- nrow(target_scores)/nrow(all_scores)    
  } else {
    inter_ratio <- 1
    random_ratio <- 1
    target_ratio <- 1
  }
  
  max_density <- max(c(all_density$y,random_density$y,target_density$y))*1.1
  
  layout(matrix(c(1,1,1,1,2,2,2),ncol=1))
  par(mar=c(2,4,3,3))
  hist(all_scores[,score_name],probability=T,50,xlab="",ylab="density",main=label,ylim=c(0,max_density))
  lines(all_density$x,all_density$y,col=thecols[1],lty=2)
  if(nrow(inter_scores)>0){
    lines(inter_density$x,inter_density$y*inter_ratio,col=thecols[2],lwd=2)
  }
  lines(random_density$x,random_density$y*random_ratio,col=thecols[3],lwd=2)
  lines(target_density$x,target_density$y*target_ratio,col=thecols[4],lwd=2)
  legend("topright",legend=c("all","inter","random","target"),col=thecols,lty=c(2,1,1,1),lwd=2,cex=0.7)
  par(mar=c(5,4,2,3))
  boxplot(at=c(4,3,2,1),xlab=score_name,horizontal=T,list(all_scores[,score_name],inter_scores[,score_name],random_scores[,score_name],target_scores[,score_name]),cex.axis=0.7,names=F,col=thecols)
  labels <- c("INTER","RANDOM","TARGET")
#   text(1:length(labels), par("usr")[3] - 0.05, cex = 0.7,srt = 45, adj = 1, labels, xpd = TRUE)
  
}

get_target_list_distances <- function(target_list,random_list=NULL,distance_matrix){
  
    per_element <- function(i){
      if(is.null(random_list)){
        get_target_distances(target_list[[i]],distance_matrix=distance_matrix)
      } else {
        get_target_distances(target_list[[i]],random_list[[i]],distance_matrix=distance_matrix)
      }
    }    
    per <- lapply(1:length(target_list),per_element)
    result <- list(target_to_target = unlist(sapply(per,function(x)x$target_to_target,simplify=T)))
    if(!is.null(random_list)){
      result[["target_to_random"]] <- unlist(sapply(per,function(x)x$target_to_random,simplify=T))
    }
    result[["target_to_genome"]] <- unlist(sapply(per,function(x)x$target_to_genome,simplify=T))
    
    return(result)
}

get_target_distances <- function(raw_target,random=NULL,distance_matrix){
  
  target <- raw_target[raw_target %in% rownames(distance_matrix)]
  
  target_to_target_distance_matrix <- distance_matrix[target,target]
  target_to_target_distances <- as.integer(target_to_target_distance_matrix)
  target_to_target_distances <- target_to_target_distances[target_to_target_distances>0]
  
  if(!is.null(random)){
    target_to_random_distance_matrix <- distance_matrix[target,random]
    target_to_random_distances <- as.integer(target_to_random_distance_matrix)
    target_to_random_distances <- target_to_random_distances[target_to_random_distances>0]
    genome <- setdiff(colnames(distance_matrix),c(target,random))
  } else {
    genome <- setdiff(colnames(distance_matrix),target) 
  }
  
  target_to_genome_distance_matrix <- distance_matrix[target,genome]
  target_to_genome_distances <- as.integer(target_to_genome_distance_matrix)
  target_to_genome_distances <- target_to_genome_distances[target_to_genome_distances>0]
  
  result <- list(target_to_target = target_to_target_distances)
  if(!is.null(random)){
    result[["target_to_random"]] <- target_to_random_distances
  }
  result[["target_to_genome"]] <- target_to_genome_distances
  
  return(result)
  
}


