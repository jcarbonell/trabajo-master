
nh_simulate_genes <- function(families, overlapping, sif, inter_label="inter",family_proportion=0.10,byseed=T,seed_radio=5) {
    
  # init
  nfams <- length(families)  
  total_fam_genes <- sum(diag(overlapping))    
  initial_size <- round(total_fam_genes/family_proportion)
    
  # init reference subnet
  all_interactors <- unique(rbind(sif$V1,sif$V3))
  
  # select connected genes
  if(byseed){
    seed <- sample(all_interactors)[1]    
    maxiter <- 100
    iter <- 0
    radio <- seed_radio
    ref_subnet_genes <- c()
    while(length(ref_subnet_genes)<initial_size | iter> maxiter){
      nhood <- nh_get_neiborhood(seed,sif,radio=radio)
      ref_subnet_genes <- nh_get_nodes(nhood)      
      iter <- iter +1
      radio = radio +1
    }
  } else {
    interactome_subset <- sample(all_interactors)[1:initial_size]
    ref_subnet <- get.all.shortest.paths.Luz(sif, interactome_subset, 0)
    ref_subnet_genes <- unique(c(ref_subnet$V1,ref_subnet$V3))   
  }   
  ng <- length(ref_subnet_genes)

  # label connected genes
  ref_subnet_genes_sample <- sample(ref_subnet_genes)    
  gene_list <- list()
  cont <- 0
  for(i in 1:nfams){
    ngenes <- overlapping[i,i]*2 - sum(overlapping[,i])
    if(ngenes>0){
      gene_list[[i]] <- ref_subnet_genes_sample[(cont+1):(cont+ngenes)]
    } else {
      gene_list[[i]] <- list()
    }
    cont <- cont + ngenes
  }
  names(gene_list) <- families

  # overlapping
  for(i in 1:nfams){
    for(j in i:nfams){
      if(i!=j){
        if(overlapping[i,j]>0){
          ngenes <- overlapping[i,j]
          overlapped_genes <- ref_subnet_genes_sample[(cont+1):(cont+ngenes)]
          gene_list[[i]] <- c(gene_list[[i]],overlapped_genes)
          gene_list[[j]] <- c(gene_list[[j]],overlapped_genes)
          cont <- cont + ngenes
        }
      }      
    }
  }

  return(gene_list)

}

paint_family_graph <- function(sif,gene_scores,gene_fams,highlighteds=NULL,min_score=NULL,max_score=NULL,paint_labels=T,inter_color="#AAAAAA",inter_label="inter",title="",layout=NULL){
  
  # prepare border colors
  all_sif_genes <- unique(c(sif$V1,sif$V3))
  
  # recover and label inter
  
  sif_scores <- gene_scores[all_sif_genes]  
  names(sif_scores) <- all_sif_genes
  sif_scores[is.na(sif_scores)] <- 0
  
  all_gene_fams <- gene_fams[all_sif_genes]
  names(all_gene_fams) <- all_sif_genes
  all_gene_fams[is.na(all_gene_fams)] <- inter_label
  
  # normalize score
  if(is.null(min_score)){
    themin <- min(sif_scores)
  } else {
    themin <- min_score
  }
  if(is.null(max_score)){
    themax <- max(sif_scores)
  } else{
    themax <- max_score
  }
  if(themin==themax & is.null(min_score) & is.null(max_score)){
    nscores <- rep(0.5,length(sif_scores))
  } else {
    centered_nscores <- sif_scores-themin
    centered_nscores[which(centered_nscores<0)] <- 0  
    nscores <- centered_nscores/(themax-themin)
  }
   
  
  # colors
  gene_score_colors <- rgb(nscores,nscores,nscores)
  names(gene_score_colors) <- names(sif_scores)
    
  fams <- unique(gene_fams)
  fam_colors <- rainbow(length(fams))
  names(fam_colors) <- fams  
  fam_colors[inter_label] <- inter_color
  gene_fam_colors <- fam_colors[all_gene_fams]
  names(gene_fam_colors) <- names(all_gene_fams)
      
  # paint graph  
  paint_igraph(sif,gene_score_colors,gene_fam_colors,fam_colors,highlighteds=highlighteds,paint_labels=paint_labels,title=title,layout=layout)
  
}

nh_simulate <- function(family_names,overlapping_matrix,radio,randoms,interactome,ninterm=2,confidence=c(0.75,0.95,0.99)) {
    
  # INIT ###############################################
  
  # working variables
  nfams <- length(family_names)  
  rownames(overlapping_matrix) <- family_names
  colnames(overlapping_matrix) <- family_names
    
  # SIMULATION #########################################
  
  # simulate genes
  cat("   Simulating genes\n")
  initial_family_genes <- nh_simulate_genes(family_names,overlapping_matrix,interactome,family_proportion=0.4,seed_radio=radio)
  initial_gene_list <- unique(unlist(initial_family_genes))
    
  # add random genes
  family_genes <- list()
  for(i in 1:nfams){
    family_genes[[family_names[i]]] <- unique(c(initial_family_genes[[i]],sample(interactome_genes,randoms[i])))
  }
  gene_list <- unique(unlist(family_genes))
  
  
  # COMPUTE SCORE #######################################
  
  cat("   Testing simulated genes\n")
  
  # recompute interactome
  cat("       Scanning simulated gene interactome\n")  
  subnet <- get.all.shortest.paths.Luz(interactome, gene_list, ninterm)
  
  # compute score
  cat("       Computing scores\n")
  scores <- nh_score(family_genes,subnet)
  
  # paint network
#   layout <- nh_get_standard_layout(subnet)
#   paint_family_graph(subnet,scores$pnh[gene_list],scores$hits,title="Simulated genes",paint_labels=T,layout=layout)
#   paint_family_graph(subnet,scores$nfc[gene_list],scores$hits,title="Simulated genes",paint_labels=T,layout=layout)
  
  # paint without inter-inter edges
#   clean_subnet <- nh_remove_inter_inter_interactions(subnet,scores$hits)
#   paint_family_graph(clean_subnet,scores$pnh,scores$hits,title="Simulated genes",paint_labels=T,layout=layout)
  
  
  # THRESHOLD SCORE ######################################
  
  cat("       Selecting best scored genes\n") 
  
  p <- length(initial_gene_list)
  n <- length(gene_list) - p
 
  nfc_validation <- nh_test_score(scores$nfc,confidence,initial_family_genes,gene_list,p,n)
  phn_validation <- nh_test_score(scores$phn,confidence,initial_family_genes,gene_list,p,n)
   
  # RESULTS ###############################################  
  gene_frame <- list(
    initial_family_genes=initial_family_genes,
    initial_gene_list=initial_gene_list,
    family_genes=family_genes,
    gene_list    
  )
  
  nconfidence <- length(confidence)
    
  nfc_validation_frame <- unlist(nfc_validation)
  names(nfc_validation_frame) <- paste(rep(names(nfc_validation),each=nconfidence),confidence,sep="_")
    
  phn_validation_frame <- unlist(phn_validation)
  names(phn_validation_frame) <- paste(rep(names(phn_validation),each=nconfidence),confidence,sep="_")
  
  validation_frame <- list(
    confidence=confidence,
    nfc=nfc_validation_frame,
    phn=phn_validation_frame
  )
    
  score_frame <- scores
  
  return(list(
    gene_frame=gene_frame,
    validation_frame=validation_frame,
    score_frame=score_frame
    )
  )
  
}


nh_test_score <- function(score,confidence,initial_family_genes,gene_list,p,n){
  
  initial_gene_list <- unique(unlist(initial_family_genes))
  nconf <- length(confidence)
    
  tp <- numeric(nconf)
  fn <- numeric(nconf)
  fp <- numeric(nconf)
  tn <- numeric(nconf)
  tpr <- numeric(nconf)
  fnr <- numeric(nconf)
  fpr <- numeric(nconf)
  tnr <- numeric(nconf)
  
  for(i in 1:nconf){
  
    nh_classif <- nh_quantile_classify(score,confidence[i])
    # nh_paint_classification(nh_classif)
    
    selected_family_genes <- gene_list[which(is.element(gene_list,nh_classif$selected_genes))]
        
    correct_genes <- selected_family_genes[which(is.element(selected_family_genes,initial_family_genes))]
    corrects <- length(correct_genes)
    correct_ratio <- round(digits=3,(corrects/length(initial_gene_list))*100)
         
    tp[i] <- corrects
    fn[i] <- p - tp[i]
    fp[i] <- length(selected_family_genes) - tp[i]
    tn[i] <- n - fp[i]    
    tpr[i] <- ( tp[i] / p ) * 100
    fnr[i] <- ( fn[i] / n ) * 100
    fpr[i] <- ( fp[i] / p ) * 100
    tnr[i] <- ( tn[i] / n ) * 100
    
  }
 
  names(tp) <- confidence
  names(fn) <- confidence
  names(fp) <- confidence
  names(tn) <- confidence
  names(tpr) <- confidence
  names(fnr) <- confidence
  names(fpr) <- confidence
  names(tnr) <- confidence
    
  list(tp=tp,fn=fn,fp=fp,tn=tn,tpr=tpr,fnr=fnr,fpr=fpr,tnr=tnr)
    
}


    