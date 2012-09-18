

create_score_data <- function(sre_file,distance,output=F,out.file="score_data.txt"){
  
  # get data
  load(sre_file)
  sre_name <- substr(sre_file,regexpr("sre",sre_file),regexpr(".rdata",sre_file)-1)
  sre <- get(sre_name)
  rm(list=sre_name)
  gc()
  
  # create data frame
  score_data <- c()
  for(i in 1:sre$params$nsims){
      
    gst <- sre$score_runs$final_score$multi_score_list[[i]]$global_score_table
    target <- rep(0,nrow(gst))
    target[which(rownames(gst) %in% sre$family_set_series[[i]]$fam_disease_genes)] <- 1
    target <- as.factor(target)
    levels(target) <- c("0"="RANDOM","1"="TARGET")
  
    this_data <- cbind(
      "gene"=rownames(gst),
      gst,
      "type"=target,
      "repetition"=i,
      "distance"=distance,
      "nfams"=sre$family_set_series[[i]]$nfams,
      "ndisease_genes"=sre$family_set_series[[i]]$ndisease_genes,
      "ngenes"=sre$family_set_series[[i]]$ngenes,
      "sre_file"=sre_file
    )
    
    score_data <- rbind(score_data,this_data)
    
  }
  if(output){
    write.table(score_data,sep="\t",quote=F,row.names=F,col.names=T,file=paste(dirname(sre_file),"/",basename(dirname(sre_file)),out.file,sep=""))
  }
  
  return(score_data)
  
}


filter_data <- function(data=NULL,nfams=NULL,ngenes=NULL,ndisease_genes=NULL,distance=NULL,sim_set=NULL,type="RANDOM"){
  
  # number of families
  if(is.null(nfams)){
    nfams_filter <- rep(TRUE,nrow(data))
  } else {
    nfams_filter <- data$nfams==nfams
  }
  
  # number of genes per family
  if(is.null(ngenes)){
    ngenes_filter <- rep(TRUE,nrow(data))
  } else {
    ngenes_filter <- data$ngenes==ngenes
  }
  
  # number of disease genes
  if(is.null(ndisease_genes)){
    ndisease_genes_filter <- rep(TRUE,nrow(data))
  } else {
    ndisease_genes_filter <- data$ndisease_genes==ndisease_genes
  }
  
  # distance
  if(is.null(distance)){
    distance_filter <- rep(TRUE,nrow(data))
  } else {
    distance_filter <- data$distance==distance
  }
  
  # sim_set
  if(is.null(sim_set)){
    sim_set_filter <- rep(TRUE,nrow(data))
  } else {
    sim_set_filter <- rep(FALSE,nrow(data))
    sim_set_filter[grep(sim_set,data$sre_file)] <- TRUE
  }
  
  # gene type
  if(is.null(type)){
    type_filter <- rep(TRUE,nrow(data))
  } else {
    type_filter <- data$type==type
  }
  
#   print(summary(nfams_filter))
#   print(summary(ngenes_filter))
#   print(summary(ndisease_genes_filter))
#   print(summary(distance_filter))
#   print(summary(sim_set_filter))
#   print(summary(type_filter))
  
  data[nfams_filter & ngenes_filter & ndisease_genes_filter & distance_filter & sim_set_filter & type_filter,]
  
  
}


compare_distributions <-function(data,nfams_vector,ngenes_vector,ndisease_genes_vector){
  
  
  
  
  
  
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, by.row=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





