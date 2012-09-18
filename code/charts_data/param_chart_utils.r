
get_quantiles <- function(files){
  
  n <- length(files)
  
  sre_quantiles <- list()
  for(i in 1:n){
    
    cat("Loading ",files[i],"...\n")
    load(files[i])
    sre_name <- substr(files[i],regexpr("sre",files[i]),regexpr(".rdata",files[i])-1)
    
    sre_quantiles[[i]] <- unlist(lapply(get(sre_name)$score_runs$final_score$global_evaluation_list,function(x) x$genes_quantile[,"no_zero_mean"]))
    
    rm(list=sre_name)
    gc()
    
  }  
  
  sre_quantiles_mean <- sapply(sre_quantiles,mean)
  sre_quantiles_sd <- sapply(sre_quantiles,sd)
  
  return(list(mean=sre_quantiles_mean,sd=sre_quantiles_sd,quantiles=sre_quantiles))
         
}


param_chart <- function(quantiles_list,values,labels,xlabel,thetitle){
  
  n <- length(quantiles_list)
  
  colors <- c("red","blue","pink")
  
  paint_area <- function(quantiles,color){
    col_values <- col2rgb(color)/255
    fill_color <- rgb(col_values[1],col_values[2],col_values[3],0.2)
    polygon(border=0,col=fill_color, y=c(quantiles$mean+quantiles$sd,rev(quantiles$mean-quantiles$sd)),x=c(values,rev(values)))        
  }
  
  plot(values,quantiles_list[[1]]$mean,col=colors[1],type="l",ylim=c(0,1),xlab=xlabel,ylab="posición relativa media",main=thetitle)
  points(values,quantiles_list[[1]]$mean,col=colors[1])  
  paint_area(quantiles_list[[1]],colors[1])
      
  if(n>1){
    for(i in 2:n){
      lines(values,quantiles_list[[i]]$mean,col=colors[i])      
      points(values,quantiles_list[[i]]$mean,col=colors[i])
  
      paint_area(quantiles_list[[i]],colors[i])
      
    }    
  }
  
  legend("bottomright",legend=labels,col=colors,lwd=1)
}

param_distribution <- function(all_data,xlabel,thetitle){
  
  p <- ggplot(all_data, aes(x=factor(values), y = rank, fill = distancia)) + geom_boxplot(outlier.size=0.1)
  p + opts(title=thetitle) + 
    xlab(xlabel) + 
    ylab("posición relativa en el ranking") + 
    opts(plot.margin = unit(1*c(3, 3, 3, 3), "lines")) +
    opts(axis.title.x=theme_text(vjust=-2,size=13)) + 
    opts(axis.title.y=theme_text(angle=90, vjust=-0.5,size=13)) + 
    opts(plot.title=theme_text(size=16, vjust=3))
  
  
    
}

