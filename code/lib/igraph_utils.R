
paint_igraph <- function(sif,node_color_ref=NULL,border_color_ref=NULL,label_colors=NULL,highlighteds=NULL,paint_labels=T,title="",layout=NULL){
  
  raw_nodes <- nh_get_nodes(sif)
  
  if(is.null(node_color_ref)){
    node_color_ref <- rep("blue",length(raw_nodes))
    names(node_color_ref) <- raw_nodes
  }
  if(is.null(border_color_ref)){
    border_color_ref <- rep("black",length(raw_nodes))
    names(border_color_ref) <- raw_nodes    
  }

  # init graph
  network <- graph.data.frame(sif[,c(1,3)])
    
  # define graph attributes
  V(network)$label <- V(network)$name
  V(network)$label.font <- 1
  V(network)$label.cex <- 0.5
  V(network)$label.color <- rgb(0,0,0,0.5)
  E(network)$arrow.mode <- ">"
  E(network)$width <- 1
  
  node_list <- V(network)$name
  node_colors <- node_color_ref[node_list]  
  border_colors <- border_color_ref[node_list]

  # layout  
  if(is.null(layout)){
    mylayout <- layout.graphopt(network,niterInteger=1000)
  } else {
    mylayout <- layout[node_list,]
  }
  
  # other params
  if(paint_labels==F){
    V(network)$label <- NA
  } 
    
  # plot graph
  plot(network,layout=mylayout, vertex.color=border_colors,vertex.label.dist=0.4,vertex.size=5,edge.arrow.size=0.05)
  plot(network,layout=mylayout, edge.color=NA,vertex.color=node_colors,vertex.label.dist=0.4,vertex.size=3,edge.arrow.size=0.05,add=T)
  
  if(!is.null(label_colors)){
    legend(par()$xaxp[1],par()$yaxp[2],legend=names(label_colors),col=label_colors,lwd=1,cex=0.55)
  }
  
  title(title)      
  
  if(!is.null(highlighteds)){
    print(highlighteds)
    nmylayout <- layout.norm(mylayout, -1, 1, -1, 1)
    rownames(nmylayout) <- V(network)$name
    radio <- 0.05
    for(i in 1:length(highlighteds)){
      coords <- nmylayout[highlighteds[i],]
      draw_circle(coords[1],coords[2],radio=radio)      
      draw_circle(coords[1],coords[2],radio=radio+0.01)      
#       draw_circle(coords[1],coords[2]-radio*1.25,radio=radio/2,col="white",fill=T)
#       text(coords[1],coords[2]-radio,label=i,col="red",cex=0.75)      
    }
  }
    
}
  
draw_circle <- function(x,y,radio=1,steps=30,col="red",fill=F){
  
  t <- seq(0,2*pi,length=steps) 
  
  coords <- t(rbind(x+sin(t)*radio,y+cos(t)*radio))
  
  if(fill){
    polygon(coords,col=col,border=F)
  } else {
    lines(coords,col=col)
  }
  
}

  