
source("lib/empirica_utils.r")
require(ggplot2)

all_experiments <- read.table("charts_data/all_experiments.txt",sep="\t",header=T,stringsAsFactors=F)


# global by distance
png("charts_data/emp_tabla_global.png",width=700,heigh=400)
ggplot(all_experiments, aes(no_zero_mean, colour = factor(distance))) + geom_density() + geom_histogram() + opts(title="Distribución global")
dev.off()

# tabla distancia
tabla_distancia <- filter_data(all_experiments,sim_set="tabla_distancia")
png("charts_data/emp_tabla_distancia.png",width=700,heigh=400)
ggplot(tabla_distancia, aes(no_zero_mean, colour = factor(distance))) +  geom_density() + geom_histogram() + opts(title="Tabla distancias")
dev.off()

# tabla familias
tabla_familias_sp <- filter_data(all_experiments,distance="SP",sim_set="tabla_familias")
tabla_familias_inter <- filter_data(all_experiments,distance="IN",sim_set="tabla_familias")
tabla_familias_rw <- filter_data(all_experiments,distance="RW",sim_set="tabla_familias")
tabla_familias_rw_all <- filter_data(all_experiments,distance="RW",sim_set="tabla_familias",type=NULL)
  # sp
png("charts_data/emp_tabla_familias_sp.png",width=700,heigh=400)
ggplot(tabla_familias_sp, aes(no_zero_mean, colour = factor(nfams))) +  geom_density() + opts(title="Shortest Path")
dev.off()
  # inter
png("charts_data/emp_tabla_familias_inter.png",width=700,heigh=800)
p1 <- ggplot(tabla_familias_inter, aes(no_zero_mean, colour = factor(nfams))) +  geom_density() + opts(title="Intermediación")
p2 <- ggplot(tabla_familias_inter[which(tabla_familias_inter$no_zero_mean<0.1),], aes(no_zero_mean, colour = factor(nfams))) +  geom_density()
multiplot(p1,p2)
dev.off()
  # rw
png("charts_data/emp_tabla_familias_rw.png",width=700,heigh=800)
p1 <- ggplot(tabla_familias_rw, aes(no_zero_mean, colour = factor(nfams))) +  geom_density() + opts(title="Random walk")
p2 <- ggplot(tabla_familias_rw[which(tabla_familias_rw$no_zero_mean<0.1),], aes(no_zero_mean, colour = factor(nfams))) +  geom_density()
multiplot(p1,p2)
dev.off()
p1 <- ggplot(tabla_familias_rw, aes(no_zero_mean, colour = factor(nfams))) +  geom_density() + opts(title="Random walk")
p1 + scale_x_log10()

# tabla genes
tabla_genes_sp <- filter_data(all_experiments,distance="SP",sim_set="tabla_genes")
tabla_genes_inter <- filter_data(all_experiments,distance="IN",sim_set="tabla_genes")
tabla_genes_rw <- filter_data(all_experiments,distance="RW",sim_set="tabla_genes")
  # sp
ggplot(tabla_genes_sp, aes(no_zero_mean, colour = factor(ngenes))) +  geom_density() + opts(title="Shortest Path")
  # inter
p1 <- ggplot(tabla_genes_inter, aes(no_zero_mean, colour = factor(ngenes))) +  geom_density() + opts(title="Intermediación")
p2 <- ggplot(tabla_genes_inter[which(tabla_genes_inter$no_zero_mean<0.1),], aes(no_zero_mean, colour = factor(ngenes))) +  geom_density()
multiplot(p1,p2)
  # rw
p1 <- ggplot(tabla_genes_rw, aes(no_zero_mean, colour = factor(ngenes))) +  geom_density() + opts(title="Intermediación")
p2 <- ggplot(tabla_genes_rw[which(tabla_genes_rw$no_zero_mean<0.02),], aes(no_zero_mean, colour = factor(ngenes))) +  geom_density()
multiplot(p1,p2)


ggplot(filter_data(all_experiments,distance="SP",sim_set="tabla_distancia",type=NULL),aes(no_zero_mean, colour = factor(type))) +  geom_density() + opts(title="Shortest path")
ggplot(filter_data(all_experiments,distance="IN",sim_set="tabla_distancia",type=NULL),aes(no_zero_mean, colour = factor(type))) +  geom_density() + opts(title="Intermediación")
ggplot(filter_data(all_experiments,distance="RW",sim_set="tabla_distancia",type=NULL),aes(no_zero_mean, colour = factor(type))) +  geom_density() + opts(title="Random walk")

ggplot(filter_data(all_experiments,distance="SP",sim_set="tabla_distancia",type=NULL),aes(no_zero_max, colour = factor(type))) +  geom_density() + opts(title="Shortest path")
ggplot(filter_data(all_experiments,distance="IN",sim_set="tabla_distancia",type=NULL),aes(no_zero_max, colour = factor(type))) +  geom_density() + opts(title="Intermediación")
ggplot(filter_data(all_experiments,distance="RW",sim_set="tabla_distancia",type=NULL),aes(no_zero_max, colour = factor(type))) +  geom_density() + opts(title="Random walk")


sp_den <- density(all_experiments[which(all_experiments$distance=="SP"),]$no_zero_mean,na.rm=T)
inter_den <- density(all_experiments[which(all_experiments$distance=="IN"),]$no_zero_mean,na.rm=T)
rw_den <- density(all_experiments[which(all_experiments$distance=="RW"),]$no_zero_mean,na.rm=T)
hist(all_experiments$no_zero_mean,100,probability=T)
lines(sp_den,col="red")
lines(inter_den,col="green")
lines(rw_den,col="blue")


require(cmaes)

paint_optimization <- function(scale,shape){
  
  model_den_vector <- dgamma(x,scale=scale,shape=shape)   
  data_den_vector <- data_den$y  
  max_den <- max(data_den_vector,model_den_vector,na.rm=T)
  data_den_norm <- data_den_vector/max_den
  model_den_norm <- model_den_vector/max_den
  
  layout(matrix(c(1,1,1,1,2,3),ncol=2,byrow=3))
  plot(x,data_den_vector,col="black",ylim=c(0,max_den),type="l")
  lines(x,model_den_vector,col="red",type="l")
  plot(x,data_den_vector,col="black",type="l")
  plot(x,model_den_vector,col="red",type="l")
}

get_diff <- function(params,paint=T) {
 
  iter <<- iter + 1
  print(iter)
  
  scale <- params["scale"]
  shape <- params["shape"]
  
  x <- data_den$x  
  
  model_den_vector <- dgamma(x,scale=scale,shape=shape)
  data_den_vector <- data_den$y

  max_den <- max(data_den_vector,model_den_vector,na.rm=T)
  data_den_norm <- data_den_vector/max_den
  model_den_norm <- model_den_vector/max_den
  
  if(paint){
    plot(x,data_den_vector,col="black",ylim=c(0,max_den),type="l")
    lines(x,model_den_vector,col="red")
    Sys.sleep(0.15)
  }
  
#   print(max_den)
#   print(data_den_vector)
#   print(data_den_norm)
  
  error <- 1- cor(data_den_norm,model_den_norm)^2 #sum(abs(data_den_norm-model_den_norm))
  infos[[iter]] <<- c(scale,shape,error=error)
  
  return(error)
}

data <- filter_data(tabla_familias_rw,nfams=5)$no_zero_mean
data_den <- density(data,na.rm=T)

dev.off()
iter <- 0
infos <- list()
system.time(opti <- cma_es(c(scale=mean(data,na.rm=T),shape=1),get_diff,lower=c(scale=0,shape=0),upper=c(scale=2,shape=2),control=list(maxit=20,sigma=0.001)))
opti$par
sorted_infos <- infos[order(sapply(infos,function(x) x["error"]))]

paint_optimization(scale=opti$par[1],shape=opti$par[2])






get_diff_beta <- function(params,paint=T,debug=F) {
  
  iter <<- iter + 1
  print(iter)
  
  shape1 <- params[1]
  shape2 <- params[2]
  
  x <- data_den$x  
  
  print("preparado")
  print(x)
  print(params)
  print(shape1)
  print(shape2)
  
  model_den_vector <- dbeta(x,shape1=shape1,shape2=shape2)   
  data_den_vector <- data_den$y
  
  max_den <- max(data_den_vector,model_den_vector,na.rm=T)
  data_den_norm <- data_den_vector/max_den
  model_den_norm <- model_den_vector/max_den
  
  if(paint){
    plot(x,data_den_vector,col="black",ylim=c(0,max_den),type="l")
    lines(x,model_den_vector,col="red")
    Sys.sleep(0.15)
  }
    
  error <- 1- cor(data_den_norm,model_den_norm)^2 #sum(abs(data_den_norm-model_den_norm))
  infos[[iter]] <<- c(shape1,shape2,error=error)
  
  if(debug){
    return(list(error=error,x=x))
  } else {
    return(error)
  }
  
}
get_diff_beta(sorted_infos[[2]])

dev.off()
iter <- 0
infos <- list()
system.time(opti <- cma_es(c(shape1=mean(data,na.rm=T),shape2=1),get_diff_beta,lower=c(shape1=0,shape2=0),upper=c(shape1=5,shape2=5),control=list(maxit=20)))
opti$par
sorted_infos <- infos[order(sapply(infos,function(x) x["error"]))]
get_diff_beta(sorted_infos[[1]])






ggplot(tabla_familias_rw_all, aes(no_zero_max, colour = factor(nfams), fill=factor(nfams))) +  geom_density(alpha=.3) + opts(title="Random walk") + facet_grid(type~.) + scale_x_log10() + theme_bw()



