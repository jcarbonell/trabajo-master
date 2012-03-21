
source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("simulation.r")

home_path <- getwd()

# input params

nfams <- 3
ngenes <- 20
nsims <- 20

# load interactomes

interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactome <- interactomes[["binding"]]
distance_matrix <- compute_distance_matrix(interactome)

# simulation

fam_simulation_list <- list()
scores_frame_list <- list()
evaluation_list <- list()

for(k in 1:nsims){

  cat("\n>>simulation step",k,"...\n")
  
  # simulate family genes
  cat("  simulating family genes...\n")
  fam_simulation_list[[k]] <- simulate_family_genes(nfams,ngenes,distance_matrix=distance_matrix)
  
  # score computation
  cat("  computing score...\n")
  scores_frame_list[[k]] <- compute_score(fam_simulation_list[[k]]$raw_gene_list,interactome,default_score2)
  
  # evaluation
  cat("  evaluating prioritization...\n")
  evaluation_list[[k]] <- evaluate_prioritization(scores_frame_list[[k]],fam_simulation_list[[k]],paint=F)
  
  cat("OK\n")
  
}

# relative position plot
mean_relative_position <- sapply(evaluation_list,function(x) x$mean_relative_position)
mean_relative_position_with_inter <- sapply(evaluation_list,function(x) x$mean_relative_position_with_inter)
plot(mean_relative_position,col="red",ylim=c(0,1),pch=20,xlab="step",ylab="relative position")
points(mean_relative_position_with_inter,col="blue",ylim=c(0,1),pch=20)
legend(nsims*0.05,0.2,legend=c("input genes","input genes with intermediates"),col=c("red","blue"),cex=0.6,lwd=1)
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

# score distribution
fam_genes_score <- as.vector(sapply(evaluation_list,function(x) x$fam_genes_score))
all_scores <- unlist(sapply(scores_frame_list, function(x) x$score_table$score))
paint_evaluation(all_scores,fam_genes_score,"All scores")

h <- hist(all_scores,35,)
hfg <- hist(fam_genes_score,35,add=T,col="red")
d <- density(all_scores)
lines(d,col="darkgray")
hfg <- hist(fam_genes_score,plot=F)
dfg <- density(fam_genes_score)
lines(dfg$x,dfg$y*max(hfg$counts),col="red")



