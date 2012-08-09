
source("shortestpath-mariangela.R")
source("score.r")
source("utils.r")
source("simulation.r")
source("evaluation.r")
source("postprocess.r")

home_path <- getwd()

# input params
nfams <- 3
ngenes <- 50
nsims <- 10

# load interactomes
raw_interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
interactomes <- list(binding=raw_interactomes$binding)
meta_interactome  <- get_meta_interactome(interactomes)

# known genes simulation
# all_genes <- unique(unique(meta_interactome$V1),unique(meta_interactome$V3))
# family_set_series <- simulate_known_family_set_series(paste(home_path,"/../misc/omim_clean_annot.txt",sep=""),all_genes,nfams=nfams,ngenes=ngenes,nsims=nsims)

# dd <- get_target_distances(target,random,distance_matrix=meta_interactome_matrix)

# pure simulated genes
meta_interactome_matrix <- compute_distance_matrix(meta_interactome[sample(1:nrow(meta_interactome),50000),])
family_set_series <- simulate_family_set_series(nfams=nfams,ngenes=ngenes,nsims=nsims,distance_matrix=meta_interactome_matrix)

sim_dis <- get_target_distances(target,random,distance_matrix=meta_interactome_matrix)


dannot <- read.delim("../misc/omim_clean_annot.txt",sep="\t",header=T,stringsAsFactors=F)
disease_genes <- lapply(by(dannot,dannot$disease,function(x) x$gene,simplify=T),function(x)x)
genes_per_disease <- sapply(disease_genes,length)
sel_disease_genes <- disease_genes[genes_per_disease>3]

tld <- get_target_list_distances(target_list=sel_disease_genes,distance_matrix=meta_interactome_matrix)
target_to_target_freqs <- table(tld$target_to_target)
target_to_genome_freqs <- table(tld$target_to_genome)
max_dis <- max(c(as.integer(names(target_to_target_freqs)),as.integer(names(target_to_genome_freqs))))
par(mfrow=c(2,1))
barplot(target_to_target_freqs[1:max_dis],names=1:max_dis,main="Target to target distance")
barplot(target_to_genome_freqs[1:max_dis],names=1:max_dis,main="Target to genome distance")



