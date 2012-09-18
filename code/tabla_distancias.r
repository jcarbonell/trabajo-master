
source("lib/shortestpath-mariangela.R")
source("lib/score.r")
source("lib/score_functions.r")
source("lib/utils.r")
source("lib/simulation.r")
source("lib/simulation_run.r")
source("lib/evaluation.r")
source("lib/postprocess.r")
source("lib/random_walk.r")

# previous to executions
home_path <- getwd()
work_path <- paste(home_path,"/work/tabla_distancia",sep="")
dir.create(work_path)

# family set series
params <- list(nsims=100,nfams=3,ndisease_genes=3,ngenes=150)
fss <- simulate_known_family_set_series(nsims=params$nsims,nfams=params$nfams,ndisease_genes=params$ndisease_genes,ngenes=params$ngenes)

# shortest_path
sre_sp <- do_shortest_path_simulation(fss,paste(work_path,"/",get_simulate_id(params,prefix="SP"),sep=""),params,get_back=T)
#sre_inter <- do_intermediation_simulation(fss,paste(work_path,"/",get_simulate_id(params,prefix="INTER"),sep=""),params,get_back=T)
sre_rw <- do_random_walk_simulation(fss,paste(work_path,"/",get_simulate_id(params,prefix="RW"),sep=""),params,get_back=T)

