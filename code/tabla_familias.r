
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
work_path <- paste(home_path,"/work/tabla_familias",sep="")
dir.create(work_path)

# family set series
params_3f <- list(nsims=100,nfams=3,ndisease_genes=3,ngenes=100)
params_4f <- list(nsims=100,nfams=4,ndisease_genes=4,ngenes=100)
params_5f <- list(nsims=100,nfams=5,ndisease_genes=5,ngenes=100)
params_7f <- list(nsims=100,nfams=7,ndisease_genes=7,ngenes=100)
params_10f <- list(nsims=100,nfams=10,ndisease_genes=10,ngenes=100)
fss_3f <- simulate_known_family_set_series(params=params_3f)
fss_4f <- simulate_known_family_set_series(params=params_4f)
fss_5f <- simulate_known_family_set_series(params=params_5f)
fss_7f <- simulate_known_family_set_series(params=params_7f)
fss_10 <- simulate_known_family_set_series(params=params_10f)

# shortest_path
do_shortest_path_simulation(fss_3f,paste(work_path,"/",get_simulate_id(params_3f,prefix="SP"),sep=""),params_3f)
do_shortest_path_simulation(fss_4f,paste(work_path,"/",get_simulate_id(params_4f,prefix="SP"),sep=""),params_4f)
do_shortest_path_simulation(fss_5f,paste(work_path,"/",get_simulate_id(params_5f,prefix="SP"),sep=""),params_5f)
do_shortest_path_simulation(fss_7f,paste(work_path,"/",get_simulate_id(params_7f,prefix="SP"),sep=""),params_7f)
do_shortest_path_simulation(fss_10,paste(work_path,"/",get_simulate_id(params_10f,prefix="SP"),sep=""),params_10f)

# # intermediation
# do_intermediation_simulation(fss_3f,paste(work_path,"/",get_simulate_id(params_3f,prefix="INTER"),sep=""),params_3f)
# do_intermediation_simulation(fss_4f,paste(work_path,"/",get_simulate_id(params_4f,prefix="INTER"),sep=""),params_4f)
# do_intermediation_simulation(fss_5f,paste(work_path,"/",get_simulate_id(params_5f,prefix="INTER"),sep=""),params_5f)
# do_intermediation_simulation(fss_7f,paste(work_path,"/",get_simulate_id(params_7f,prefix="INTER"),sep=""),params_7f)
# do_intermediation_simulation(fss_10,paste(work_path,"/",get_simulate_id(params_10f,prefix="INTER"),sep=""),params_10f)

# random walk
do_random_walk_simulation(fss_3f,paste(work_path,"/",get_simulate_id(params_3f,prefix="RW"),sep=""),params_3f)
do_random_walk_simulation(fss_4f,paste(work_path,"/",get_simulate_id(params_4f,prefix="RW"),sep=""),params_4f)
do_random_walk_simulation(fss_5f,paste(work_path,"/",get_simulate_id(params_5f,prefix="RW"),sep=""),params_5f)
do_random_walk_simulation(fss_7f,paste(work_path,"/",get_simulate_id(params_7f,prefix="RW"),sep=""),params_7f)
do_random_walk_simulation(fss_10,paste(work_path,"/",get_simulate_id(params_10f,prefix="RW"),sep=""),params_10f)
