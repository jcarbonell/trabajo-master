
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
work_path <- paste(home_path,"/work/tabla_disease",sep="")
dir.create(work_path)

# family set series
params_5d <- list(nsims=100,nfams=5,ndisease_genes=5,ngenes=100)
params_4d <- list(nsims=100,nfams=5,ndisease_genes=4,ngenes=100)
params_3d <- list(nsims=100,nfams=5,ndisease_genes=3,ngenes=100)
params_2d <- list(nsims=100,nfams=5,ndisease_genes=2,ngenes=100)
params_1d <- list(nsims=100,nfams=5,ndisease_genes=1,ngenes=100)
fss_5d <- simulate_known_family_set_series(params=params_5d)
fss_4d <- simulate_known_family_set_series(params=params_4d)
fss_3d <- simulate_known_family_set_series(params=params_3d)
fss_2d <- simulate_known_family_set_series(params=params_2d)
fss_1d <- simulate_known_family_set_series(params=params_1d)

# shortest_path
do_shortest_path_simulation(fss_5d,paste(work_path,"/",get_simulate_id(params_5d,prefix="SP"),sep=""),params_5d)
do_shortest_path_simulation(fss_4d,paste(work_path,"/",get_simulate_id(params_4d,prefix="SP"),sep=""),params_4d)
do_shortest_path_simulation(fss_3d,paste(work_path,"/",get_simulate_id(params_3d,prefix="SP"),sep=""),params_3d)
do_shortest_path_simulation(fss_2d,paste(work_path,"/",get_simulate_id(params_2d,prefix="SP"),sep=""),params_2d)
do_shortest_path_simulation(fss_1d,paste(work_path,"/",get_simulate_id(params_1d,prefix="SP"),sep=""),params_1d)

# # intermediation
# do_intermediation_simulation(fss_5d,paste(work_path,"/",get_simulate_id(params_5d,prefix="INTER"),sep=""),params_5d)
# do_intermediation_simulation(fss_4d,paste(work_path,"/",get_simulate_id(params_4d,prefix="INTER"),sep=""),params_4d)
# do_intermediation_simulation(fss_3d,paste(work_path,"/",get_simulate_id(params_3d,prefix="INTER"),sep=""),params_3d)
# do_intermediation_simulation(fss_2d,paste(work_path,"/",get_simulate_id(params_2d,prefix="INTER"),sep=""),params_2d)
# do_intermediation_simulation(fss_1d,paste(work_path,"/",get_simulate_id(params_1d,prefix="INTER"),sep=""),params_1d)

# random walk
do_random_walk_simulation(fss_5d,paste(work_path,"/",get_simulate_id(params_5d,prefix="RW"),sep=""),params_5d)
do_random_walk_simulation(fss_4d,paste(work_path,"/",get_simulate_id(params_4d,prefix="RW"),sep=""),params_4d)
do_random_walk_simulation(fss_3d,paste(work_path,"/",get_simulate_id(params_3d,prefix="RW"),sep=""),params_3d)
do_random_walk_simulation(fss_2d,paste(work_path,"/",get_simulate_id(params_2d,prefix="RW"),sep=""),params_2d)
do_random_walk_simulation(fss_1d,paste(work_path,"/",get_simulate_id(params_1d,prefix="RW"),sep=""),params_1d)
