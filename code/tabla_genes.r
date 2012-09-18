
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
work_path <- paste(home_path,"/work/tabla_genes",sep="")
dir.create(work_path)

# family set series
params_25g <- list(nsims=100,nfams=3,ndisease_genes=3,ngenes=25)
params_50g <- list(nsims=100,nfams=3,ndisease_genes=3,ngenes=50)
params_100g <- list(nsims=100,nfams=3,ndisease_genes=3,ngenes=100)
params_200g <- list(nsims=100,nfams=3,ndisease_genes=3,ngenes=200)
params_500g <- list(nsims=100,nfams=3,ndisease_genes=3,ngenes=500)
fss_25g <- simulate_known_family_set_series(params=params_25g)
fss_50g <- simulate_known_family_set_series(params=params_50g)
fss_100g <- simulate_known_family_set_series(params=params_100g)
fss_200g <- simulate_known_family_set_series(params=params_200g)
fss_500g <- simulate_known_family_set_series(params=params_500g)

# shortest_path
do_shortest_path_simulation(fss_25g,paste(work_path,"/",get_simulate_id(params_25g,prefix="SP"),sep=""),params_25g)
do_shortest_path_simulation(fss_50g,paste(work_path,"/",get_simulate_id(params_50g,prefix="SP"),sep=""),params_50g)
do_shortest_path_simulation(fss_100g,paste(work_path,"/",get_simulate_id(params_100g,prefix="SP"),sep=""),params_100g)
do_shortest_path_simulation(fss_200g,paste(work_path,"/",get_simulate_id(params_200g,prefix="SP"),sep=""),params_200g)
do_shortest_path_simulation(fss_500g,paste(work_path,"/",get_simulate_id(params_500g,prefix="SP"),sep=""),params_500g)

# # intermediation
# do_intermediation_simulation(fss_25g,paste(work_path,"/",get_simulate_id(params_25g,prefix="INTER"),sep=""),params_25g)
# do_intermediation_simulation(fss_50g,paste(work_path,"/",get_simulate_id(params_50g,prefix="INTER"),sep=""),params_50g)
# do_intermediation_simulation(fss_100g,paste(work_path,"/",get_simulate_id(params_100g,prefix="INTER"),sep=""),params_100g)
# do_intermediation_simulation(fss_200g,paste(work_path,"/",get_simulate_id(params_200g,prefix="INTER"),sep=""),params_200g)
# do_intermediation_simulation(fss_500g,paste(work_path,"/",get_simulate_id(params_500g,prefix="INTER"),sep=""),params_500g)

# random walk
do_random_walk_simulation(fss_25g,paste(work_path,"/",get_simulate_id(params_25g,prefix="RW"),sep=""),params_25g)
do_random_walk_simulation(fss_50g,paste(work_path,"/",get_simulate_id(params_50g,prefix="RW"),sep=""),params_50g)
do_random_walk_simulation(fss_100g,paste(work_path,"/",get_simulate_id(params_100g,prefix="RW"),sep=""),params_100g)
do_random_walk_simulation(fss_200g,paste(work_path,"/",get_simulate_id(params_200g,prefix="RW"),sep=""),params_200g)
do_random_walk_simulation(fss_500g,paste(work_path,"/",get_simulate_id(params_500g,prefix="RW"),sep=""),params_500g)
