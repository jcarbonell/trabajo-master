
do_shortest_path_simulation <- function(family_set_series,outdir,params,get_back=F){
  
  print_title("SHORTEST PATH")
  
  if(file.exists(outdir)){
    
    cat("ERROR: output folder already exists!!!!")
    
  } else{
      
    dir.create(outdir)
    
    # load interactomes
    interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
    
    # run prioritization
    sre_sp <- run_and_evaluate(family_set_series,interactomes,score_function_names=c("final_score"),global_score_methods=c("no_zero_mean","no_zero_max"),test.inter=F)

    # save
    sre_sp$params <- params
    sre_sp$family_set_series <- family_set_series
    save(sre_sp,file=paste(outdir,"/sre_sp.rdata",sep=""))
    
    # plot
    png(paste(outdir,"/summary.png",sep=""),width=700,height=500)
    paint_global_score(sre_sp$score_runs[[1]],label="Shortest path")
    dev.off()    
    
  }
    
  gc()
  
  if(get_back){
    return(sre_sp)
  }
  
}



do_intermediation_simulation <- function(family_set_series,outdir,params,get_back=F){
  
  print_title("INTERMEDIATION")
  
  if(file.exists(outdir)){
    
    cat("ERROR: output folder already exists!!!!")
    
  } else{
  
    dir.create(outdir)
    
    # load interactomes
    interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
    load("../interactomes/rwd/intermediation_matrices.rdata")
        
    # run prioritization
    sre_inter <- run_and_evaluate(family_set_series,interactomes,probability_matrices=intermediation_matrices,score_function_names=c("final_score"),global_score_methods=c("no_zero_mean","no_zero_max"),test.inter=F)

    # save
    sre_inter$params <- params
    sre_inter$family_set_series <- family_set_series
    save(sre_inter,file=paste(outdir,"/sre_inter.rdata",sep=""))
    
    # plot
    png(paste(outdir,"/summary.png",sep=""),width=700,height=500)
    paint_global_score(sre_inter$score_runs[[1]],label="Intermediation")
    dev.off()
    
  }
    
  gc()
  
  if(get_back){
    return(sre_inter)
  }
  
}

do_random_walk_simulation <- function(family_set_series,outdir,params,get_back=F){
  
  print_title("RANDOM WALK")
    
  if(file.exists(outdir)){
    
    cat("ERROR: output folder already exists!!!!")
    
  } else{
  
    dir.create(outdir)
    
    # load interactomes
    interactomes <- load_interactomes(paste(home_path,"/../interactomes/",sep=""))
    load("../interactomes/rwd/all_rwd.rdata")
    interactomes <- interactomes[names(probability_matrices)]
    
    # run prioritization
    sre_rw <- run_and_evaluate(family_set_series,interactomes,probability_matrices=probability_matrices,score_function_names=c("final_score"),global_score_methods=c("no_zero_mean","no_zero_max"),test.inter=F)

    # save
    sre_rw$params <- params
    sre_rw$family_set_series <- family_set_series
    save(sre_rw,file=paste(outdir,"/sre_rw.rdata",sep=""))
    
    # plot
    png(paste(outdir,"/summary.png",sep=""),width=700,height=500)
    paint_global_score(sre_rw$score_runs[[1]],label="Random walk")
    dev.off()
    
  }
    
  gc()
  
  if(get_back){
    return(sre_inter)
  }
  
}

print_title <- function(thetitle){
  cat("\n>>>> ",thetitle," simulation run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
}

get_simulate_id <- function(params,prefix="sim",sufix=""){  
  paste(prefix,"__s",params$nsims,"__f",params$nfams,"__d",params$ndisease_genes,"__g",params$ngenes,"__",format(Sys.time(), "%d_%m_%Y__%H:%M:%S"),"__",sufix,sep="")  
}