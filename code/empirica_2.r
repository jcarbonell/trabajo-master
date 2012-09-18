source("lib/empirica_utils.r")

all_data <- read.table("~/Pensacola/master/trabajo_master_backup/work2/alldata.txt",stringsAsFactors=F)$V1

n <- length(all_data)
distance <- substr(basename(dirname(all_data)),1,2)
for(i in 1:n){
  cat("Processing file ",i," de ", n, "...\n",sep="") 
  create_score_data(all_data[i],distance[i],output=T)
  gc()
}
cat("finished")