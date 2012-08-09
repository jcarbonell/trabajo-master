
args <- commandArgs(trailingOnly = TRUE)

home_path <- args[1]

source(paste(home_path,"/","shortestpath-mariangela.R",sep=""))
source(paste(home_path,"/","score.r",sep=""))
source(paste(home_path,"/","utils.r",sep=""))
source(paste(home_path,"/","postprocess.r",sep=""))
require(WriteXLS)

#
# INPUT PARAMS
#

gene_list_files <- unlist(strsplit(args[2],","))
family_names <- unlist(strsplit(args[3],","))
outdir <- args[4]

#
# PREPARE DATA
#

dir.create(outdir)

raw_gene_list <- load_family_genes(gene_list_files=gene_list_files,group_names=family_names,hasHeader=T)

interactomes <- load_interactomes(paste(home_path,"/interactomes/",sep=""))

global_score_methods <- c("no_zero_mean","no_zero_max")


#
# COMPUTE SCORE
#

multi_score <- compute_multi_score(raw_gene_list,interactomes,default_score,global_score_methods=global_score_methods,verbose=T)


#
# OUTPUT
#
header <- paste("gene",paste(colnames(multi_score$global_score_frame),collapse="\t"),sep="\t")
header <- gsub("no_zero_mean","global_mean",header)
header <- gsub("no_zero_max","global_max",header)
header <- gsub("gene_cluster_score","cluster_score",header)

global_score_frame_file <- paste(outdir,"/global_score_frame.txt",sep="")

write(header,file=global_score_frame_file)
write.table(format(digits=5,multi_score$global_score_frame),file=global_score_frame_file,sep="\t",quote=F,col.names=F,append=T)
