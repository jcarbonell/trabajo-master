
rm(list=ls())

source("igraph_utils.R")
source("nh_simulation_lib.R")
source("shortestpath-mariangela.R")
source("nh_lib.R")


# load interactome
interactome <- read.table(file="reference_lists/translated_hsa_alldb_curatedphysical_genes_interactome_nr.sif", header = FALSE, sep = "\t", quote = "",stringsAsFactors = FALSE, as.is=TRUE)

# prepare family configuration
families <- c("fam1","fam2","fam3")
# overlapping <- matrix(c( 10,3,0, 3,10,0, 0,0,15 ),ncol=3,byrow=T); rownames(overlapping) <- families; colnames(overlapping) <- families
# overlapping <- matrix(c( 5,2,0, 2,4,0, 0,0,7 ),ncol=3,byrow=T); rownames(overlapping) <- families; colnames(overlapping) <- families
overlapping <- matrix(c( 2,0,0, 0,1,0, 0,0,2 ),ncol=3,byrow=T); rownames(overlapping) <- families; colnames(overlapping) <- families

# simulate genes
raw_gene_list <- nh_simulate_genes(families,overlapping,interactome,family_proportion=0.4)
gene_list <- unique(unlist(raw_gene_list))

# recompute interactome
ninterm <- 5
subnet <- get.all.shortest.paths.Luz(interactome, gene_list, ninterm)

# compute score
scores <- nh_score(raw_gene_list,subnet)
clean_subnet <- nh_remove_inter_inter_interactions(subnet,scores$hits)

# paint network
layout <- nh_get_standard_layout(subnet)
paint_family_graph(subnet,scores$pnh,scores$hits,title="Simulated genes",paint_labels=T,layout=layout)
paint_family_graph(clean_subnet,scores$pnh,scores$hits,title="Simulated genes",paint_labels=T,layout=layout)

# selected neighbors
# neighborhood <- nh_get_neiborhood("RB1",subnet,radio=1)
# paint_family_graph(neighborhood,scores$pnh,scores$hits,title="RB1 neighborhood",paint_labels=T,layout=layout)
# 
# best_nfc <- names(which(scores$nfc==max(scores$nfc)))[1]
# best_nfc_nhood <- nh_get_neiborhood(best_nfc,subnet,radio=1)
# paint_family_graph(best_nfc_nhood,scores$pnh,scores$hits,title=paste(best_nfc,"neighborhood (best NFC)"),paint_labels=T,layout=layout)
# 
# best_pnh <- names(which(scores$pnh==max(scores$pnh)))[1]
# best_pnh_nhood <- nh_get_neiborhood(best_pnh,subnet,radio=2)
# paint_family_graph(best_pnh_nhood,scores$pnh,scores$hits,highlighteds=best_pnh,title=paste(best_pnh,"neighborhood (best PNH)"),paint_labels=T,layout=layout)

# rank
#nh_classif <- nh_kmeans_classify(scores$pnh)
nh_classif <- nh_quantile_classify(scores$pnh,0.95)

nh_paint_classification(nh_classif)

selected_family_genes <- gene_list[which(is.element(gene_list,nh_classif$selected_genes))]

for(node in selected_family_genes){
  nhood <- nh_get_neiborhood(node,clean_subnet,radio=3)
  paint_family_graph(nhood,scores$pnh,scores$hits,highlighteds=node,title=paste(node,"neighborhood"),paint_labels=T,layout=layout)
  scan()
}

