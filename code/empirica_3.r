
source("lib/empirica_utils.r")
require(ggplot2)

all_experiments <- read.table("charts_data/all_experiments.txt",sep="\t",header=T,stringsAsFactors=F)


# global by distance
png("charts_data/emp_tabla_global.png",width=700,heigh=400)
ggplot(all_experiments, aes(no_zero_mean, colour = factor(distance))) + geom_density() + geom_histogram() + opts(title="Distribución global")
dev.off()

# tabla RW varianción en genes y nfams
tabla_familias_rw <- filter_data(all_experiments,distance="RW",sim_set="tabla_familias")
tabla_genes_rw <- filter_data(all_experiments,distance="RW",sim_set="tabla_genes")

p1 <- ggplot(tabla_familias_rw, aes(no_zero_mean, colour = factor(nfams))) +  geom_density() + opts(title="Random walk") +
  scale_x_log10() + ylab("densidad") + 
  opts(plot.margin = unit(1*c(2, 2, 2, 2), "lines")) +
  opts(axis.title.x=theme_text(vjust=-1,size=13)) + 
  opts(axis.title.y=theme_text(angle=90, vjust=-0.2,size=13)) + 
  opts(plot.title=theme_text(size=16, vjust=3))
p2 <- ggplot(tabla_genes_rw, aes(no_zero_mean, colour = factor(ngenes))) +  geom_density() +
  scale_x_log10() + ylab("densidad") + 
  opts(plot.margin = unit(1*c(0, 2, 4, 2), "lines")) +
  opts(axis.title.x=theme_text(vjust=-1,size=13)) + 
  opts(axis.title.y=theme_text(angle=90, vjust=-0.1,size=13)) + 
  opts(plot.title=theme_text(size=16, vjust=3))
png("charts_data/emp_rw_nfams_ngenes",width=700,heigh=500)
multiplot(p1,p2)
dev.off()