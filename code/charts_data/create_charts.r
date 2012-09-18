
source("~/")
load("tabla_distancia/SP__s100__f3__d3__g150__09_09_2012__10:22:42__/sre_sp.rdata")
load("tabla_distancia/INTER__s100__f3__d3__g150__14_09_2012__20:44:17__/sre_inter.rdata")
load("tabla_distancia/RW__s100__f3__d3__g150__09_09_2012__12:18:58__/sre_rw.rdata")



require(ggplot2)
source("param_chart_utils.r")

load("tabla_familias_datos.rdata")
png("tabla_familias.png",width=700,height=500)
param_distribution(tabla_familias_datos,"número de familias","Número de familias")
dev.off()

load("tabla_disease_datos.rdata")
png("tabla_disease.png",width=700,height=500)
param_distribution(tabla_disease_datos,"número de genes causantes","Número de genes causantes")
dev.off()
