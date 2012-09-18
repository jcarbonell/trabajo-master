source("param_chart_utils.r")

values <- c(1,2,3,4,5)

#
# SHORTEST PATH
#

sp_files <- c(
  "tabla_disease/SP__s100__f5__d1__g100__09_09_2012__15:58:21__/sre_sp.rdata",
  "tabla_disease/SP__s100__f5__d2__g100__09_09_2012__14:39:12__/sre_sp.rdata",
  "tabla_disease/SP__s100__f5__d3__g100__09_09_2012__13:18:35__/sre_sp.rdata",
  "tabla_disease/SP__s100__f5__d4__g100__09_09_2012__11:57:22__/sre_sp.rdata",
  "tabla_disease/SP__s100__f5__d5__g100__09_09_2012__10:38:14__/sre_sp.rdata"
)
sp_quantiles <- get_quantiles(sp_files)
sp_data <- do.call("rbind",lapply(1:length(sp_quantiles$quantiles), function(x) cbind(sp_quantiles$quantiles[[x]],values[x],1)))



#
# INTERMEDIATION
#


inter_files <- c(
  "tabla_disease/INTER__s100__f5__d1__g100__15_09_2012__00:30:03__/sre_inter.rdata",
  "tabla_disease/INTER__s100__f5__d2__g100__14_09_2012__23:35:15__/sre_inter.rdata",
  "tabla_disease/INTER__s100__f5__d3__g100__14_09_2012__22:39:06__/sre_inter.rdata",
  "tabla_disease/INTER__s100__f5__d4__g100__14_09_2012__21:45:07__/sre_inter.rdata",
  "tabla_disease/INTER__s100__f5__d5__g100__14_09_2012__20:49:20__/sre_inter.rdata"
)
inter_quantiles <- get_quantiles(inter_files)
inter_data <- do.call("rbind",lapply(1:length(inter_quantiles$quantiles), function(x) cbind(inter_quantiles$quantiles[[x]],values[x],2)))


#
# RANDOM WALK
#


rw_files <- c(
  "tabla_disease/RW__s100__f5__d1__g100__09_09_2012__20:53:49__/sre_rw.rdata",
  "tabla_disease/RW__s100__f5__d2__g100__09_09_2012__19:59:41__/sre_rw.rdata",
  "tabla_disease/RW__s100__f5__d3__g100__09_09_2012__19:04:52__/sre_rw.rdata",
  "tabla_disease/RW__s100__f5__d4__g100__09_09_2012__18:09:36__/sre_rw.rdata",
  "tabla_disease/RW__s100__f5__d5__g100__09_09_2012__17:14:49__/sre_rw.rdata"
)
rw_quantiles <- get_quantiles(rw_files)
rw_data <- do.call("rbind",lapply(1:length(rw_quantiles$quantiles), function(x) cbind(rw_quantiles$quantiles[[x]],values[x],3)))


#
# ALL DATA TOGETHER
#

# get data
tabla_disease_datos <- as.data.frame(rbind(sp_data,inter_data,rw_data))
colnames(tabla_disease_datos) <- c("rank","values","distancia")
tabla_disease_datos$distancia <- factor(tabla_disease_datos$distancia)
levels(tabla_disease_datos$distancia) <- list(SP=1,ID=2,RW=3)

# save data
save(tabla_disease_datos,file="tabla_disease_datos.rdata")






