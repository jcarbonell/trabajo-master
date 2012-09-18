source("param_chart_utils.r")

fams <- c(3,4,5,7,10)

#
# SHORTEST PATH
#

sp_files <- c(
  "tabla_familias/SP__s100__f3__d3__g100__09_09_2012__10:32:25__/sre_sp.rdata",
  "tabla_familias/SP__s100__f4__d4__g100__09_09_2012__10:57:03__/sre_sp.rdata",
  "tabla_familias/SP__s100__f5__d5__g100__09_09_2012__11:44:16__/sre_sp.rdata",
  "tabla_familias/SP__s100__f7__d7__g100__09_09_2012__13:05:14__/sre_sp.rdata",
  "tabla_familias/SP__s100__f10__d10__g100__09_09_2012__16:12:20__/sre_sp.rdata"  
)
sp_quantiles <- get_quantiles(sp_files)
sp_data <- do.call("rbind",lapply(1:length(sp_quantiles$quantiles), function(x) cbind(sp_quantiles$quantiles[[x]],fams[x],1)))



#
# INTERMEDIATION
#

inter_files <- c(
  "tabla_familias/INTER__s100__f3__d3__g100__14_09_2012__20:48:50__/sre_inter.rdata",
  "tabla_familias/INTER__s100__f4__d4__g100__14_09_2012__21:03:34__/sre_inter.rdata",
  "tabla_familias/INTER__s100__f5__d5__g100__14_09_2012__21:34:00__/sre_inter.rdata",
  "tabla_familias/INTER__s100__f7__d7__g100__14_09_2012__22:29:40__/sre_inter.rdata",
  "tabla_familias/INTER__s100__f10__d10__g100__15_09_2012__00:53:33__/sre_inter.rdata"  
)
inter_quantiles <- get_quantiles(inter_files)
inter_data <- do.call("rbind",lapply(1:length(inter_quantiles$quantiles), function(x) cbind(inter_quantiles$quantiles[[x]],fams[x],2)))


#
# RANDOM WALK
#


rw_files <- c(
  "tabla_familias/RW__s100__f3__d3__g100__10_09_2012__00:22:10__/sre_rw.rdata",
  "tabla_familias/RW__s100__f4__d4__g100__10_09_2012__00:37:44__/sre_rw.rdata",
  "tabla_familias/RW__s100__f5__d5__g100__10_09_2012__01:09:13__/sre_rw.rdata",
  "tabla_familias/RW__s100__f7__d7__g100__10_09_2012__02:06:05__/sre_rw.rdata",
  "tabla_familias/RW__s100__f10__d10__g100__10_09_2012__04:27:53__/sre_rw.rdata"
)
rw_quantiles <- get_quantiles(rw_files)
rw_data <- do.call("rbind",lapply(1:length(rw_quantiles$quantiles), function(x) cbind(rw_quantiles$quantiles[[x]],fams[x],3)))


#
# ALL DATA TOGETHER
#

# get data
tabla_familias_datos <- as.data.frame(rbind(sp_data,inter_data,rw_data))
colnames(tabla_familias_datos) <- c("rank","values","distancia")
tabla_familias_datos$distancia <- factor(tabla_familias_datos$distancia)
levels(tabla_familias_datos$distancia) <- list(SP=1,ID=2,RW=3)

# save data
save(tabla_familias_datos,file="tabla_familias_datos.rdata")






