# Now for each base population conduct OP testing ####
source('/mnt/results/functions/OP_testing.R')
load("/mnt/results/3k/all.seeds.RData")
lapply(1:30,function(pop){
  s <- (seed.vector[pop])+100
  load.pop.name <- paste0("/mnt/results/3k/base_populations/base_population_",pop,".RData")
  load(load.pop.name)
  op.test <- OP_testing(map.info=base_pop_data$genetic.map,parent.info=base_pop_data$parents,
                        parent.phenos=base_pop_data$parents.phenos,parents.TGV=base_pop_data$parents.TGV,
                        cross.prog = 1,dom.coeff = 1,A =1,a = -100,h2 = .3,n.cores= 5,use.seed = T,starting_seed = s)
  
  save.pop.name <- paste0("/mnt/results/3k/base_populations_op_test/base_population_op_test_",pop,".RData")
  save(base_pop_data,op.test,file=save.pop.name,compress=T)
})
