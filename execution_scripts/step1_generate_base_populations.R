# Run this on Domino ####
source("/mnt/results/functions/create_map.R")
source("/mnt/results/functions/create_parents.R")
source("/mnt/results/functions/calc_TGV.R")
source("/mnt/results/functions/sim_phenos.R")
set.seed(2323523)
seed.vector <- sample(1235:4536345,size = 300,replace =F)
save(seed.vector,file="/mnt/results/3k/all.seeds.RData")
rm(.Random.seed);library(parallel)
mclapply(1:30,function(pop){
  set.seed(seed.vector[pop])
  genetic.map <- create_genetic_map(num.chromos = 12,map.length = 1800,snp.markers = 3000,fullyinformed_markers =  120,total.QTL = 2640,
                                    total.SNPQTL = 1960,snp.marker.maf = c(.2,.49),snp.qtl.maf = c(.01,.49),shape_param = c(1,5),
                                    distribute.loci = "even",snp.qtl.distribution = "random",snp.marker.distribution = "even")
  
  parents <- create_parents(map.info = genetic.map,num.parents=280,max.delt.allele=14)
  parents.TGV <- calc_TGV(  geno.info = parents,map.info = genetic.map,A = 1,a = -100,dom.coeff = 1,founder = T)
  parents.phenos <- sim_phenos(TGV.object = parents.TGV,h2 = .3)
  founder.h2 <- var(parents.phenos$genetic.values)/var(parents.phenos$phenos)
  E.sd <- sqrt((var(parents.phenos$genetic.values)/founder.h2) - var(parents.phenos$genetic.values))
  
  this.pop <- paste0("/mnt/results/3k/base_populations/","base_population_",pop,".RData")
  base_pop_data <- list(
    "genetic.map" = genetic.map,
    "parents" = parents,
    "parents.TGV"=parents.TGV,
    "parents.phenos"=parents.phenos,
    "founder.h2"=founder.h2,
    "E.sd"=E.sd)
  save(base_pop_data,file = this.pop,compress = T)
},mc.cores=5)

