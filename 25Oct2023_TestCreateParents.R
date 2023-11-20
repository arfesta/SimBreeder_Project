#25 Oct 2023 - simulation using create_parents.R from github repo

# Start by creating genetic map using create_map.R function, check MAF of SNPQTLs
 source('functions/create_map.R')
 source('functions/create_parents_v3.R') # Copied from Github repo on 12 Oct, not tested yet
 source('functions/calc_TGV.R')
 source('functions/sim_phenos.R')
 source('functions/OP_testing.R')
 load("25Sept_seeds.RData")

   genetic.map <- create_genetic_map(num.chromos = 12,map.length = 1800,num.markers = 120,total.QTL = 2952,num.SNPQTL = 1960,
                                    distribute.loci = "even",marker.distribution = "equally-spaced")

 summary(genetic.map$genetic.map$MAF)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.00000 0.01269 0.01647 0.06817 0.03624 0.49970 

 snpqtls <- which(genetic.map$genetic.map$types=="snpqtl")
 summary(genetic.map$genetic.map$MAF[snpqtls])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01001 0.01244 0.01483 0.01490 0.01737 0.02000 

# Change MAFs to be between 0.1 and 0.5
  y_beta <- rbeta(5000, shape1 = 1, shape2 = 5)
  MAFS <- y_beta[intersect(which(y_beta>=0.1),which(y_beta<0.5))]
  genetic.map$genetic.map$MAF[snpqtls] <- sample(MAFS, size = 1960, replace = TRUE)
# create_parents 
  parents <- create_parents(map.info = genetic.map$genetic.map,num.parents = 300,max.delt.allele = 0,par.markers.unique = T,heterozygous.markers = T)
# Check MAFs of SNPQTLs in parents to see if they match the expected 0.1 - 0.5 range
	results.maf <- vector("list")
	for(N in 1:nrow(parents$parent.SNPQTL.matrix)){
		if(length(table(parents$parent.SNPQTL.matrix[N,]))==4){
			results.maf[[N]] <- (2*(table(parents$parent.SNPQTL.matrix[N,]))[1] + 
			table(parents$parent.SNPQTL.matrix[N,])[2])/(2*sum(table(parents$parent.SNPQTL.matrix[N,])))
		} else {
			results.maf[[N]] <- (table(parents$parent.SNPQTL.matrix[N,])[2])/(2*sum(table(parents$parent.SNPQTL.matrix[N,])))
		}
	}
	summary(unlist(results.maf))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03333 0.08167 0.12000 0.14402 0.19167 0.40333		- Not quite as high as hoped, but better than before
	
#calc_TGV 
  parents.TGV <- calc_TGV(geno.info = parents,map.info = genetic.map$genetic.map,A = 0.001,a = -0.001,dom.coeff = 0,founder = T)
# sim_phenos
  parents.phenos <- sim_phenos(TGV.object = parents.TGV,h2 = .3)

  founder.h2 <- var(parents.phenos$genetic.values)/var(parents.phenos$phenos)
  E.sd <- sqrt((var(parents.phenos$genetic.values)/founder.h2) - var(parents.phenos$genetic.values))

  base_pop_data <- list(
   "genetic.map" = genetic.map$genetic.map,
    "parents" = parents,
    "parents.TGV"=parents.TGV,
    "parents.phenos"=parents.phenos,
    "founder.h2"=founder.h2,
    "E.sd"=E.sd)
  s <- all_seeds[1]+10  
  op.test <- OP_testing(map.info=base_pop_data$genetic.map,parent.info=base_pop_data$parents,
                      parent.phenos=base_pop_data$parents.phenos,parents.TGV=base_pop_data$parents.TGV,
                      cross.prog = 1,dom.coeff = 0,A =0.001,a = -0.001,h2 = .3,use.seed = T,starting_seed = s)

# Try running Titan03_ibdblup_NoClust.R code to see if it can finish 10 generations of simulation with this dataset
# First have to source remaining function files
source('functions/create_cross_design.R')
source('functions/make_crosses.R')
source('functions/extract_selections_ibd.R')

# Define variables used in various steps of simulation
## Run Values ####
NumParents=64
prog.per.cross = 60
af.selection = "snpibd"
rel.mat.cross = "snpibd"
num.sel.af = 64
num.selections.within.family = 1
num.of.cores = 30
mt.design="AssortiveMating"

#Values for alleles and trait variation
Major.value =0.001
Minor.value =-0.001
Dominance.Coeff=0
indiv.tree.h2 = .3
pop = 1
this.gen <- all_seeds[[pop]]

# Run first generation of simulation
  ### Create matrices to hold outputs ####
  simulations=1; generations=10
  genetic.gain.mine <- matrix(data=NA,nrow=simulations,ncol=generations)
  genotypic.variance <- matrix(data=NA,nrow=simulations,ncol=generations)
  phenotypic.gain.mine <- matrix(data=NA,nrow=simulations,ncol=generations)
  mean.select.inbreeding<- matrix(data=NA,nrow=simulations,ncol=generations)
  mean.pop.inbreeding<- matrix(data=NA,nrow=simulations,ncol=generations)
  bulmer.effect <- matrix(data=NA,nrow=simulations,ncol=generations)
  coancest.threshold <- matrix(data=NA,nrow=simulations,ncol=generations)
  delt.alleles <- matrix(data=NA,nrow=simulations,ncol=generations)
  max.delt.alleles <- matrix(data=NA,nrow=simulations,ncol=generations)
  min.delt.alleles <- matrix(data=NA,nrow=simulations,ncol=generations)
  selection.snp.genos <- list()
  
  ### Subset the open pollinated parent output to use as input for 1st Generation ####
  op.test$mean.parent.phenos <- op.test$mean.parent.phenos[1:NumParents]
  op.test$mean.parent.tgv <- op.test$mean.parent.tgv[1:NumParents]
  op.test$delt.alleles <- op.test$delt.alleles[1:NumParents]
  op.test$genetic.values <- op.test$genetic.values[1:NumParents]
  op.test$phenos <- op.test$phenos[1:NumParents]
  op.test$genos.3d <- op.test$genos.3d[,1:NumParents,]
  op.test$marker.matrix <- op.test$marker.matrix[1:NumParents,]
  
  #### Create 1st gen progeny#######
  set.seed(this.gen+1)
  cross.design <- create_cross_design(parent.info = op.test,prog.percross = prog.per.cross,generation  = 1,use.op.par.phenos = T,
					coancest.thresh = T)
  set.seed(this.gen+2)
  progeny1 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = base_pop_data$parents,
				num.cores = num.of.cores)
  progeny1.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny1,cross.design = cross.design,A = Major.value,a = Minor.value,
					dom.coeff = Dominance.Coeff,par.markers.unique=T)
  set.seed(this.gen); progeny1.phenos <- sim_phenos(TGV.object = progeny1.TGV,h2 = indiv.tree.h2)
  
  progeny1.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = base_pop_data$parents.TGV,
                                             past.phenos = base_pop_data$parents.phenos,parent.info = op.test,progeny.info = progeny1,
                                             progeny.TGV = progeny1.TGV,progeny.phenos = progeny1.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)

# Completed with no errors - run generation 2 steps
  
  
  genetic.gain.mine[1,1] <- mean(progeny1.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,1] <- var(progeny1.extractions$select.genval)
  phenotypic.gain.mine[1,1] <-mean(progeny1.extractions$selection.phenos)-mean(op.test$phenos)
  mean.select.inbreeding[1,1] <- mean(progeny1.extractions$select.inbred.level)
  mean.pop.inbreeding[1,1]<- mean(progeny1.extractions$prog.inbred.level)
  bulmer.effect[1,1] <- progeny1.extractions$bulmer.effect
  coancest.threshold[1,1] <- 0
  delt.alleles[1,1] <- sum(progeny1.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,1] <- max(progeny1.extractions$delt.allele)
  min.delt.alleles[1,1] <- min(progeny1.extractions$delt.allele)
  selection.snp.genos[1] <- progeny1.extractions$selection.snp.genos
  
  ####Generation 2#####
  set.seed(this.gen*2+1);
  cross.design <- create_cross_design(parent.info = progeny1.extractions,prog.percross = prog.per.cross,generation  = 2,coancest.thresh=T)
  set.seed(this.gen*2+2);progeny2 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny1.extractions,num.cores = num.of.cores)
  progeny2.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny2,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*2); progeny2.phenos <- sim_phenos(TGV.object = progeny2.TGV,h2 = indiv.tree.h2)
  progeny2.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny1.TGV,
                                             past.phenos = progeny1.phenos,parent.info = progeny1.extractions,progeny.info = progeny2,
                                             progeny.TGV = progeny2.TGV,progeny.phenos = progeny2.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = F,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  # Error at this point
  # num.of.cores,useSNPs="subset", useSNPfxn="0.25")
#	Error in length(snpqtl.loci) * useSNPfxn : 
#	  non-numeric argument to binary operator ** SOlution - remove quotes from around 0.25 so it is numeric
  
  genetic.gain.mine[1,2] <- mean(progeny2.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,2]<- var(progeny2.extractions$select.genval)
  phenotypic.gain.mine[1,2] <-mean(progeny2.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,2]<- mean(progeny2.extractions$select.inbred.level)
  mean.pop.inbreeding[1,2]<- mean(progeny2.extractions$prog.inbred.level)
  bulmer.effect[1,2]<- progeny2.extractions$bulmer.effect
  coancest.threshold[1,2] <- cross.design$coancestry.threshold
  delt.alleles[1,2] <- sum(progeny2.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,2] <- max(progeny2.extractions$delt.allele)
  min.delt.alleles[1,2] <- min(progeny2.extractions$delt.allele)
  selection.snp.genos[2] <- progeny2.extractions$selection.snp.genos

  
  ####Generations 3####
  set.seed(this.gen*3+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny2.extractions,prog.percross = prog.per.cross,generation  = 3)
  set.seed(this.gen*3+2); progeny3 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny2.extractions,num.cores = num.of.cores)
  progeny3.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny3,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*3); progeny3.phenos <- sim_phenos(TGV.object = progeny3.TGV,h2 = indiv.tree.h2)
  progeny3.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny2.TGV,
                                             past.phenos = progeny2.phenos,parent.info = progeny2.extractions,progeny.info = progeny3,
                                             progeny.TGV = progeny3.TGV,progeny.phenos = progeny3.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = F,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  
  genetic.gain.mine[1,3] <- mean(progeny3.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,3]<- var(progeny3.extractions$select.genval)
  phenotypic.gain.mine[1,3] <-mean(progeny3.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,3]<- mean(progeny3.extractions$select.inbred.level)
  mean.pop.inbreeding[1,3]<- mean(progeny3.extractions$prog.inbred.level)
  bulmer.effect[1,3]<- progeny3.extractions$bulmer.effect
  coancest.threshold[1,3] <- cross.design$coancestry.threshold
  delt.alleles[1,3] <- sum(progeny3.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,3] <- max(progeny3.extractions$delt.allele)
  min.delt.alleles[1,3] <- min(progeny3.extractions$delt.allele)
  selection.snp.genos[3] <- progeny3.extractions$selection.snp.genos

  
  ####Generations 3####
  set.seed(this.gen*3+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny2.extractions,prog.percross = prog.per.cross,generation  = 3)
  set.seed(this.gen*3+2); progeny3 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny2.extractions,num.cores = num.of.cores)
  progeny3.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny3,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*3); progeny3.phenos <- sim_phenos(TGV.object = progeny3.TGV,h2 = indiv.tree.h2)
  progeny3.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny2.TGV,
                                             past.phenos = progeny2.phenos,parent.info = progeny2.extractions,progeny.info = progeny3,
                                             progeny.TGV = progeny3.TGV,progeny.phenos = progeny3.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = F,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  
  genetic.gain.mine[1,3] <- mean(progeny3.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,3]<- var(progeny3.extractions$select.genval)
  phenotypic.gain.mine[1,3] <-mean(progeny3.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,3]<- mean(progeny3.extractions$select.inbred.level)
  mean.pop.inbreeding[1,3]<- mean(progeny3.extractions$prog.inbred.level)
  bulmer.effect[1,3]<- progeny3.extractions$bulmer.effect
  coancest.threshold[1,3] <- cross.design$coancestry.threshold
  delt.alleles[1,3] <- sum(progeny3.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,3] <- max(progeny3.extractions$delt.allele)
  min.delt.alleles[1,3] <- min(progeny3.extractions$delt.allele)
  selection.snp.genos[3] <- progeny3.extractions$selection.snp.genos
 
  ####Generations 3####
  set.seed(this.gen*3+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny2.extractions,prog.percross = prog.per.cross,generation  = 3)
  set.seed(this.gen*3+2); progeny3 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny2.extractions,num.cores = num.of.cores)
  progeny3.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny3,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*3); progeny3.phenos <- sim_phenos(TGV.object = progeny3.TGV,h2 = indiv.tree.h2)
  progeny3.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny2.TGV,
                                             past.phenos = progeny2.phenos,parent.info = progeny2.extractions,progeny.info = progeny3,
                                             progeny.TGV = progeny3.TGV,progeny.phenos = progeny3.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = F,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  
  genetic.gain.mine[1,3] <- mean(progeny3.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,3]<- var(progeny3.extractions$select.genval)
  phenotypic.gain.mine[1,3] <-mean(progeny3.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,3]<- mean(progeny3.extractions$select.inbred.level)
  mean.pop.inbreeding[1,3]<- mean(progeny3.extractions$prog.inbred.level)
  bulmer.effect[1,3]<- progeny3.extractions$bulmer.effect
  coancest.threshold[1,3] <- cross.design$coancestry.threshold
  delt.alleles[1,3] <- sum(progeny3.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,3] <- max(progeny3.extractions$delt.allele)
  min.delt.alleles[1,3] <- min(progeny3.extractions$delt.allele)
  selection.snp.genos[3] <- progeny3.extractions$selection.snp.genos

 
  ####Generation 4####
  set.seed(this.gen*4+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny3.extractions,prog.percross = prog.per.cross,generation  = 4)
  set.seed(this.gen*4+2); progeny4 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny3.extractions,num.cores = num.of.cores)
  progeny4.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny4,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*4); progeny4.phenos <- sim_phenos(TGV.object = progeny4.TGV,h2 = indiv.tree.h2)
  progeny4.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny3.TGV,
                                             past.phenos = progeny3.phenos,parent.info = progeny3.extractions,progeny.info = progeny4,
                                             progeny.TGV = progeny4.TGV,progeny.phenos = progeny4.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  genetic.gain.mine[1,4] <- mean(progeny4.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,4]<- var(progeny4.extractions$select.genval)
  phenotypic.gain.mine[1,4] <-mean(progeny4.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,4]<- mean(progeny4.extractions$select.inbred.level)
  mean.pop.inbreeding[1,4]<- mean(progeny4.extractions$prog.inbred.level)
  bulmer.effect[1,4]<- progeny4.extractions$bulmer.effect
  coancest.threshold[1,4] <- cross.design$coancestry.threshold
  delt.alleles[1,4] <- sum(progeny4.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,4] <- max(progeny4.extractions$delt.allele)
  min.delt.alleles[1,4] <- min(progeny4.extractions$delt.allele)
  selection.snp.genos[4] <- progeny4.extractions$selection.snp.genos
  
  
  ####Generation 5####
  set.seed(this.gen*5+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny4.extractions,prog.percross = prog.per.cross,generation  = 5)
  set.seed(this.gen*5+2); progeny5 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny4.extractions,num.cores = num.of.cores)
  progeny5.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny5,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*5); progeny5.phenos <- sim_phenos(TGV.object = progeny5.TGV,h2 = indiv.tree.h2)
  progeny5.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny4.TGV,
                                             past.phenos = progeny4.phenos,parent.info = progeny4.extractions,progeny.info = progeny5,
                                             progeny.TGV = progeny5.TGV,progeny.phenos = progeny5.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  genetic.gain.mine[1,5] <- mean(progeny5.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,5]<- var(progeny5.extractions$select.genval)
  phenotypic.gain.mine[1,5] <-mean(progeny5.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,5]<- mean(progeny5.extractions$select.inbred.level)
  mean.pop.inbreeding[1,5]<- mean(progeny5.extractions$prog.inbred.level)
  bulmer.effect[1,5]<- progeny5.extractions$bulmer.effect
  coancest.threshold[1,5] <- cross.design$coancestry.threshold
  delt.alleles[1,5] <- sum(progeny5.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,5] <- max(progeny5.extractions$delt.allele)
  min.delt.alleles[1,5] <- min(progeny5.extractions$delt.allele)
  selection.snp.genos[5] <- progeny5.extractions$selection.snp.genos
  
  ####Generation 6####
  set.seed(this.gen*6+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny5.extractions,prog.percross = prog.per.cross,generation  = 6)
  set.seed(this.gen*6+2); progeny6 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny5.extractions,num.cores = num.of.cores)
  progeny6.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny6,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*6); progeny6.phenos <- sim_phenos(TGV.object = progeny6.TGV,h2 = indiv.tree.h2)
  progeny6.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny5.TGV,
                                             past.phenos = progeny5.phenos,parent.info = progeny5.extractions,progeny.info = progeny6,
                                             progeny.TGV = progeny6.TGV,progeny.phenos = progeny6.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  genetic.gain.mine[1,6] <- mean(progeny6.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,6]<- var(progeny6.extractions$select.genval)
  phenotypic.gain.mine[1,6] <-mean(progeny6.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,6]<- mean(progeny6.extractions$select.inbred.level)
  mean.pop.inbreeding[1,6]<- mean(progeny6.extractions$prog.inbred.level)
  bulmer.effect[1,6]<- progeny6.extractions$bulmer.effect
  coancest.threshold[1,6] <- cross.design$coancestry.threshold
  delt.alleles[1,6] <- sum(progeny6.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,6] <- max(progeny6.extractions$delt.allele)
  min.delt.alleles[1,6] <- min(progeny6.extractions$delt.allele)
  selection.snp.genos[6] <- progeny6.extractions$selection.snp.genos
  
  ####Generation 7####
  set.seed(this.gen*7+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny6.extractions,prog.percross = prog.per.cross,generation  = 7)
  set.seed(this.gen*7+2); progeny7 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny6.extractions,num.cores = num.of.cores)
  progeny7.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny7,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*7); progeny7.phenos <- sim_phenos(TGV.object = progeny7.TGV,h2 = indiv.tree.h2)
  progeny7.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny6.TGV,
                                             past.phenos = progeny6.phenos,parent.info = progeny6.extractions,progeny.info = progeny7,
                                             progeny.TGV = progeny7.TGV,progeny.phenos = progeny7.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  genetic.gain.mine[1,7] <- mean(progeny7.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,7]<- var(progeny7.extractions$select.genval)
  phenotypic.gain.mine[1,7] <-mean(progeny7.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,7]<- mean(progeny7.extractions$select.inbred.level)
  mean.pop.inbreeding[1,7]<- mean(progeny7.extractions$prog.inbred.level)
  bulmer.effect[1,7]<- progeny7.extractions$bulmer.effect
  coancest.threshold[1,7] <- cross.design$coancestry.threshold
  delt.alleles[1,7] <- sum(progeny7.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,7] <- max(progeny7.extractions$delt.allele)
  min.delt.alleles[1,7] <- min(progeny7.extractions$delt.allele)
  selection.snp.genos[7] <- progeny7.extractions$selection.snp.genos
  
  ####Generation 8####
  set.seed(this.gen*8+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny7.extractions,prog.percross = prog.per.cross,generation  = 8)
  set.seed(this.gen*8+2); progeny8 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny7.extractions,num.cores = num.of.cores)
  progeny8.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny8,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*8); progeny8.phenos <- sim_phenos(TGV.object = progeny8.TGV,h2 = indiv.tree.h2)
  progeny8.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny7.TGV,
                                             past.phenos = progeny7.phenos,parent.info = progeny7.extractions,progeny.info = progeny8,
                                             progeny.TGV = progeny8.TGV,progeny.phenos = progeny8.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  genetic.gain.mine[1,8] <- mean(progeny8.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,8]<- var(progeny8.extractions$select.genval)
  phenotypic.gain.mine[1,8] <-mean(progeny8.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,8]<- mean(progeny8.extractions$select.inbred.level)
  mean.pop.inbreeding[1,8]<- mean(progeny8.extractions$prog.inbred.level)
  bulmer.effect[1,8]<- progeny8.extractions$bulmer.effect
  coancest.threshold[1,8] <- cross.design$coancestry.threshold
  delt.alleles[1,8] <- sum(progeny8.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,8] <- max(progeny8.extractions$delt.allele)
  min.delt.alleles[1,8] <- min(progeny8.extractions$delt.allele)
  selection.snp.genos[8] <- progeny8.extractions$selection.snp.genos
  
  ####Generation 9####
  set.seed(this.gen*9+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny8.extractions,prog.percross = prog.per.cross,generation  = 9)
  set.seed(this.gen*9+2); progeny9 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny8.extractions,num.cores = num.of.cores)
  progeny9.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny9,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*9); progeny9.phenos <- sim_phenos(TGV.object = progeny9.TGV,h2 = indiv.tree.h2)
  progeny9.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny8.TGV,
                                             past.phenos = progeny8.phenos,parent.info = progeny8.extractions,progeny.info = progeny9,
                                             progeny.TGV = progeny9.TGV,progeny.phenos = progeny9.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  genetic.gain.mine[1,9] <- mean(progeny9.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,9]<- var(progeny9.extractions$select.genval)
  phenotypic.gain.mine[1,9] <-mean(progeny9.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,9]<- mean(progeny9.extractions$select.inbred.level)
  mean.pop.inbreeding[1,9]<- mean(progeny9.extractions$prog.inbred.level)
  bulmer.effect[1,9]<- progeny9.extractions$bulmer.effect
  coancest.threshold[1,9] <- cross.design$coancestry.threshold
  delt.alleles[1,9] <- sum(progeny9.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,9] <- max(progeny9.extractions$delt.allele)
  min.delt.alleles[1,9] <- min(progeny9.extractions$delt.allele)
  selection.snp.genos[9] <- progeny9.extractions$selection.snp.genos
    
  ####Generation 10####
  set.seed(this.gen*10+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny9.extractions,prog.percross = prog.per.cross,generation  = 10)
  set.seed(this.gen*10+2); progeny10 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny9.extractions,num.cores = num.of.cores)
  progeny10.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny10,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*10); progeny10.phenos <- sim_phenos(TGV.object = progeny10.TGV,h2 = indiv.tree.h2)
  progeny10.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny9.TGV,
                                              past.phenos = progeny9.phenos,parent.info = progeny9.extractions,progeny.info = progeny10,
                                              progeny.TGV = progeny10.TGV,progeny.phenos = progeny10.phenos,among.family.selection = af.selection,
                                              num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores,useSNPs="subset", useSNPfxn=0.25)
  
  genetic.gain.mine[1,10] <- mean(progeny10.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,10]<- var(progeny10.extractions$select.genval)
  phenotypic.gain.mine[1,10] <-mean(progeny10.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,10]<- mean(progeny10.extractions$select.inbred.level)
  mean.pop.inbreeding[1,10]<- mean(progeny10.extractions$prog.inbred.level)
  bulmer.effect[1,10]<- progeny10.extractions$bulmer.effect
  coancest.threshold[1,10] <- cross.design$coancestry.threshold
  delt.alleles[1,10] <- sum(progeny10.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,10] <- max(progeny10.extractions$delt.allele)
  min.delt.alleles[1,10] <- min(progeny10.extractions$delt.allele)
  selection.snp.genos[10] <- progeny10.extractions$selection.snp.genos
  
  # All 10 generations completed successfully - copy results to list
   ibd_rand25pct_List <- list("g"=genetic.gain.mine,"gv"=genotypic.variance, "pg"=phenotypic.gain.mine,"si"=mean.select.inbreeding,
                  "pi"=mean.pop.inbreeding,"be"=bulmer.effect, "ct"=coancest.threshold, "maxd"=max.delt.alleles,"mind"=min.delt.alleles,
                  "da"=delt.alleles,"snp.genos"=selection.snp.genos)
 
 # Now run 10 generations of ABLUP simulation
 ## Run Values ####
NumParents=64
prog.per.cross = 60
af.selection = "ABLUP"
rel.mat.cross = "pedigree"
num.sel.af = 64
num.of.cores = 1
mt.design="AssortiveMating"

#Values for alleles and trait variation
Major.value =0.001
Minor.value =-0.001
Dominance.Coeff=0
indiv.tree.h2 = .3
  ###Create matrices to hold outputs####
  simulations=1; generations=10
  genetic.gain.mine <- matrix(data=NA,nrow=simulations,ncol=generations)
  genotypic.variance <- matrix(data=NA,nrow=simulations,ncol=generations)
  phenotypic.gain.mine <- matrix(data=NA,nrow=simulations,ncol=generations)
  mean.select.inbreeding<- matrix(data=NA,nrow=simulations,ncol=generations)
  mean.pop.inbreeding<- matrix(data=NA,nrow=simulations,ncol=generations)
  bulmer.effect <- matrix(data=NA,nrow=simulations,ncol=generations)
  coancest.threshold <- matrix(data=NA,nrow=simulations,ncol=generations)
  delt.alleles <- matrix(data=NA,nrow=simulations,ncol=generations)
  max.delt.alleles <- matrix(data=NA,nrow=simulations,ncol=generations)
  min.delt.alleles <- matrix(data=NA,nrow=simulations,ncol=generations)
  ### Create list to hold selected progeny snpqtl genotype matrices from even-numbered generations
  selection.snp.genos <- list()
  
  # Don't repeat OP testing - use same OP-selected parents as for ibd simulation
    #### Create 1st gen progeny#######
  set.seed(this.gen+1)
  cross.design <- create_cross_design(parent.info = op.test,prog.percross = prog.per.cross,generation  = 1,use.op.par.phenos = T,coancest.thresh = T)
  set.seed(this.gen+2)
  progeny1 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = base_pop_data$parents,num.cores = num.of.cores)
  progeny1.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny1,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen); progeny1.phenos <- sim_phenos(TGV.object = progeny1.TGV,h2 = indiv.tree.h2)
  
  progeny1.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = base_pop_data$parents.TGV,
                                             past.phenos = base_pop_data$parents.phenos,parent.info = op.test,progeny.info = progeny1,
                                             progeny.TGV = progeny1.TGV,progeny.phenos = progeny1.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  
  genetic.gain.mine[1,1] <- mean(progeny1.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,1] <- var(progeny1.extractions$select.genval)
  phenotypic.gain.mine[1,1] <-mean(progeny1.extractions$selection.phenos)-mean(op.test$phenos)
  mean.select.inbreeding[1,1] <- mean(progeny1.extractions$select.inbred.level)
  mean.pop.inbreeding[1,1]<- mean(progeny1.extractions$prog.inbred.level)
  bulmer.effect[1,1] <- progeny1.extractions$bulmer.effect
  coancest.threshold[1,1] <- 0
  delt.alleles[1,1] <- sum(progeny1.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,1] <- max(progeny1.extractions$delt.allele)
  min.delt.alleles[1,1] <- min(progeny1.extractions$delt.allele)
  selection.snp.genos[1] <- progeny1.extractions$selection.snp.genos

  rm(progeny1); gc()
  
  ####Generation 2#####
  set.seed(this.gen*2+1);
  cross.design <- create_cross_design(parent.info = progeny1.extractions,prog.percross = prog.per.cross,generation  = 2,coancest.thresh=T)
  set.seed(this.gen*2+2);progeny2 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny1.extractions,num.cores = num.of.cores)
  progeny2.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny2,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*2); progeny2.phenos <- sim_phenos(TGV.object = progeny2.TGV,h2 = indiv.tree.h2)
  progeny2.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny1.TGV,
                                             past.phenos = progeny1.phenos,parent.info = progeny1.extractions,progeny.info = progeny2,
                                             progeny.TGV = progeny2.TGV,progeny.phenos = progeny2.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  
  genetic.gain.mine[1,2] <- mean(progeny2.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,2]<- var(progeny2.extractions$select.genval)
  phenotypic.gain.mine[1,2] <-mean(progeny2.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,2]<- mean(progeny2.extractions$select.inbred.level)
  mean.pop.inbreeding[1,2]<- mean(progeny2.extractions$prog.inbred.level)
  bulmer.effect[1,2]<- progeny2.extractions$bulmer.effect
  coancest.threshold[1,2] <- cross.design$coancestry.threshold
  delt.alleles[1,2] <- sum(progeny2.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,2] <- max(progeny2.extractions$delt.allele)
  min.delt.alleles[1,2] <- min(progeny2.extractions$delt.allele)
  selection.snp.genos[2] <- progeny2.extractions$selection.snp.genos

  rm(progeny1.TGV,progeny1.extractions, progeny1.phenos,progeny2);gc()
  
  ####Generations 3####
  set.seed(this.gen*3+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny2.extractions,prog.percross = prog.per.cross,generation  = 3)
  set.seed(this.gen*3+2); progeny3 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny2.extractions,num.cores = num.of.cores)
  progeny3.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny3,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*3); progeny3.phenos <- sim_phenos(TGV.object = progeny3.TGV,h2 = indiv.tree.h2)
  progeny3.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny2.TGV,
                                             past.phenos = progeny2.phenos,parent.info = progeny2.extractions,progeny.info = progeny3,
                                             progeny.TGV = progeny3.TGV,progeny.phenos = progeny3.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  
  genetic.gain.mine[1,3] <- mean(progeny3.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,3]<- var(progeny3.extractions$select.genval)
  phenotypic.gain.mine[1,3] <-mean(progeny3.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,3]<- mean(progeny3.extractions$select.inbred.level)
  mean.pop.inbreeding[1,3]<- mean(progeny3.extractions$prog.inbred.level)
  bulmer.effect[1,3]<- progeny3.extractions$bulmer.effect
  coancest.threshold[1,3] <- cross.design$coancestry.threshold
  delt.alleles[1,3] <- sum(progeny3.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,3] <- max(progeny3.extractions$delt.allele)
  min.delt.alleles[1,3] <- min(progeny3.extractions$delt.allele)
  selection.snp.genos[3] <- progeny3.extractions$selection.snp.genos
  
  rm(progeny2.TGV,progeny2.extractions, progeny2.phenos,progeny3);gc()
  
  ####Generations 4####
  set.seed(this.gen*4+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny3.extractions,prog.percross = prog.per.cross,generation  = 4)
  set.seed(this.gen*4+2); progeny4 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny3.extractions,num.cores = num.of.cores)
  progeny4.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny4,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*4); progeny4.phenos <- sim_phenos(TGV.object = progeny4.TGV,h2 = indiv.tree.h2)
  progeny4.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny3.TGV,
                                             past.phenos = progeny3.phenos,parent.info = progeny3.extractions,progeny.info = progeny4,
                                             progeny.TGV = progeny4.TGV,progeny.phenos = progeny4.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  genetic.gain.mine[1,4] <- mean(progeny4.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,4]<- var(progeny4.extractions$select.genval)
  phenotypic.gain.mine[1,4] <-mean(progeny4.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,4]<- mean(progeny4.extractions$select.inbred.level)
  mean.pop.inbreeding[1,4]<- mean(progeny4.extractions$prog.inbred.level)
  bulmer.effect[1,4]<- progeny4.extractions$bulmer.effect
  coancest.threshold[1,4] <- cross.design$coancestry.threshold
  delt.alleles[1,4] <- sum(progeny4.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,4] <- max(progeny4.extractions$delt.allele)
  min.delt.alleles[1,4] <- min(progeny4.extractions$delt.allele)
  selection.snp.genos[4] <- progeny4.extractions$selection.snp.genos
  
  rm(progeny3.TGV,progeny3.extractions, progeny3.phenos,progeny4);gc()
  
  
  ####Generations 5####
  set.seed(this.gen*5+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny4.extractions,prog.percross = prog.per.cross,generation  = 5)
  set.seed(this.gen*5+2); progeny5 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny4.extractions,num.cores = num.of.cores)
  progeny5.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny5,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*5); progeny5.phenos <- sim_phenos(TGV.object = progeny5.TGV,h2 = indiv.tree.h2)
  progeny5.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny4.TGV,
                                             past.phenos = progeny4.phenos,parent.info = progeny4.extractions,progeny.info = progeny5,
                                             progeny.TGV = progeny5.TGV,progeny.phenos = progeny5.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  genetic.gain.mine[1,5] <- mean(progeny5.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,5]<- var(progeny5.extractions$select.genval)
  phenotypic.gain.mine[1,5] <-mean(progeny5.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,5]<- mean(progeny5.extractions$select.inbred.level)
  mean.pop.inbreeding[1,5]<- mean(progeny5.extractions$prog.inbred.level)
  bulmer.effect[1,5]<- progeny5.extractions$bulmer.effect
  coancest.threshold[1,5] <- cross.design$coancestry.threshold
  delt.alleles[1,5] <- sum(progeny5.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,5] <- max(progeny5.extractions$delt.allele)
  min.delt.alleles[1,5] <- min(progeny5.extractions$delt.allele)
  selection.snp.genos[5] <- progeny5.extractions$selection.snp.genos
  
  rm(progeny4.TGV,progeny4.extractions, progeny4.phenos,progeny5);gc()
  
  ####Generations 6####
  set.seed(this.gen*6+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny5.extractions,prog.percross = prog.per.cross,generation  = 6)
  set.seed(this.gen*6+2); progeny6 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny5.extractions,num.cores = num.of.cores)
  progeny6.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny6,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*6); progeny6.phenos <- sim_phenos(TGV.object = progeny6.TGV,h2 = indiv.tree.h2)
  progeny6.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny5.TGV,
                                             past.phenos = progeny5.phenos,parent.info = progeny5.extractions,progeny.info = progeny6,
                                             progeny.TGV = progeny6.TGV,progeny.phenos = progeny6.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  genetic.gain.mine[1,6] <- mean(progeny6.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,6]<- var(progeny6.extractions$select.genval)
  phenotypic.gain.mine[1,6] <-mean(progeny6.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,6]<- mean(progeny6.extractions$select.inbred.level)
  mean.pop.inbreeding[1,6]<- mean(progeny6.extractions$prog.inbred.level)
  bulmer.effect[1,6]<- progeny6.extractions$bulmer.effect
  coancest.threshold[1,6] <- cross.design$coancestry.threshold
  delt.alleles[1,6] <- sum(progeny6.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,6] <- max(progeny6.extractions$delt.allele)
  min.delt.alleles[1,6] <- min(progeny6.extractions$delt.allele)
  selection.snp.genos[6] <- progeny6.extractions$selection.snp.genos
  
  rm(progeny5.TGV,progeny5.extractions, progeny5.phenos,progeny6);gc()
  
  ####Generations 7####
  set.seed(this.gen*7+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny6.extractions,prog.percross = prog.per.cross,generation  = 7)
  set.seed(this.gen*7+2); progeny7 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny6.extractions,num.cores = num.of.cores)
  progeny7.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny7,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*7); progeny7.phenos <- sim_phenos(TGV.object = progeny7.TGV,h2 = indiv.tree.h2)
  progeny7.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny6.TGV,
                                             past.phenos = progeny6.phenos,parent.info = progeny6.extractions,progeny.info = progeny7,
                                             progeny.TGV = progeny7.TGV,progeny.phenos = progeny7.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  genetic.gain.mine[1,7] <- mean(progeny7.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,7]<- var(progeny7.extractions$select.genval)
  phenotypic.gain.mine[1,7] <-mean(progeny7.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,7]<- mean(progeny7.extractions$select.inbred.level)
  mean.pop.inbreeding[1,7]<- mean(progeny7.extractions$prog.inbred.level)
  bulmer.effect[1,7]<- progeny7.extractions$bulmer.effect
  coancest.threshold[1,7] <- cross.design$coancestry.threshold
  delt.alleles[1,7] <- sum(progeny7.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,7] <- max(progeny7.extractions$delt.allele)
  min.delt.alleles[1,7] <- min(progeny7.extractions$delt.allele)
  selection.snp.genos[7] <- progeny7.extractions$selection.snp.genos
  
  rm(progeny6.TGV,progeny6.extractions, progeny6.phenos,progeny7);gc()
  
  ####Generations 8####
  set.seed(this.gen*8+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny7.extractions,prog.percross = prog.per.cross,generation  = 8)
  set.seed(this.gen*8+2); progeny8 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny7.extractions,num.cores = num.of.cores)
  progeny8.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny8,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*8); progeny8.phenos <- sim_phenos(TGV.object = progeny8.TGV,h2 = indiv.tree.h2)
  progeny8.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny7.TGV,
                                             past.phenos = progeny7.phenos,parent.info = progeny7.extractions,progeny.info = progeny8,
                                             progeny.TGV = progeny8.TGV,progeny.phenos = progeny8.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  genetic.gain.mine[1,8] <- mean(progeny8.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,8]<- var(progeny8.extractions$select.genval)
  phenotypic.gain.mine[1,8] <-mean(progeny8.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,8]<- mean(progeny8.extractions$select.inbred.level)
  mean.pop.inbreeding[1,8]<- mean(progeny8.extractions$prog.inbred.level)
  bulmer.effect[1,8]<- progeny8.extractions$bulmer.effect
  coancest.threshold[1,8] <- cross.design$coancestry.threshold
  delt.alleles[1,8] <- sum(progeny8.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,8] <- max(progeny8.extractions$delt.allele)
  min.delt.alleles[1,8] <- min(progeny8.extractions$delt.allele)
  selection.snp.genos[8] <- progeny8.extractions$selection.snp.genos
  
  rm(progeny7.TGV,progeny7.extractions, progeny7.phenos,progeny8);gc()
  
  ####Generations 9####
  set.seed(this.gen*9+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny8.extractions,prog.percross = prog.per.cross,generation  = 9)
  set.seed(this.gen*9+2); progeny9 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny8.extractions,num.cores = num.of.cores)
  progeny9.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny9,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*9); progeny9.phenos <- sim_phenos(TGV.object = progeny9.TGV,h2 = indiv.tree.h2)
  progeny9.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny8.TGV,
                                             past.phenos = progeny8.phenos,parent.info = progeny8.extractions,progeny.info = progeny9,
                                             progeny.TGV = progeny9.TGV,progeny.phenos = progeny9.phenos,among.family.selection = af.selection,
                                             num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  genetic.gain.mine[1,9] <- mean(progeny9.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,9]<- var(progeny9.extractions$select.genval)
  phenotypic.gain.mine[1,9] <-mean(progeny9.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,9]<- mean(progeny9.extractions$select.inbred.level)
  mean.pop.inbreeding[1,9]<- mean(progeny9.extractions$prog.inbred.level)
  bulmer.effect[1,9]<- progeny9.extractions$bulmer.effect
  coancest.threshold[1,9] <- cross.design$coancestry.threshold
  delt.alleles[1,9] <- sum(progeny9.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,9] <- max(progeny9.extractions$delt.allele)
  min.delt.alleles[1,9] <- min(progeny9.extractions$delt.allele)
  selection.snp.genos[9] <- progeny9.extractions$selection.snp.genos
  
  rm(progeny8.TGV,progeny8.extractions, progeny8.phenos,progeny9);gc()
  
  ####Generations 10####
  set.seed(this.gen*10+1); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny9.extractions,prog.percross = prog.per.cross,generation  = 10)
  set.seed(this.gen*10+2); progeny10 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny9.extractions,num.cores = num.of.cores)
  progeny10.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny10,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
  set.seed(this.gen*10); progeny10.phenos <- sim_phenos(TGV.object = progeny10.TGV,h2 = indiv.tree.h2)
  progeny10.extractions <- extract_selections_ibd(map.info = base_pop_data$genetic.map,cross.design = cross.design,past.tgv = progeny9.TGV,
                                              past.phenos = progeny9.phenos,parent.info = progeny9.extractions,progeny.info = progeny10,
                                              progeny.TGV = progeny10.TGV,progeny.phenos = progeny10.phenos,among.family.selection = af.selection,
                                              num.selections.among.family = num.sel.af,reduced = T,num.cores = num.of.cores)
  
  genetic.gain.mine[1,10] <- mean(progeny10.extractions$select.genval) - mean(op.test$genetic.values)
  genotypic.variance[1,10]<- var(progeny10.extractions$select.genval)
  phenotypic.gain.mine[1,10] <-mean(progeny10.extractions$selection.phenos) - mean(op.test$phenos)
  mean.select.inbreeding[1,10]<- mean(progeny10.extractions$select.inbred.level)
  mean.pop.inbreeding[1,10]<- mean(progeny10.extractions$prog.inbred.level)
  bulmer.effect[1,10]<- progeny10.extractions$bulmer.effect
  coancest.threshold[1,10] <- cross.design$coancestry.threshold
  delt.alleles[1,10] <- sum(progeny10.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,10] <- max(progeny10.extractions$delt.allele)
  min.delt.alleles[1,10] <- min(progeny10.extractions$delt.allele)
  selection.snp.genos[10] <- progeny10.extractions$selection.snp.genos
  
  rm(progeny9.TGV,progeny9.extractions, progeny9.phenos,progeny10,progeny10.extractions,progeny10.phenos,progeny10.TGV);gc()
  
  
  #### List for output objects####
  ablup_List <- list("g"=genetic.gain.mine,"gv"=genotypic.variance, "pg"=phenotypic.gain.mine,"si"=mean.select.inbreeding,
                  "pi"=mean.pop.inbreeding,"be"=bulmer.effect, "ct"=coancest.threshold, "maxd"=max.delt.alleles,"mind"=min.delt.alleles,
                  "da"=delt.alleles,"snp.genos"=selection.snp.genos)

# Recover the results from ablup_List and ibd_rand25pct_List and plot genetic gain over 10 generations
results <- matrix(NA,nrow=1,ncol=12)
results.df <- as.data.frame(results)
names(results.df) <- c("select","result","gen01","gen02","gen03","gen04","gen05","gen06","gen07","gen08","gen09","gen10")
results.df <- rbind(results.df,c("ablup","ggain",round(ablup_List$g,0))) 
results.df <- rbind(results.df,c("ibdblup","ggain",round(ibd_rand25pct_List$g,0)))                 
  
