#Source file for breeding simulator####
source('/media/titan/ARF/Simulator.Package/breeding.simulator.R')

#Create a list of 1000 different sets of 10 seeds for each base population####
the.seed <- vector("list")
first <- 1000
for(i in 1:1000){
  the.seed[[i]] <- c(seq(first,first+9,1))
  first <- first + 9
}

#Inputs for simulations####
  NumParents=64
  prog.per.cross = 60
  af.selection = "ABLUP"
  wf.selection = "NO"
  sel.strategy = "within.family"
  rel.mat.cross = "pedigree"
  num.sel.af = 64
  run.in.parallel = T
  num.of.cores = 20

#Values for alleles and trait variation
  Major.value =1
  Minor.value =-100
  Dominance.Coeff=1
  indiv.tree.h2 = .3

#For each base population run the simulator for 10 generations #####
for(pop.num in 1:300){ #pop.num specifies the base population to use
###Load base population####
  pop.name <- paste("/media/titan/ARF/base.population/base.population.",pop.num,".RData",sep="")
  load(pop.name)

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
  
  
###Subset the open pollinated parent output to use as input for 1st Generation####
  OP.parent.phenos$pars <- OP.parent.phenos$pars[1:NumParents]
  OP.parent.phenos$delt.alleles <- OP.parent.phenos$delt.alleles[1:NumParents]
  OP.parent.phenos$genetic.values <- OP.parent.phenos$genetic.values[1:NumParents]
  OP.parent.phenos$phenos <- OP.parent.phenos$phenos[1:NumParents]
  OP.parent.phenos$genos.3d <- OP.parent.phenos$genos.3d[,1:NumParents,]
  OP.parent.phenos$marker.matrix <- OP.parent.phenos$marker.matrix[1:NumParents,]
  
  
  #64 parents ABLUP among, ABLUP within family
#### Create 1st gen progeny#######
  set.seed(the.seed[[pop.num]][1])
  cross.design <- create.cross.design(parentinfo = OP.parent.phenos,prog.percross = prog.per.cross,gen = 1,use.op.par.phenos = T)
  progeny1 <- make.crosses(crossdesign = cross.design, parentinfo = parents, run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny1.TGV <- create.progeny.TGV(crossdesign = cross.design, proginfo = progeny1)
  set.seed(NULL); progeny1.phenos <- create.progeny.phenos(crossdesign = cross.design, progeny.TGV = progeny1.TGV, h2 = indiv.tree.h2)
  progeny1.extractions <- extract.selections(crossdesign = cross.design, parentinfo = OP.parent.phenos, past.tgv = parents.TGV,past.phenos = parents.phenos, relmatrix = rel.mat.cross,
                                             progenyinfo = progeny1,progenyTGV = progeny1.TGV,progenyphenos = progeny1.phenos,selection.strategy = sel.strategy,
                                             among.family.selection = af.selection, within.family.selection = wf.selection, reduced = T,numSelections.among.family = num.sel.af)
  
  genetic.gain.mine[1,1] <- mean(progeny1.extractions$select.genval)-founder.TGV
  genotypic.variance[1,1] <- var(progeny1.extractions$select.genval)
  phenotypic.gain.mine[1,1] <-mean(progeny1.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,1] <- mean(progeny1.extractions$select.inbred.level)
  mean.pop.inbreeding[1,1]<- mean(progeny1.extractions$prog.inbred.level)
  bulmer.effect[1,1] <- progeny1.extractions$bulmer.effect
  coancest.threshold[1,1] <- 0
  delt.alleles[1,1] <- sum(progeny1.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,1] <- max(progeny1.extractions$delt.allele)
  min.delt.alleles[1,1] <- min(progeny1.extractions$delt.allele)
  rm(progeny1); gc()
  
####Generation 2#####
  set.seed(the.seed[[pop.num]][2])
  cross.design <- create.cross.design(parentinfo = progeny1.extractions,prog.percross = prog.per.cross,gen = 2)
  progeny2 <- make.crosses(crossdesign=cross.design,parentinfo = progeny1.extractions,run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny2.TGV <- create.progeny.TGV(crossdesign=cross.design,proginfo = progeny2)
  set.seed(NULL); progeny2.phenos <- create.progeny.phenos(crossdesign=cross.design,progeny.TGV = progeny2.TGV, h2 = indiv.tree.h2)
  progeny2.extractions <- extract.selections(crossdesign=cross.design, past.tgv = progeny1.TGV, relmatrix = rel.mat.cross,
                                             past.phenos = progeny1.phenos, parentinfo = progeny1.extractions,progenyinfo = progeny2,
                                             progenyTGV = progeny2.TGV,progenyphenos = progeny2.phenos,selection.strategy = sel.strategy,
                                             among.family.selection = af.selection,within.family.selection = wf.selection,numSelections.among.family = num.sel.af)
  
  genetic.gain.mine[1,2] <- mean(progeny2.extractions$select.genval)-founder.TGV
  genotypic.variance[1,2]<- var(progeny2.extractions$select.genval)
  phenotypic.gain.mine[1,2] <-mean(progeny2.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,2]<- mean(progeny2.extractions$select.inbred.level)
  mean.pop.inbreeding[1,2]<- mean(progeny2.extractions$prog.inbred.level)
  bulmer.effect[1,2]<- progeny2.extractions$bulmer.effect
  coancest.threshold[1,2] <- cross.design$coancestry.threshold
  delt.alleles[1,2] <- sum(progeny2.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,2] <- max(progeny2.extractions$delt.allele)
  min.delt.alleles[1,2] <- min(progeny2.extractions$delt.allele)
  
  
  rm(progeny1.TGV,progeny1.extractions, progeny1.phenos,progeny2);gc()
####Generations 3####
  set.seed(the.seed[[pop.num]][3])
  cross.design <- create.cross.design(parentinfo = progeny2.extractions,prog.percross = prog.per.cross,gen = 3)
  progeny3 <- make.crosses(crossdesign=cross.design,parentinfo = progeny2.extractions,run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny3.TGV <- create.progeny.TGV(crossdesign=cross.design,proginfo = progeny3)
  set.seed(NULL); progeny3.phenos <- create.progeny.phenos(crossdesign=cross.design,progeny.TGV = progeny3.TGV,h2 = indiv.tree.h2)
  progeny3.extractions <- extract.selections(crossdesign=cross.design, past.tgv = progeny2.TGV, relmatrix = rel.mat.cross,
                                             past.phenos = progeny2.phenos, parentinfo = progeny2.extractions,progenyinfo = progeny3,
                                             progenyTGV = progeny3.TGV,progenyphenos = progeny3.phenos,selection.strategy = sel.strategy, 
                                             among.family.selection = af.selection,within.family.selection = wf.selection,numSelections.among.family = num.sel.af)
  
  
  genetic.gain.mine[1,3] <- mean(progeny3.extractions$select.genval)-founder.TGV
  genotypic.variance[1,3]<- var(progeny3.extractions$select.genval)
  phenotypic.gain.mine[1,3] <-mean(progeny3.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,3]<- mean(progeny3.extractions$select.inbred.level)
  mean.pop.inbreeding[1,3]<- mean(progeny3.extractions$prog.inbred.level)
  bulmer.effect[1,3]<- progeny3.extractions$bulmer.effect
  coancest.threshold[1,3] <- cross.design$coancestry.threshold
  delt.alleles[1,3] <- sum(progeny3.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,3] <- max(progeny3.extractions$delt.allele)
  min.delt.alleles[1,3] <- min(progeny3.extractions$delt.allele)
  
  rm(progeny2.TGV,progeny2.extractions, progeny2.phenos,progeny3); gc()
  
####Generation 4####
  set.seed(the.seed[[pop.num]][4])
  cross.design <- create.cross.design(parentinfo = progeny3.extractions,prog.percross = prog.per.cross,gen = 4)
  progeny4 <- make.crosses(crossdesign=cross.design,parentinfo = progeny3.extractions,run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny4.TGV <- create.progeny.TGV(crossdesign=cross.design,proginfo = progeny4)
  set.seed(NULL); progeny4.phenos <- create.progeny.phenos(crossdesign=cross.design,progeny.TGV = progeny4.TGV, h2 = indiv.tree.h2)
  progeny4.extractions <- extract.selections(crossdesign=cross.design, relmatrix = rel.mat.cross,
                                             past.tgv = progeny3.TGV, past.phenos = progeny3.phenos, parentinfo = progeny3.extractions,
                                             progenyinfo = progeny4, progenyTGV = progeny4.TGV,progenyphenos = progeny4.phenos,
                                             selection.strategy = sel.strategy, 
                                             among.family.selection = af.selection,within.family.selection = wf.selection,numSelections.among.family = num.sel.af)
  
  
  genetic.gain.mine[1,4] <- mean(progeny4.extractions$select.genval)-founder.TGV
  genotypic.variance[1,4]<- var(progeny4.extractions$select.genval)
  phenotypic.gain.mine[1,4] <-mean(progeny4.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,4]<- mean(progeny4.extractions$select.inbred.level)
  mean.pop.inbreeding[1,4]<- mean(progeny4.extractions$prog.inbred.level)
  bulmer.effect[1,4]<- progeny4.extractions$bulmer.effect
  coancest.threshold[1,4] <- cross.design$coancestry.threshold
  delt.alleles[1,4] <- sum(progeny4.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,4] <- max(progeny4.extractions$delt.allele)
  min.delt.alleles[1,4] <- min(progeny4.extractions$delt.allele)
  
  rm(progeny3.TGV,progeny3.extractions, progeny3.phenos,progeny4)
  
####Generation 5####
  set.seed(the.seed[[pop.num]][5])
  cross.design <- create.cross.design(parentinfo = progeny4.extractions,prog.percross = prog.per.cross,gen = 5)
  progeny5 <- make.crosses(crossdesign=cross.design,parentinfo = progeny4.extractions,run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny5.TGV <- create.progeny.TGV(crossdesign=cross.design,proginfo = progeny5)
  set.seed(NULL); progeny5.phenos <- create.progeny.phenos(crossdesign=cross.design,progeny.TGV = progeny5.TGV, h2 = indiv.tree.h2)
  progeny5.extractions <- extract.selections(crossdesign=cross.design, relmatrix = rel.mat.cross,
                                             past.tgv = progeny4.TGV, past.phenos = progeny4.phenos, parentinfo = progeny4.extractions,
                                             progenyinfo = progeny5, progenyTGV = progeny5.TGV,progenyphenos = progeny5.phenos,
                                             selection.strategy = sel.strategy, 
                                             among.family.selection = af.selection,within.family.selection = wf.selection,numSelections.among.family = num.sel.af)
  
  
  genetic.gain.mine[1,5] <- mean(progeny5.extractions$select.genval)-founder.TGV
  genotypic.variance[1,5]<- var(progeny5.extractions$select.genval)
  phenotypic.gain.mine[1,5] <-mean(progeny5.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,5]<- mean(progeny5.extractions$select.inbred.level)
  mean.pop.inbreeding[1,5]<- mean(progeny5.extractions$prog.inbred.level)
  bulmer.effect[1,5]<- progeny5.extractions$bulmer.effect
  coancest.threshold[1,5] <- cross.design$coancestry.threshold
  delt.alleles[1,5] <- sum(progeny5.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,5] <- max(progeny5.extractions$delt.allele)
  min.delt.alleles[1,5] <- min(progeny5.extractions$delt.allele)
  
  rm(progeny4.TGV,progeny4.extractions, progeny4.phenos,progeny5); gc()
  
####Generation 6####
  set.seed(the.seed[[pop.num]][6])
  cross.design <- create.cross.design(parentinfo = progeny5.extractions,prog.percross = prog.per.cross,gen = 6)
  progeny6 <- make.crosses(crossdesign=cross.design,parentinfo = progeny5.extractions,run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny6.TGV <- create.progeny.TGV(crossdesign=cross.design,proginfo = progeny6)
  set.seed(NULL); progeny6.phenos <- create.progeny.phenos(crossdesign=cross.design,progeny.TGV = progeny6.TGV, h2 = indiv.tree.h2 )
  progeny6.extractions <- extract.selections(crossdesign=cross.design,relmatrix = rel.mat.cross,
                                             past.tgv = progeny5.TGV, past.phenos = progeny5.phenos, parentinfo = progeny5.extractions,
                                             progenyinfo = progeny6, progenyTGV = progeny6.TGV,progenyphenos = progeny6.phenos,
                                             selection.strategy = sel.strategy, 
                                             among.family.selection = af.selection,within.family.selection = wf.selection,numSelections.among.family = num.sel.af)
  
  
  genetic.gain.mine[1,6] <- mean(progeny6.extractions$select.genval)-founder.TGV
  genotypic.variance[1,6]<- var(progeny6.extractions$select.genval)
  phenotypic.gain.mine[1,6] <-mean(progeny6.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,6]<- mean(progeny6.extractions$select.inbred.level)
  mean.pop.inbreeding[1,6]<- mean(progeny6.extractions$prog.inbred.level)
  bulmer.effect[1,6]<- progeny6.extractions$bulmer.effect
  coancest.threshold[1,6] <- cross.design$coancestry.threshold
  delt.alleles[1,6] <- sum(progeny6.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,6] <- max(progeny6.extractions$delt.allele)
  min.delt.alleles[1,6] <- min(progeny6.extractions$delt.allele)
  
  rm(progeny5.TGV,progeny5.extractions, progeny5.phenos,progeny6); gc()
  
####Generation 7####
  set.seed(the.seed[[pop.num]][7])
  cross.design <- create.cross.design(parentinfo = progeny6.extractions,prog.percross = prog.per.cross,gen = 7)
  progeny7 <- make.crosses(crossdesign=cross.design,parentinfo = progeny6.extractions,run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny7.TGV <- create.progeny.TGV(crossdesign=cross.design,proginfo = progeny7)
  set.seed(NULL); progeny7.phenos <- create.progeny.phenos(crossdesign=cross.design,progeny.TGV = progeny7.TGV, h2 = indiv.tree.h2 )
  progeny7.extractions <- extract.selections(crossdesign=cross.design,relmatrix = rel.mat.cross,
                                             past.tgv = progeny6.TGV, past.phenos = progeny6.phenos, parentinfo = progeny6.extractions,
                                             progenyinfo = progeny7, progenyTGV = progeny7.TGV,progenyphenos = progeny7.phenos,
                                             selection.strategy = sel.strategy, 
                                             among.family.selection = af.selection,within.family.selection = wf.selection,numSelections.among.family = num.sel.af)
  
  
  genetic.gain.mine[1,7] <- mean(progeny7.extractions$select.genval)-founder.TGV
  genotypic.variance[1,7]<- var(progeny7.extractions$select.genval)
  phenotypic.gain.mine[1,7] <-mean(progeny7.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,7]<- mean(progeny7.extractions$select.inbred.level)
  mean.pop.inbreeding[1,7]<- mean(progeny7.extractions$prog.inbred.level)
  bulmer.effect[1,7]<- progeny7.extractions$bulmer.effect
  coancest.threshold[1,7] <- cross.design$coancestry.threshold
  delt.alleles[1,7] <- sum(progeny7.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,7] <- max(progeny7.extractions$delt.allele)
  min.delt.alleles[1,7] <- min(progeny7.extractions$delt.allele)
  
  rm(progeny6.TGV,progeny6.extractions, progeny6.phenos,progeny7); gc()
  
####Generation 8####
  set.seed(the.seed[[pop.num]][8])
  cross.design <- create.cross.design(parentinfo = progeny7.extractions,prog.percross = prog.per.cross,gen = 8)
  progeny8 <- make.crosses(crossdesign=cross.design,parentinfo = progeny7.extractions,run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny8.TGV <- create.progeny.TGV(crossdesign=cross.design,proginfo = progeny8)
  set.seed(NULL); progeny8.phenos <- create.progeny.phenos(crossdesign=cross.design,progeny.TGV = progeny8.TGV, h2 = indiv.tree.h2)
  progeny8.extractions <- extract.selections(crossdesign=cross.design,relmatrix = rel.mat.cross,
                                             past.tgv = progeny7.TGV, past.phenos = progeny7.phenos, parentinfo = progeny7.extractions,
                                             progenyinfo = progeny8, progenyTGV = progeny8.TGV,progenyphenos = progeny8.phenos,
                                             selection.strategy = sel.strategy, 
                                             among.family.selection = af.selection,within.family.selection = wf.selection,numSelections.among.family = num.sel.af)
  
  
  genetic.gain.mine[1,8] <- mean(progeny8.extractions$select.genval)-founder.TGV
  genotypic.variance[1,8]<- var(progeny8.extractions$select.genval)
  phenotypic.gain.mine[1,8] <-mean(progeny8.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,8]<- mean(progeny8.extractions$select.inbred.level)
  mean.pop.inbreeding[1,8]<- mean(progeny8.extractions$prog.inbred.level)
  bulmer.effect[1,8]<- progeny8.extractions$bulmer.effect
  coancest.threshold[1,8] <- cross.design$coancestry.threshold
  delt.alleles[1,8] <- sum(progeny8.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,8] <- max(progeny8.extractions$delt.allele)
  min.delt.alleles[1,8] <- min(progeny8.extractions$delt.allele)
  
  rm(progeny7.TGV,progeny7.extractions, progeny7.phenos,progeny8); gc()
  
####Generation 9####
  set.seed(the.seed[[pop.num]][9])
  cross.design <- create.cross.design(parentinfo = progeny8.extractions,prog.percross = prog.per.cross,gen = 9)
  progeny9 <- make.crosses(crossdesign=cross.design,parentinfo = progeny8.extractions,run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny9.TGV <- create.progeny.TGV(crossdesign=cross.design,proginfo = progeny9)
  set.seed(NULL); progeny9.phenos <- create.progeny.phenos(crossdesign=cross.design,progeny.TGV = progeny9.TGV, h2 = indiv.tree.h2)
  progeny9.extractions <- extract.selections(crossdesign=cross.design,relmatrix = rel.mat.cross,
                                             past.tgv = progeny8.TGV, past.phenos = progeny8.phenos, parentinfo = progeny8.extractions,
                                             progenyinfo = progeny9, progenyTGV = progeny9.TGV,progenyphenos = progeny9.phenos,
                                             selection.strategy = sel.strategy, 
                                             among.family.selection = af.selection,within.family.selection = wf.selection,numSelections.among.family = num.sel.af)
  
  
  genetic.gain.mine[1,9] <- mean(progeny9.extractions$select.genval)-founder.TGV
  genotypic.variance[1,9]<- var(progeny9.extractions$select.genval)
  phenotypic.gain.mine[1,9] <-mean(progeny9.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,9]<- mean(progeny9.extractions$select.inbred.level)
  mean.pop.inbreeding[1,9]<- mean(progeny9.extractions$prog.inbred.level)
  bulmer.effect[1,9]<- progeny9.extractions$bulmer.effect
  coancest.threshold[1,9] <- cross.design$coancestry.threshold
  delt.alleles[1,9] <- sum(progeny9.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,9] <- max(progeny9.extractions$delt.allele)
  min.delt.alleles[1,9] <- min(progeny9.extractions$delt.allele)
  rm(progeny8.TGV,progeny8.extractions, progeny8.phenos,progeny9); gc()
  
####Generation 10####
  set.seed(the.seed[[pop.num]][10])
  cross.design <- create.cross.design(parentinfo = progeny9.extractions,prog.percross = prog.per.cross,gen = 10)
  progeny10 <- make.crosses(crossdesign=cross.design,parentinfo = progeny9.extractions,run.parallel = run.in.parallel,num.cores = num.of.cores)
  progeny10.TGV <- create.progeny.TGV(crossdesign=cross.design,proginfo = progeny10)
  set.seed(NULL); progeny10.phenos <- create.progeny.phenos(crossdesign=cross.design,progeny.TGV = progeny10.TGV, h2 = indiv.tree.h2 )
  progeny10.extractions <- extract.selections(crossdesign=cross.design, relmatrix = rel.mat.cross,
                                              past.tgv = progeny9.TGV, past.phenos = progeny9.phenos, parentinfo = progeny9.extractions,
                                              progenyinfo = progeny10, progenyTGV = progeny10.TGV,progenyphenos = progeny10.phenos,
                                              selection.strategy = sel.strategy, 
                                              among.family.selection = af.selection,within.family.selection = wf.selection,numSelections.among.family = num.sel.af)
  
  
  genetic.gain.mine[1,10] <- mean(progeny10.extractions$select.genval)-founder.TGV
  genotypic.variance[1,10]<- var(progeny10.extractions$select.genval)
  phenotypic.gain.mine[1,10] <-mean(progeny10.extractions$selection.phenos)-founder.phenos
  mean.select.inbreeding[1,10]<- mean(progeny10.extractions$select.inbred.level)
  mean.pop.inbreeding[1,10]<- mean(progeny10.extractions$prog.inbred.level)
  bulmer.effect[1,10]<- progeny10.extractions$bulmer.effect
  coancest.threshold[1,10] <- cross.design$coancestry.threshold
  delt.alleles[1,10] <- sum(progeny10.extractions$delt.alleles)/NumParents
  max.delt.alleles[1,10] <- max(progeny10.extractions$delt.allele)
  min.delt.alleles[1,10] <- min(progeny10.extractions$delt.allele)
  
  rm(progeny9.TGV,progeny9.extractions, progeny9.phenos,progeny10); gc()
  
  
#### List for output objects####
  newList <- list("g"=genetic.gain.mine,"gv"=genotypic.variance, "pg"=phenotypic.gain.mine,"si"=mean.select.inbreeding,
                  "pi"=mean.pop.inbreeding,"be"=bulmer.effect, "ct"=coancest.threshold, "maxd"=max.delt.alleles,"mind"=min.delt.alleles,
                  "da"=delt.alleles)
  the.name <- paste("AGA.64.OP.pop", pop.num, sep="")
  assign(x=the.name,value=newList)
  save(list=the.name,file=paste("/media/titan/ARF/results/AAA.64/output/",the.name,".RDA",sep=""))
}
