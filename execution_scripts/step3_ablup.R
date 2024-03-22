# Run this on Domino ####
source('/mnt/results/functions/calc_TGV.R')
source('/mnt/results/functions/sim_phenos.R')
source('/mnt/results/functions/create_cross_design.R')
source('/mnt/results/functions/make_crosses.R')
source('/mnt/results/functions/calc_TGV.R')
source('/mnt/results/functions/extract_selections_ibd.R')

## Run Values ####
NumParents=64
prog.per.cross = 60
af.selection = "ABLUP" #ABLUP GBLUP_FI snpibd
rel.mat.cross = "pedigree"
num.sel.af = 64
num.of.cores = 5
mt.design="AssortiveMating"

#Values for alleles and trait variation
Major.value =1
Minor.value =-100
Dominance.Coeff=1
indiv.tree.h2 = .3
load("/mnt/results/3k/all.seeds.RData")
for(pop.num in 1:10){
 # res <- mclapply(111:118,function(pop.num){
#pt <- proc.time()
#pop.num=1
    this.gen <- seed.vector[pop.num]
    
    load.pop <- paste0("/mnt/results/3k/base_populations_op_test/base_population_op_test_",pop.num,".RData")
    load(load.pop)
    
    ###Create matrices to hold outputs####
    simulations=1; generations=5
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
    op.test$mean.parent.phenos <- op.test$mean.parent.phenos[1:NumParents]
    op.test$mean.parent.tgv <- op.test$mean.parent.tgv[1:NumParents]
    op.test$delt.alleles <- op.test$delt.alleles[1:NumParents]
    op.test$genetic.values <- op.test$genetic.values[1:NumParents]
    op.test$phenos <- op.test$phenos[1:NumParents]
    op.test$genos.3d <- op.test$genos.3d[,1:NumParents,]
    op.test$marker.matrix <- op.test$marker.matrix[1:NumParents,]
    
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
    rm(progeny1); gc()
    
    ####Generation 2#####
    set.seed(this.gen+10);
    cross.design <- create_cross_design(parent.info = progeny1.extractions,prog.percross = prog.per.cross,generation  = 2,coancest.thresh=T)
    set.seed(this.gen+12);progeny2 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny1.extractions,num.cores = num.of.cores)
    progeny2.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny2,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
    set.seed(this.gen+13); progeny2.phenos <- sim_phenos(TGV.object = progeny2.TGV,h2 = indiv.tree.h2)
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
    
    rm(progeny1.TGV,progeny1.extractions, progeny1.phenos,progeny2);gc()
    
    ####Generations 3####
    set.seed(this.gen+21); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny2.extractions,prog.percross = prog.per.cross,generation  = 3)
    set.seed(this.gen+22); progeny3 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny2.extractions,num.cores = num.of.cores)
    progeny3.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny3,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
    set.seed(this.gen+23); progeny3.phenos <- sim_phenos(TGV.object = progeny3.TGV,h2 = indiv.tree.h2)
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
    
    rm(progeny2.TGV,progeny2.extractions, progeny2.phenos,progeny3);gc()
    
    ####Generations 4####
    set.seed(this.gen+31); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny3.extractions,prog.percross = prog.per.cross,generation  = 4)
    set.seed(this.gen+32); progeny4 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny3.extractions,num.cores = num.of.cores)
    progeny4.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny4,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
    set.seed(this.gen+33); progeny4.phenos <- sim_phenos(TGV.object = progeny4.TGV,h2 = indiv.tree.h2)
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
    
    rm(progeny3.TGV,progeny3.extractions, progeny3.phenos,progeny4);gc()
    
    
    ####Generations 5####
    set.seed(this.gen+41); cross.design <- create_cross_design(coancest.thresh=T,parent.info = progeny4.extractions,prog.percross = prog.per.cross,generation  = 5)
    set.seed(this.gen+42); progeny5 <- make_crosses(map.info=base_pop_data$genetic.map, cross.design = cross.design, parent.info = progeny4.extractions,num.cores = num.of.cores)
    progeny5.TGV <- calc_TGV(map.info = base_pop_data$genetic.map,geno.info = progeny5,cross.design = cross.design,A = Major.value,a = Minor.value,dom.coeff = Dominance.Coeff)
    set.seed(this.gen+43); progeny5.phenos <- sim_phenos(TGV.object = progeny5.TGV,h2 = indiv.tree.h2)
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
    
    rm(progeny4.TGV,progeny4.extractions, progeny4.phenos,progeny5);gc()
    
    #### List for output objects####
    newList <- list("g"=genetic.gain.mine,"gv"=genotypic.variance, "pg"=phenotypic.gain.mine,"si"=mean.select.inbreeding,
                    "pi"=mean.pop.inbreeding,"be"=bulmer.effect, "ct"=coancest.threshold, "maxd"=max.delt.alleles,"mind"=min.delt.alleles,
                    "da"=delt.alleles)
    
    the.name <- paste("ABLUP.nodom.64.OP.pop_", pop.num, sep="")
    assign(x=the.name,value=newList)
    save(list=the.name,file=paste("/mnt/results/3k/64/ABLUP/",the.name,".RData",sep=""))
    print(pop.num)
}   
#    proc.time() - pt
#  },mc.cores=8)
  
  
#  the.name <- paste("ABLUP.nodom.64.OP.pop", pop.num, sep="")
#  assign(x=the.name,value=newList)
#  save(list=the.name,file=paste("./data/results3/64/ABLUP/",the.name,".RData",sep=""))
#}
