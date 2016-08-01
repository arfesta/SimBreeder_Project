###Set Working Directory, Save loctation, & Source file####
setwd("~/Desktop/SIMULATOR") #MAC
Rep="~/Desktop/SIMULATOR/data/miniTest/" #SBMP or SELF
source("~/Desktop/SIMULATOR/FINAL.create.mapV9.R")
Save=FALSE

###Create matrices to hold outputs####
simulations=1; generations=11
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

###User Defined Inputs for Map####
Haldane=TRUE        # Use either Haldane Function or 
map.length=1800     # Total map length in centimorgans
num.chromos=12      # Number of chromosomes
total.loci= 60000   # Total number of loci to be distributed to equally distributed to chromosomes (must be divisble by num.chromos)
rQTL=1000           # Total number of random QTL involved in a trait 
dom.coeff=1         # Dominance level: Range (0-1)
NumParents=300      # Number of parents which are created to be inital founder population
MaxDeltAlleles=14   # The max number of deleterious alleles per parent
NumSNPQTL=(NumParents*2)*(MaxDeltAlleles/2)      # Total number of QTL which are affected by dominance (must be less than TotalQTL)
TotalQTL= rQTL+NumSNPQTL
NumMarkers=2500     # Number of markers that are used in breeding strategies 
Heritability=c(.3)  # Heritablity of the trait under consideration
prog.percross=20    # Number of progeny to be made for each cross
CrossFile1="~/Desktop/SIMULATOR/prototype.txt" #Input crossfile for 1st generation
#SNPeffects= rep(0, 100) # This can be set to specify the effects of SNPS (length SNP effects must match NumSNPQTL)
#useSNPeffects=T         # Operator switch for using the SNP effects

####CrossDesign Options for 1st Generation####
crossfileinput=T; random.mating=F; self=F; sp.mateselect=F; mp.mateselect=F

####Selection Options for all generations####
P=F; Ablup=T; Gblup=F; hblup=F        # Method for calculating BV
mateselect=T; sbpm=F; selfp=F; sbpw=F # Method for making selections
#Should Selfing test be conducted each generation to identify deleterious alleles?
SELFTEST=F 

###Create Founder Population####
map<- create.Map()
parents <- create.Parents()
parents.TGV <- create.TGV(prefix="parents.TGV", gen=0)
parents.phenos <- create.Unscaled.Phenos(TGV=parents.TGV, gen=0)

Env.var <- parents.phenos$E.var # store multipier for creating on future generation phenos
a <- sqrt(50)/sqrt(var(parents.TGV$genetic.values)) #scalar factor to set genetic varaince to 50
b <- 50/(var(parents.TGV$genetic.values))  #scalar factor to set mean genetic values to 50
founder.TGV <- mean(parents.TGV$genetic.values) #mean of raw founder genetic values
founder.phenos <- mean(parents.phenos$phenos) # mean of raw founder phenotypes

genetic.gain.mine[1,1] <- 0
genotypic.variance[1,1] <- var(parents.TGV$genetic.values) 
phenotypic.gain.mine[1,1] <- 0
mean.select.inbreeding[1,1] <- 0 
mean.pop.inbreeding[1,1]<- 0
bulmer.effect[1,1] <- 0
coancest.threshold[1,1] <- 0
delt.alleles[1,1] <- sum(parents$delt.allele)/300
max.delt.alleles[1,1] <- max(parents$delt.allele)
min.delt.alleles[1,1] <- min(parents$delt.allele)

#### Create 1st gen progeny#######
cross.design <- createCrossDesign()
if (SELFTEST){selfed.progeny <- mclapply(rep(1:NumParents,1),createSelfProgeny,mc.cores=18); names(selfed.progeny) <- paste0(1:NumParents)}
#no_cores <- 3; CL <- makeCluster(4); clusterExport(CL,list("map","prog.percross","parents"),envir=environment())
progeny1 <- apply(cross.design$crossdesign,1,makecrosses); #stopCluster(CL)
progeny1 <- create.genos()
progeny1.TGV <- create.TGV(prefix="progeny1.TGV",data1=progeny1, gen=1)
rm(.Random.seed); progeny1.phenos <- create.Unscaled.Phenos(TGV=progeny1.TGV, gen=1)
progeny1.extractions <- Extract.Selections(parents.tgv=parents.TGV, gen=1,
                                         progenyinfo=progeny1,progenyTGV=progeny1.TGV, prefix="selection1.pop")

genetic.gain.mine[1,2] <- mean(progeny1.extractions$select.genval)-founder.TGV
genotypic.variance[1,2] <- var(progeny1.TGV$genetic.values) 
phenotypic.gain.mine[1,2] <-mean(progeny1.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,2] <- mean(progeny1.extractions$select.inbred.level)  
mean.pop.inbreeding[1,2]<- mean(progeny1.extractions$prog.inbred.level)
bulmer.effect[1,2] <- progeny1.extractions$bulmer.effect
coancest.threshold[1,2] <- 0
delt.alleles[1,2] <- sum(progeny1.extractions$delt.alleles)/300
max.delt.alleles[1,2] <- max(progeny1.extractions$delt.allele)
min.delt.alleles[1,2] <- min(progeny1.extractions$delt.allele)
rm(progeny1,progeny1.phenos)

####CrossDesign Options for Future Generations####
crossfileinput=F; random.mating=F; self=F; sp.mateselect=T; mp.mateselect=F

####Generation 2#####
cross.design <- createCrossDesign(parentinfo=progeny1.extractions, gen=2)
if (SELFTEST==T) {selfed.progeny <- mclapply(rep(1:NumParents,1),createSelfProgeny,MoreArgs=c(parentinfo=progeny1.extractions,phenos=progeny1.phenos$phenos, TGV=progeny1.TGV) ,mc.cores=3);names(selfed.progeny) <- paste(progeny1.extractions$select.ped.ids)}
#CL <- makeCluster(no_cores);clusterExport(CL,list("map","prog.percross","progeny1.extractions"),envir=environment())
progeny2 <- apply(cross.design$crossdesign,1,makecrosses, gen=2,parentinfo=progeny1.extractions); #stopCluster(CL)
progeny2 <- create.genos(progeny=progeny2)
progeny2.TGV <- create.TGV(prefix="progeny2.TGV",data1=progeny2, gen=2)
rm(.Random.seed); progeny2.phenos <- create.Unscaled.Phenos(TGV=progeny2.TGV, gen=2)
progeny2.extractions <- Extract.Selections(parents.tgv=progeny1.TGV, gen=2, past.phenos=progeny1.extractions,
                                           progenyinfo=progeny2,progenyTGV=progeny2.TGV, progenyphenos= progeny2.phenos, prefix="selection2.pop")

genetic.gain.mine[1,3] <- mean(progeny2.extractions$select.genval)-founder.TGV
genotypic.variance[1,3]<- var(progeny2.TGV$genetic.values) 
phenotypic.gain.mine[1,3] <-mean(progeny2.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,3]<- mean(progeny2.extractions$select.inbred.level)  
mean.pop.inbreeding[1,3]<- mean(progeny2.extractions$prog.inbred.level)
bulmer.effect[1,3]<- progeny2.extractions$bulmer.effect
coancest.threshold[1,3] <- cross.design$coancestry.threshold
delt.alleles[1,3] <- sum(progeny2.extractions$delt.alleles)/300
max.delt.alleles[1,3] <- max(progeny2.extractions$delt.allele)
min.delt.alleles[1,3] <- min(progeny2.extractions$delt.allele)

rm(progeny1.TGV,progeny1.extractions, progeny2, progeny2.phenos)

####Generations 3####
cross.design <- createCrossDesign(parentinfo=progeny2.extractions, gen=3)
if (SELFTEST==T) {selfed.progeny <- mclapply(rep(1:NumParents,1),createSelfProgeny,MoreArgs=c(parentinfo=progeny2.extractions,phenos=progeny2.phenos$phenos, TGV=progeny2.TGV),mc.cores=3) ;names(selfed.progeny) <- paste(progeny2.extractions$select.ped.ids)}
#CL <- makeCluster(no_cores); clusterExport(CL,list("map","prog.percross","progeny2.extractions"),envir=environment())
progeny3 <- apply(cross.design$crossdesign,1,makecrosses, gen=3,parentinfo=progeny2.extractions); #stopCluster(CL)
progeny3 <- create.genos(progeny=progeny3)
progeny3.TGV <- create.TGV(prefix="progeny3.TGV",data1=progeny3, gen=3)
rm(.Random.seed); progeny3.phenos <- create.Unscaled.Phenos(TGV=progeny3.TGV, gen=3)
progeny3.extractions <- Extract.Selections(parents.tgv=progeny2.TGV, gen=3, past.phenos=progeny2.extractions,
                                           progenyinfo=progeny3,progenyTGV=progeny3.TGV, progenyphenos= progeny3.phenos, prefix="selection3.pop")


genetic.gain.mine[1,4] <- mean(progeny3.extractions$select.genval)-founder.TGV
genotypic.variance[1,4]<- var(progeny3.TGV$genetic.values) 
phenotypic.gain.mine[1,4] <-mean(progeny3.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,4]<- mean(progeny3.extractions$select.inbred.level)  
mean.pop.inbreeding[1,4]<- mean(progeny3.extractions$prog.inbred.level)
bulmer.effect[1,4]<- progeny3.extractions$bulmer.effect
coancest.threshold[1,4] <- cross.design$coancestry.threshold
delt.alleles[1,4] <- sum(progeny3.extractions$delt.alleles)/300
max.delt.alleles[1,4] <- max(progeny3.extractions$delt.allele)
min.delt.alleles[1,4] <- min(progeny3.extractions$delt.allele)

rm(progeny2.extractions,progeny2.TGV,progeny3,progeny3.phenos)

####Generation 4####
cross.design <- createCrossDesign(parentinfo=progeny3.extractions, gen=4)
if (SELFTEST==T) {selfed.progeny <- mclapply(rep(1:NumParents,1),createSelfProgeny,MoreArgs=c(parentinfo=progeny3.extractions,phenos=progeny3.phenos$phenos, TGV=progeny3.TGV),mc.cores=3);names(selfed.progeny) <- paste(progeny3.extractions$select.ped.ids)}
#CL <- makeCluster(no_cores); clusterExport(CL,list("map","prog.percross","progeny3.extractions"),envir=environment())
progeny4 <- apply(cross.design$crossdesign,1,makecrosses, gen=4,parentinfo=progeny3.extractions);# stopCluster(CL)
progeny4 <- create.genos(progeny=progeny4)
progeny4.TGV <- create.TGV(prefix="progeny4.TGV",data1=progeny4, gen=4)
rm(.Random.seed); progeny4.phenos <- create.Unscaled.Phenos(TGV=progeny4.TGV, gen=4)
progeny4.extractions <- Extract.Selections(parents.tgv=progeny3.TGV, gen=4, past.phenos=progeny3.extractions,
                                           progenyinfo=progeny4,progenyTGV=progeny4.TGV, progenyphenos= progeny4.phenos, prefix="selection2.pop")


genetic.gain.mine[1,5] <- mean(progeny4.extractions$select.genval)-founder.TGV
genotypic.variance[1,5]<- var(progeny4.TGV$genetic.values) 
phenotypic.gain.mine[1,5] <-mean(progeny4.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,5]<- mean(progeny4.extractions$select.inbred.level)  
mean.pop.inbreeding[1,5]<- mean(progeny4.extractions$prog.inbred.level)
bulmer.effect[1,5]<- progeny4.extractions$bulmer.effect
coancest.threshold[1,5] <- cross.design$coancestry.threshold
delt.alleles[1,5] <- sum(progeny4.extractions$delt.alleles)/300
max.delt.alleles[1,5] <- max(progeny4.extractions$delt.allele)
min.delt.alleles[1,5] <- min(progeny4.extractions$delt.allele)
rm(progeny3.extractions,progeny3.TGV,progeny4,progeny4.phenos)

####Generation 5####
if (SELFTEST==T) {system.time(selfed.progeny <- sapply(1:NumParents,createSelfProgeny,parentinfo=progeny4.extractions,phenos=progeny4.phenos$phenos, TGV=progeny4.TGV));names(selfed.progeny) <- paste(progeny4.extractions$select.ped.ids)}
cross.design <- createCrossDesign(parentinfo=progeny4.extractions, gen=5)
#CL <- makeCluster(no_cores); clusterExport(CL,list("map","prog.percross","progeny4.extractions"),envir=environment())
progeny5 <- apply(cross.design$crossdesign,1,makecrosses, gen=5,parentinfo=progeny4.extractions); #stopCluster(CL)
progeny5 <- create.genos(progeny=progeny5)
progeny5.TGV <- create.TGV(prefix="progeny5.TGV",data1=progeny5, gen=5)
rm(.Random.seed); progeny5.phenos <- create.Unscaled.Phenos(TGV=progeny5.TGV, gen=5)
progeny5.extractions <- Extract.Selections(parents.tgv=progeny4.TGV, gen=5, past.phenos=progeny4.extractions,
                                           progenyinfo=progeny5,progenyTGV=progeny5.TGV, progenyphenos= progeny5.phenos, prefix="selection2.pop")


genetic.gain.mine[1,6] <- mean(progeny5.extractions$select.genval)-founder.TGV
genotypic.variance[1,6]<- var(progeny5.TGV$genetic.values) 
phenotypic.gain.mine[1,6] <-mean(progeny5.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,6]<- mean(progeny5.extractions$select.inbred.level)  
mean.pop.inbreeding[1,6]<- mean(progeny5.extractions$prog.inbred.level)
bulmer.effect[1,6]<- progeny5.extractions$bulmer.effect
coancest.threshold[1,6] <- cross.design$coancestry.threshold
delt.alleles[1,6] <- sum(progeny5.extractions$delt.alleles)/300
max.delt.alleles[1,6] <- max(progeny5.extractions$delt.allele)
min.delt.alleles[1,6] <- min(progeny5.extractions$delt.allele)
rm(progeny4.extractions,progeny4.TGV,progeny5,progeny5.phenos)

####Generation 6####
if (SELFTEST==T) {selfed.progeny <- mclapply(rep(1:NumParents,1),createSelfProgeny,MoreArgs=c(parentinfo=progeny5.extractions,phenos=progeny5.phenos$phenos, TGV=progeny5.TGV),mc.cores=3);names(selfed.progeny) <- paste(progeny5.extractions$select.ped.ids)}
cross.design <- createCrossDesign(parentinfo=progeny5.extractions, gen=6)
#CL <- makeCluster(no_cores); clusterExport(CL,list("map","prog.percross","progeny5.extractions"),envir=environment())
progeny6 <- apply(cross.design$crossdesign,1,makecrosses, gen=6,parentinfo=progeny5.extractions); #stopCluster(CL)
progeny6 <- create.genos(progeny=progeny6)
progeny6.TGV <- create.TGV(prefix="progeny6.TGV",data1=progeny6, gen=6)
rm(.Random.seed); progeny6.phenos <- create.Unscaled.Phenos(TGV=progeny6.TGV, gen=6)
progeny6.extractions <- Extract.Selections(parents.tgv=progeny5.TGV, gen=6, past.phenos=progeny5.extractions,
                                           progenyinfo=progeny6,progenyTGV=progeny6.TGV, progenyphenos= progeny6.phenos, prefix="selection6.pop")


genetic.gain.mine[1,7] <- mean(progeny6.extractions$select.genval)-founder.TGV
genotypic.variance[1,7]<- var(progeny6.TGV$genetic.values) 
phenotypic.gain.mine[1,7] <-mean(progeny6.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,7]<- mean(progeny6.extractions$select.inbred.level)  
mean.pop.inbreeding[1,7]<- mean(progeny6.extractions$prog.inbred.level)
bulmer.effect[1,7]<- progeny6.extractions$bulmer.effect
coancest.threshold[1,7] <- cross.design$coancestry.threshold
delt.alleles[1,7] <- sum(progeny6.extractions$delt.alleles)/300
max.delt.alleles[1,7] <- max(progeny6.extractions$delt.allele)
min.delt.alleles[1,7] <- min(progeny6.extractions$delt.allele)

rm(progeny5.extractions,progeny5.TGV,progeny6,progeny6.phenos)

####Generation 7####
if (SELFTEST==T) {system.time(selfed.progeny <- sapply(1:NumParents,createSelfProgeny,parentinfo=progeny6.extractions,phenos=progeny6.phenos$phenos, TGV=progeny6.TGV));names(selfed.progeny) <- paste(progeny6.extractions$select.ped.ids)}
cross.design <- createCrossDesign(parentinfo=progeny6.extractions, gen=7)
#CL <- makeCluster(no_cores); clusterExport(CL,list("map","prog.percross","progeny6.extractions"),envir=environment())
progeny7 <- apply(cross.design$crossdesign,1,makecrosses, gen=7,parentinfo=progeny6.extractions); #stopCluster(CL)
progeny7 <- create.genos(progeny=progeny7)
progeny7.TGV <- create.TGV(prefix="progeny7.TGV",data1=progeny7, gen=7)
rm(.Random.seed); progeny7.phenos <- create.Unscaled.Phenos(TGV=progeny7.TGV, gen=7)
progeny7.extractions <- Extract.Selections(parents.tgv=progeny6.TGV, gen=7, past.phenos=progeny6.extractions,
                                           progenyinfo=progeny7,progenyTGV=progeny7.TGV, progenyphenos= progeny7.phenos, prefix="selection7.pop")


genetic.gain.mine[1,8] <- mean(progeny7.extractions$select.genval)-founder.TGV
genotypic.variance[1,8]<- var(progeny7.TGV$genetic.values) 
phenotypic.gain.mine[1,8] <-mean(progeny7.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,8]<- mean(progeny7.extractions$select.inbred.level)  
mean.pop.inbreeding[1,8]<- mean(progeny7.extractions$prog.inbred.level)
bulmer.effect[1,8]<- progeny7.extractions$bulmer.effect
coancest.threshold[1,8] <- cross.design$coancestry.threshold
delt.alleles[1,8] <- sum(progeny7.extractions$delt.alleles)/300
max.delt.alleles[1,8] <- max(progeny7.extractions$delt.allele)
min.delt.alleles[1,8] <- min(progeny7.extractions$delt.allele)

rm(progeny6.extractions,progeny6.TGV,progeny7,progeny7.phenos)

####Generation 8####
if (SELFTEST==T) {system.time(selfed.progeny <- sapply(1:NumParents,createSelfProgeny,parentinfo=progeny7.extractions,phenos=progeny7.phenos$phenos, TGV=progeny7.TGV));names(selfed.progeny) <- paste(progeny7.extractions$select.ped.ids)}
cross.design <- createCrossDesign(parentinfo=progeny7.extractions, gen=8)
#CL <- makeCluster(no_cores); clusterExport(CL,list("map","prog.percross","progeny7.extractions"),envir=environment())
progeny8 <- apply(cross.design$crossdesign,1,makecrosses, gen=8,parentinfo=progeny7.extractions); #stopCluster(CL)
progeny8 <- create.genos(progeny=progeny8)
progeny8.TGV <- create.TGV(prefix="progeny8.TGV",data1=progeny8, gen=8)
rm(.Random.seed); progeny8.phenos <- create.Unscaled.Phenos(TGV=progeny8.TGV, gen=8)
progeny8.extractions <- Extract.Selections(parents.tgv=progeny7.TGV, gen=8, past.phenos=progeny7.extractions,
                                           progenyinfo=progeny8,progenyTGV=progeny8.TGV, progenyphenos= progeny8.phenos, prefix="selection2.pop")


genetic.gain.mine[1,9] <- mean(progeny8.extractions$select.genval)-founder.TGV
genotypic.variance[1,9]<- var(progeny8.TGV$genetic.values) 
phenotypic.gain.mine[1,9] <-mean(progeny8.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,9]<- mean(progeny8.extractions$select.inbred.level)  
mean.pop.inbreeding[1,9]<- mean(progeny8.extractions$prog.inbred.level)
bulmer.effect[1,9]<- progeny8.extractions$bulmer.effect
coancest.threshold[1,9] <- cross.design$coancestry.threshold
delt.alleles[1,9] <- sum(progeny8.extractions$delt.alleles)/300
max.delt.alleles[1,9] <- max(progeny8.extractions$delt.allele)
min.delt.alleles[1,9] <- min(progeny8.extractions$delt.allele)

rm(progeny7.extractions,progeny7.TGV,progeny8,progeny8.phenos)

####Generation 9####
if (SELFTEST==T) {system.time(selfed.progeny <- sapply(1:NumParents,createSelfProgeny,parentinfo=progeny8.extractions,phenos=progeny8.phenos$phenos, TGV=progeny8.TGV));names(selfed.progeny) <- paste(progeny8.extractions$select.ped.ids)}
cross.design <- createCrossDesign(parentinfo=progeny8.extractions, gen=9)
#CL <- makeCluster(no_cores); clusterExport(CL,list("map","prog.percross","progeny8.extractions"),envir=environment())
progeny9 <- apply(cross.design$crossdesign,1,makecrosses, gen=9,parentinfo=progeny8.extractions); #stopCluster(CL)
progeny9 <- create.genos(progeny=progeny9)
progeny9.TGV <- create.TGV(prefix="progeny9.TGV",data1=progeny9, gen=9)
rm(.Random.seed); progeny9.phenos <- create.Unscaled.Phenos(TGV=progeny9.TGV, gen=9)
progeny9.extractions <- Extract.Selections(parents.tgv=progeny8.TGV, gen=9, past.phenos=progeny8.extractions,
                                           progenyinfo=progeny9,progenyTGV=progeny9.TGV, progenyphenos= progeny9.phenos, prefix="selection9.pop")


genetic.gain.mine[1,10] <- mean(progeny9.extractions$select.genval)-founder.TGV
genotypic.variance[1,10]<- var(progeny9.TGV$genetic.values) 
phenotypic.gain.mine[1,10] <-mean(progeny9.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,10]<- mean(progeny9.extractions$select.inbred.level)  
mean.pop.inbreeding[1,10]<- mean(progeny9.extractions$prog.inbred.level)
bulmer.effect[1,10]<- progeny9.extractions$bulmer.effect
coancest.threshold[1,10] <- cross.design$coancestry.threshold
delt.alleles[1,10] <- sum(progeny9.extractions$delt.alleles)/300
max.delt.alleles[1,10] <- max(progeny9.extractions$delt.allele)
min.delt.alleles[1,10] <- min(progeny9.extractions$delt.allele)
rm(progeny8.extractions,progeny8.TGV,progeny9,progeny9.phenos)

####Generation 10####
if (SELFTEST==T) {system.time(selfed.progeny <- sapply(1:NumParents,createSelfProgeny,parentinfo=progeny9.extractions,phenos=progeny9.phenos$phenos, TGV=progeny9.TGV));names(selfed.progeny) <- paste(progeny9.extractions$select.ped.ids)}
cross.design <- createCrossDesign(parentinfo=progeny9.extractions, gen=10)
#CL <- makeCluster(no_cores); clusterExport(CL,list("map","prog.percross","progeny9.extractions"),envir=environment())
progeny10 <- apply(cross.design$crossdesign,1,makecrosses, gen=10,parentinfo=progeny9.extractions); #stopCluster(CL)
progeny10 <- create.genos(progeny=progeny10)
progeny10.TGV <- create.TGV(prefix="progeny10.TGV",data1=progeny10, gen=10)
rm(.Random.seed); progeny10.phenos <- create.Unscaled.Phenos(TGV=progeny10.TGV, gen=10)
progeny10.extractions <- Extract.Selections(parents.tgv=progeny9.TGV, gen=10, past.phenos=progeny9.extractions,
                                           progenyinfo=progeny10,progenyTGV=progeny10.TGV, progenyphenos= progeny10.phenos, prefix="selection10.pop")


genetic.gain.mine[1,11] <- mean(progeny10.extractions$select.genval)-founder.TGV
genotypic.variance[1,11]<- var(progeny10.TGV$genetic.values) 
phenotypic.gain.mine[1,11] <-mean(progeny10.extractions$selection.phenos)-founder.phenos
mean.select.inbreeding[1,11]<- mean(progeny10.extractions$select.inbred.level)  
mean.pop.inbreeding[1,11]<- mean(progeny10.extractions$prog.inbred.level)
bulmer.effect[1,11]<- progeny10.extractions$bulmer.effect
coancest.threshold[1,11] <- cross.design$coancestry.threshold
delt.alleles[1,11] <- sum(progeny10.extractions$delt.alleles)/300
max.delt.alleles[1,11] <- max(progeny10.extractions$delt.allele)
min.delt.alleles[1,11] <- min(progeny10.extractions$delt.allele)

rm(progeny9.extractions,progeny9.TGV,progeny10,progeny10.phenos)

#### List for output objects####
newList <- list("g"=genetic.gain.mine,"gv"=genotypic.variance, "pg"=phenotypic.gain.mine,"si"=mean.select.inbreeding,
                "pi"=mean.pop.inbreeding,"be"=bulmer.effect, "ct"=coancest.threshold, "maxd"=max.delt.alleles,"mind"=min.delt.alleles,
                "da"=delt.alleles)
newList
#save(newList,file="~/Desktop/SP.MateSelect.GBLUP.2.RDA")
