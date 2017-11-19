for(pop.num in 301:1000){
set.seed(Sys.time())
#Create base population for comparisons
#Load package####
source('/media/titan/ARF/Simulator.Package/breeding.simulator.R')
#Create Genetic Map####
map <- create.map(num.markers=120,numSNPQTL = 1960, numchromos = 12, maplength = 1800, totalQTL=2500,
                  totalloci = signif((120+2500)/12,digits = 2)*12,chromosomesizerange = .2)
parents <- create.parents(mapinfo=map, numparents = 280,max.delt.allele = 14,heterozygous.markers = T)
parents.TGV <- create.parents.TGV(parents=parents,mapinfo=map,A=1,a=-100,dom.coeff = 1)
rm(.Random.seed);parents.phenos <- create.parents.phenos(parents.TGV = parents.TGV, h2 = .3)
Env.var <- sqrt(var(parents.TGV$genetic.values)/.3)
founder.gvar <- var(parents.TGV$genetic.values)
founder.TGV <- mean(parents.TGV$genetic.values)
founder.phenos <- mean(parents.phenos$phenos)


OP.parent.phenos <- select.OP.Phenotype.Values(mapinfo=map, parentinfo = parents, parents.TGV = parents.TGV,
                                               parent.phenos = parents.phenos,num.select = 280,crossprog=1, dom.coeff = 1,
                                               A=1,a=-100, h2=.3 , run.parallel = T, n.cores = 20)
pop.name <- paste("/media/titan/ARF/base.population/base.population.",pop.num,".RData",sep="")
save.image(file = pop.name)}
