# Run this on Domino ####
source('/mnt/simulator/create_map.R')
source('/mnt/simulator/create_parents.R')
source('/mnt/simulator/calc_TGV.R')
source('/mnt/simulator/sim_phenos.R')
source('/mnt/simulator/create_cross_design.R')
source('/mnt/simulator/make_crosses.R')
source('/mnt/simulator/extract_selections.R')
set.seed(45342)
#set.seed(345342)
#seed.list <- c(sample(2342:2352345,20,replace=F))

#g_maps <- lapply(1:20,function(each){
#  set.seed(seed.list[each])
genetic.map <- create_genetic_map(num.chromos = 12,map.length = 1800,num.markers = (10*12),total.qtl = (500*12), save=F,map.dist = "haldane",
                                  num.snpqtl = 500,distribute.loci = "even",marker.distribution = "equally-spaced",snp.qtl.maf = c(0.01,0.02))
parents <- create_parents(map.info = genetic.map,num.parents = 2,max.delt.allele = 0,inbred.parents = T,par.markers.unique = T)
parents.TGV <- calc_TGV(geno.info = parents,map.info = genetic.map,A = 1,a = -100,dom.coeff = 0,founder = T,par.markers.unique = T)
parents.phenos <- sim_phenos(TGV.object = parents.TGV,h2 = .3)

first.cross <- create_cross_design(parent.info = parents,prog.percross = 200,generation = 1,mating.design = "RandomMating")
progeny <- make_crosses(parent.info = parents,map.info = genetic.map,cross.design = first.cross,num.cores = 2)
progeny.TGV <- calc_TGV(geno.info = progeny,map.info = genetic.map,cross.design = first.cross,A = 1,a = -100,dom.coeff = 0,par.markers.unique = T)
progeny.Phenos <- sim_phenos(h2 = .3,TGV.object = progeny.TGV)
Selections <- extract_selections(map.info = genetic.map,cross.design = first.cross,relationship.matrix.type="pedigree",
                                 past.tgv = parents.TGV,past.phenos = parents.phenos,
                                 parent.info = parents,progeny.info = progeny,progeny.TGV = progeny.TGV,
                                 progeny.phenos = progeny.Phenos,among.family.selection = "Phenotype",
                                 num.selections.within.family = 1)


second.cross <- create_cross_design(prog.percross = 300,generation = 2,mating.design = "Self",parent.info = Selections)
progeny <- make_crosses(parent.info = Selections,map.info = genetic.map,cross.design = second.cross,num.cores = 2)
progeny.TGV <- calc_TGV(geno.info = progeny,map.info = genetic.map,cross.design = second.cross,A = 1,a = -100,dom.coeff = 0,par.markers.unique = T)

test = data.frame(t(progeny.TGV$markers.matrix),stringsAsFactors = F)
for(i in 1:ncol(test)){
  test[which(test[,i] == 0),i] <- "A"
  test[which(test[,i] == 2),i] <- "B"
  test[which(test[,i] == 1),i] <- "X"
}
#install.packages("ASMap")
library(ASMap)
ut <- mstmap.data.frame(object = test,pop.type = "RIL2",as.cross = T,trace = F)

#})
geno <- ut$geno

these.chr <- unique(genetic.map$chr)
chr.lengths <- c()
for(each in 1:12){chr.lengths <- c(chr.lengths,genetic.map$pos[tail(which(genetic.map$chr == these.chr[each]),1)])}
sum(chr.lengths)
linkage_map_chr_lengths <- unlist(lapply(geno,function(x)x$map[length(x$map)]))
sum(linkage_map_chr_lengths)
mean(chr.lengths - linkage_map_chr_lengths)

?heatMap
m <- heatMap(x=ut,what="both")

map2 <- quickEst(ut, map.function = "haldane")
plot.map(map2)
map2 <- subset(map2, chr = names(nmar(map2))[6:15])
R> plot.map(map1, map2)

these.chr <- unique(genetic.map$chr)
m_dist <- c()
for(each in 1:12){
m_dist <- c(m_dist,dist(genetic.map$pos[which(genetic.map$chr == these.chr[each] & genetic.map$types=="m")])[c(1,10,18,25,31,36,40,43,45)])
}
mean(m_dist)
sd(m_dist)
gm_dist <- unlist(lapply(geno,function(x)dist(x$map)[c(1,10,18,25,31,36,40,43,45)]))
sd(gm_dist)
mean(gm_dist)
