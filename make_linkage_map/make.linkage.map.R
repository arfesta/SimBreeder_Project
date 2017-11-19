library(SimBreeder)

yo <- create_genetic_map(num.chromos = 12,map.length = 1800, num.markers = 120 ,total.QTL = 120,marker.maf = .5, 
                         marker.distribution = "equally-spaced",
                         num.SNPQTL = 120, map.dist = "kosambi",chromosome.size.range = .2)

g.map <-yo$genetic.map
the.parent <- create_parents(map.info = yo,num.parents = 1,max.delt.allele = 8,inbred.parents = T)
parent.genos = data.frame(the.parent$genos.3d)
cross <- create_cross_design(parent.info = the.parent,prog.percross = 200,generation = 1,mating.design = "Self" )
cc = cross$cross.design

cross.object <- make_crosses(parent.info = the.parent,map.info = yo,cross.design = cross,run.parallel = T,num.cores = 2)

check = cross.object$genos.3d[,,1]

out <- calc_progeny_TGV(prog.info = cross.object,map.info = yo,cross.design = cross,A = 1,a = -100,dom.coeff = 1)

test = data.frame(t(out$markers.matrix),stringsAsFactors = F)
for(i in 1:ncol(test)){
  test[which(test[,i] == 0),i] <- "A"
  test[which(test[,i] == 2),i] <- "B"
  test[which(test[,i] == 1),i] <- "X"
}

library(ASMap)
ut <- mstmap.data.frame(object = test,pop.type = "RIL2",as.cross = T,trace = F)

geno <- ut$geno

heatMap(x=ut)
#statGen(cross = ut)
plot.geno(x = ut,chr = "L1")
plot.info(x=ut)
plot.map(x = ut)
plot.qtl(x=ut)
qtl::plotplot.cross(ut)
