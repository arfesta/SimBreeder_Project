# Run this on Domino ####
source('/mnt/simulator/sim_phenos.R')
source('/mnt/simulator/create_cross_design.R')
source('/mnt/simulator/make_crosses.R')
source('/mnt/simulator/calc_TGV.R')
source('/mnt/extract_selections.new2.R')
source('/mnt/simulator/create_map.R')
source('/mnt/simulator/create_parents.R')
source('/mnt/simulator/OP_testing.R')

## Height
D=1
minor=-100
set.seed(435634)
map <- create_genetic_map(num.chromos = 12,map.length = 1800,num.markers = 540,total.qtl = 2520,num.snpqtl = 1960,
                          distribute.loci ="even",marker.distribution = "equally-spaced",snp.qtl.maf = c(0.01,0.02))
parents <- create_parents(map.info = map,num.parents = 50,max.delt.allele = 14,heterozygous.markers = F)
parent.phenos <- sim_phenos(TGV.object = parent.tgv,h2 = .3)
first.gen.cross <- create_cross_design(parent.info = parents,mating.design = "cross.file.input",cross.file = "./graham_study/g.study.1stgen.txt",generation = 1)
first.prog <- make_crosses(parent.info = parents,map.info = map,cross.design = first.gen.cross,num.cores = 30)
first.prog.tgv <- calc_TGV(geno.info = first.prog,map.info  = map,cross.design = first.gen.cross,A = 1,a = minor,dom.coeff = D)
first.prog.phenos <- sim_phenos(TGV.object = first.prog.tgv,h2 = .3)

prog.ex1 <- extract_selections(among.family.selection = "Phenotype",map.info = map,cross.design =  first.gen.cross,past.tgv = parent.tgv,past.phenos = parent.phenos,parent.info = parents,
                               progeny.info = first.prog,progeny.TGV = first.prog.tgv,progeny.phenos = first.prog.phenos,
                               num.selections.within.family = 1, 
                               num.selections.among.family = 40)
selections <- sort(as.numeric(names(prog.ex1$selection.phenos)))
second.gen.cross <- read.table("./graham_study/g.2nd.gen.ped2.txt")
second.gen.cross <- cbind(selections[second.gen.cross[,1]],selections[second.gen.cross[,2]],second.gen.cross[,3])
write.table(second.gen.cross,"./g.2nd.gen.ped.txt",row.names = F,col.names = F)
second.gen.cross <-create_cross_design(parent.info = prog.ex1,mating.design = "cross.file.input",cross.file = "g.2nd.gen.ped.txt",generation = 2) 
second.prog <- make_crosses(parent.info = prog.ex1,map.info = map,cross.design = second.gen.cross,num.cores = 30)
second.prog.tgv <- calc_TGV(geno.info = second.prog,map.info = map,cross.design = second.gen.cross,A = 1,a = minor,dom.coeff = D)
second.prog.phenos <- sim_phenos(TGV.object = second.prog.tgv,h2 = .3)

{
  first.prog <- seq(1,2000,50)
  length(first.prog)
  OP <- c(seq(1,length(first.prog),by = 4))
  Self <- c(seq(2,length(first.prog),by = 4))
  FS <- c(seq(3,length(first.prog),by = 4))
  HS <- c(seq(4,length(first.prog),by = 4))
  mean.phenos <- vector()
  for(each in 1:length(first.prog)) { mean.phenos <- c(mean.phenos,mean(second.prog.phenos$phenos[first.prog[each]:(first.prog[each]+49)]))}
  
  
  OP.phenos <- (mean.phenos[OP])
  hist(OP.phenos)
  Self.phenos <- mean.phenos[Self]
  hist(Self.phenos)
  FS.phenos <- mean.phenos[FS]
  HS.phenos <- mean.phenos[HS]
  outcross <- rep(0,10)
  half<- rep(.125,10)
  full<- rep(.25,10)
  self<- rep(.5,10)
  coastal.ht <- data.frame(par=rep(1:10,times=4),x=c(outcross,half,full,self),y=c(OP.phenos,HS.phenos,FS.phenos,Self.phenos))
  # Create Line Chart
  a<-6.75
  b<-11
  range01 <- function(x){(((b-a)*(x-min(x)))/(max(x)-min(x)))+a}
  cy <- range01(coastal.ht$y)
  coastal.ht$y <- cy
  # convert factor to numeric for convenience 
  ntrees <- max(coastal.ht$par)
  
  # get the range for the x and y axis 
  xrange <- c(0.0,0.5)
  yrange <- c(6.5,11)
  
  # set up the plot 
  plot(xrange, yrange, axes=F,type="n", xlab="Inbreeding Coefficient",
       ylab="Height (m)", bty='L', main= "Height") 
  color <- rainbow(ntrees)
  color[3] <- colors()[641]
  linetype <- c(1:ntrees) 
  plotchar <- seq(1,1+ntrees,1)
  
  # add lines 
  for (i in 1:ntrees) { 
    tree <- subset(coastal.ht, par==i) 
    lines(tree$x, tree$y, type="b", lwd=2,
          lty=linetype[i], col=color[i], pch=plotchar[i], cex=.5) 
  } 
  
  xticks <- c(0,.125,.25,.5)
  yticks <- seq(7,11,1)
  axis(2, at = yticks, labels = yticks, col.axis="black",las=1)
  axis(1, at = xticks, labels = xticks, col.axis="black",las=1,lwd.ticks = .5)
  par(xpd=TRUE)
}

library(heritability)
h2 <- data.frame(second.prog.tgv$markers.matrix,stringsAsFactors = F)
h2 <- apply(h2,MARGIN = 2, function(x) as.numeric(x))

library(rrBLUP)
out <- A.mat((h2)); colnames(out) <- names(second.prog.phenos$phenos) ; rownames(out) <- names(second.prog.phenos$phenos)

sample.num <- sample(x = 1:ncol(out),size = 100,replace = F)
sample.num <- 1:2000
out.test <- out[sample.num,sample.num]; out.phenos <- second.prog.phenos$phenos[sample.num]

the.out.data <- marker_h2(data.vector = out.phenos,geno.vector = rownames(out.test),K = out.test,fix.h2 = T,h2 = .22)
the.out.data
