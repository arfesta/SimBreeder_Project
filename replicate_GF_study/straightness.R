# Run this on Domino ####
source('/mnt/simulator/sim_phenos.R')
source('/mnt/simulator/create_cross_design.R')
source('/mnt/simulator/make_crosses.R')
source('/mnt/simulator/calc_TGV.R')
source('/mnt/simulator/extract_selections.R')
source('/mnt/simulator/create_map.R')
source('/mnt/simulator/create_parents.R')
source('/mnt/simulator/OP_testing.R')

## Striaghtness ####
D=.5
minor=-1
#set.seed(435634)
set.seed(3345345)
map <- create_genetic_map(num.chromos = 12,map.length = 1800,num.markers = 120,total.qtl = 2640,num.snpqtl = 1960,
                          distribute.loci ="even",marker.distribution = "equally-spaced",snp.qtl.maf = c(0.01,0.011))
parents <- create_parents(map.info = map,num.parents = 50,max.delt.allele = 14,heterozygous.markers = F)
parent.tgv <- calc_TGV(geno.info = parents,map.info = map,founder = T,A = 1,a = minor,dom.coeff = D)
set.seed(435634)
parent.phenos <- sim_phenos(TGV.object = parent.tgv,h2 = .3)
first.gen.cross <- create_cross_design(parent.info = parents,mating.design = "cross.file.input",cross.file = "/mnt/graham_study/g.study.1stgen.txt",generation = 1)
first.prog <- make_crosses(parent.info = parents,map.info = map,cross.design = first.gen.cross,num.cores = parallel::detectCores())
first.prog.tgv <- calc_TGV(geno.info = first.prog,map.info  = map,cross.design = first.gen.cross,A = 1,a = minor,dom.coeff = D)
first.prog.phenos <- sim_phenos(TGV.object = first.prog.tgv,h2 = .3)

prog.ex1 <- extract_selections(among.family.selection = "Phenotype",relationship.matrix.type = "pedigree",map.info = map,
                               cross.design =  first.gen.cross,past.tgv = parent.tgv,past.phenos = parent.phenos,parent.info = parents,
                               progeny.info = first.prog,progeny.TGV = first.prog.tgv,progeny.phenos = first.prog.phenos,
                               num.selections.within.family = 1, 
                               num.selections.among.family = 40)
selections <- sort(as.numeric(names(prog.ex1$selection.phenos)))
second.gen.cross <- read.table("/mnt/graham_study/g.2nd.gen.ped2.txt")
second.gen.cross <- cbind(selections[second.gen.cross[,1]],selections[second.gen.cross[,2]],second.gen.cross[,3])
write.table(second.gen.cross,"./g.2nd.gen.pedv2.txt",row.names = F,col.names = F)
second.gen.cross <-create_cross_design(parent.info = prog.ex1,mating.design = "cross.file.input",cross.file = "./g.2nd.gen.pedv2.txt",generation = 2) 
second.prog <- make_crosses(parent.info = prog.ex1,map.info = map,cross.design = second.gen.cross,num.cores = parallel::detectCores())
second.prog.tgv <- calc_TGV(geno.info = second.prog,map.info = map,cross.design = second.gen.cross,A = 1,a = minor,dom.coeff = D)
evar <- var(parent.phenos$phenos) - var(parent.phenos$genetic.values)
second.prog.phenos <- sim_phenos(TGV.object = second.prog.tgv,E.var = evar)

#### make plot ####
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
  Self.phenos <- mean.phenos[Self]
  FS.phenos <- mean.phenos[FS]
  HS.phenos <- mean.phenos[HS]
  outcross <- rep(0,10)
  half<- rep(.125,10)
  full<- rep(.25,10)
  self<- rep(.5,10)
  coastal.ht <- data.frame(par=rep(1:10,times=4),x=c(outcross,half,full,self),y=c(OP.phenos,HS.phenos,FS.phenos,Self.phenos))
  # Create Line Chart
  a<-.3
  b<-.6
  range01 <- function(x){(((b-a)*(x-min(x)))/(max(x)-min(x)))+a}
  cy <- range01(coastal.ht$y)
  coastal.ht$y <- cy
  
  
  coastal.ht$par <- as.factor(coastal.ht$par)
  coastal.ht$Parents <- coastal.ht$par
 strt <-  ggplot(coastal.ht, aes(x=x, y=y,shape=Parents,group=Parents)) + 
    geom_line() + 
    geom_point() +
    scale_shape_manual(values = c(0,9, 16, 3,11,8,4,3,2,6)) +
    scale_y_continuous(breaks=seq(0, 1, .2),limits=c(0,1)) +
    scale_x_continuous(breaks=c(0,.125,.25,.5),limits=c(0,.5)) +
    xlab("Inbreeding Coefficient") +
    ylab("Probability of Straightness") +
    labs(title = "Straightness") +
    theme_pubr()
  
  
  # convert factor to numeric for convenience 
  ntrees <- max(coastal.ht$par)
  
  # get the range for the x and y axis 
  xrange <- c(0.0,0.5)
  yrange <- c(0,1)
  
  # set up the plot 
  plot(xrange, yrange, axes=F,type="n", xlab="Inbreeding Coefficient",
       ylab="Probability of Straightness", bty='L', main= "Straighness") 
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
  yticks <- seq(0,1,.2)
  axis(2, at = yticks, labels = yticks, col.axis="black",las=1)
  axis(1, at = xticks, labels = xticks, col.axis="black",las=1,lwd.ticks = .5)
  par(xpd=TRUE)
  legend()
}

op_mean <- mean(coastal.ht$y[1:10])
hs_mean <- mean(coastal.ht$y[11:20])
fs_mean <- mean(coastal.ht$y[21:30])
s_mean <- mean(coastal.ht$y[31:40])
i_levels <- c(0,.125,.25,.5)
i_means <- c(op_mean,hs_mean,fs_mean,s_mean)
#summary(lm(coastal.ht$y ~ 0 +coastal.ht$x))
summary(lm(i_means ~ i_levels))

1-(hs_mean/op_mean)
1-(fs_mean/op_mean)
1-(s_mean/op_mean)

t.test(coastal.ht$y[11:20],coastal.ht$y[31:40],paired = T)

op_hs <- coastal.ht$y[1:10] - coastal.ht$y[11:20]
hs_fs <- coastal.ht$y[11:20] - coastal.ht$y[21:30]
fs_s <- coastal.ht$y[21:30] - coastal.ht$y[31:40]
cbind(op_hs,hs_fs,fs_s)


op_hs <- coastal.ht$y[1:10] - coastal.ht$y[11:20]
hs_fs <- coastal.ht$y[1:10] - coastal.ht$y[21:30]
fs_s <- coastal.ht$y[1:10] - coastal.ht$y[31:40]
cbind(op_hs,hs_fs,fs_s)

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
