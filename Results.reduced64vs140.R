###AAA.64####
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop1/AAA.64.OP.pop1.1.5.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop1/AAA.64.OP.pop1.6.10.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop1/AAA.64.OP.pop1.11.20.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop1/AAA.64.OP.pop1.21.30.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop3/AAApop3.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop2/AAApop2.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop4/AAApop4.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop5/AAApop5.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop6/AAApop6.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop7/AAApop7.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop8/AAApop8.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop9/AAApop9.RDA")
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.64/pop10/AAApop10.RDA")
AAA.64 <- c(AAA.64.OP.pop1.1.5,AAA.64.OP.pop1.6.10,AAA.64.OP.pop1.11.20,AAA.64.OP.pop1.21.30,AAApop2,AAApop3,AAApop4,AAApop5,
         AAApop6, AAApop7, AAApop8,AAApop9,AAApop10)
rm(list=setdiff(ls(), "AAA.64"))

####AGA.64#####
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop1/AGA.64.OP.pop1.1.5.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop1/AGA.64.OP.pop1.6.10.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop1/AGA.64.OP.pop1.11.15.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop1/AGA.64.OP.pop1.16.20.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop1/AGA.64.OP.pop1.21.25.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop1/AGA.64.OP.pop1.26.30.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop2/AGApop2.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop3/AGApop3.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop8/AGApop8.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop9/AGApop9.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop4/AGApop4.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop5/AGApop5.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop6/AGApop6.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop7/AGApop7.RDA")
load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.64/pop10/AGApop10.RDA")

AGA.64 <- c(AGA.64.OP.pop1.1.5,AGA.64.OP.pop1.6.10,AGA.64.OP.pop1.11.15,AGA.64.OP.pop1.16.20, AGA.64.OP.pop1.21.25, 
         AGA.64.OP.pop1.26.30,AGApop2,AGApop3, AGApop8, AGApop9, AGApop10,
         AGApop4, AGApop5, AGApop6, AGApop7)

rm(list=setdiff(ls(), c("AAA.64","AGA.64")))

###AAA.140####
load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.140/pop10/AAApop10.RDA")

AAA.140 <- c(AAApop10)
rm(list=setdiff(ls(), c("AAA.140","AAA.64","AGA.64")))

####AGA.140#####

load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.140/pop10/AGApop10.RDA")

AGA.140 <- c(AGApop10)

rm(list=setdiff(ls(), c("AAA.140","AGA.140","AAA.64","AGA.64")))


####AAA.64 results####
n = length(AAA.64)
AAA.64.g.gain <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.64.g.gain[i,] <- AAA.64[[i]]$g}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.64.g.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.64.g.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.64.g.gain,2,mean)

AAA.64.gv.gain <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.64.gv.gain[i,] <- AAA.64[[i]]$gv}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.64.gv.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.64.gv.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.64.gv.gain,2,mean)

AAA.64.delt.allele <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.64.delt.allele[i,] <- AAA.64[[i]]$da}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.64.delt.allele[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.64.delt.allele[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.64.delt.allele,2,mean)


AAA.64.ct  <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.64.ct[i,] <- AAA.64[[i]]$ct}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.64.ct[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.64.ct[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.64.ct,2,mean)

##AGA.64 RESULTS####
n = length(AGA.64)
AGA.64.g.gain <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AGA.64.g.gain[i,] <- AGA.64[[i]]$g}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AGA.64.g.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AGA.64.g.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AGA.64.g.gain,2,mean)

AGA.64.gv.gain <- matrix(,ncol=10,nrow=n)

for(i in 1:n){AGA.64.gv.gain[i,] <- AGA.64[[i]]$gv}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AGA.64.gv.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AGA.64.gv.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AGA.64.gv.gain,2,mean)

AGA.64.delt.allele <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AGA.64.delt.allele[i,] <- AGA.64[[i]]$da}

AGA.64.ct <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AGA.64.ct[i,] <- AGA.64[[i]]$ct}

####AAA.140 results####
n = length(AAA.140)
AAA.140.g.gain <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.140.g.gain[i,] <- AAA.140[[i]]$g}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.140.g.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.140.g.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.140.g.gain,2,mean)

AAA.140.gv.gain <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.140.gv.gain[i,] <- AAA.140[[i]]$gv}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.140.gv.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.140.gv.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.140.gv.gain,2,mean)

AAA.140.delt.allele <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.140.delt.allele[i,] <- AAA.140[[i]]$da}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.140.delt.allele[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.140.delt.allele[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.140.delt.allele,2,mean)


AAA.140.ct  <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.140.ct[i,] <- AAA.140[[i]]$ct}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.140.ct[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.140.ct[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.140.ct,2,mean)

##AGA.140 RESULTS####
n = length(AGA.140)
AGA.140.g.gain <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AGA.140.g.gain[i,] <- AGA.140[[i]]$g}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AGA.140.g.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AGA.140.g.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AGA.140.g.gain,2,mean)

AGA.140.gv.gain <- matrix(,ncol=10,nrow=n)

for(i in 1:n){AGA.140.gv.gain[i,] <- AGA.140[[i]]$gv}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AGA.140.gv.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AGA.140.gv.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AGA.140.gv.gain,2,mean)

AGA.140.delt.allele <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AGA.140.delt.allele[i,] <- AGA.140[[i]]$da}

AGA.140.ct <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AGA.140.ct[i,] <- AGA.140[[i]]$ct}


###Genetic variance####
G.var.mat <- rbind(apply(AAA.64.gv.gain,2,mean),apply(AAA.140.gv.gain,2,mean),
                   apply(AGA.64.gv.gain,2,mean),apply(AGA.140.gv.gain,2,mean))
stdev <- rbind(apply(AAA.64.gv.gain,2,sd),apply(AAA.140.gv.gain,2,sd),
               apply(AGA.64.gv.gain,2,sd),apply(AGA.140.gv.gain,2,sd))
G.var.mat <- t(G.var.mat)
stdev <- t(stdev)
low <- G.var.mat - stdev  ##  lower bound
up <- G.var.mat + stdev  ## upper bound

x <- seq(1:10)  ## x-locations to plot against
## your colour setting; should have `ncol(mat)` colours
## as an example I just use `cols = 1:ncol(mat)`
#cols <- 1:ncol((G.var.mat))
cols <- c(2,2,4,4)
## plot each column of `mat` one by one (set y-axis limit appropriately)
#matplot(x, (G.var.mat), col = cols, pch = 1:2, type = "o", ylim = c(min(low), max(up)),xlab = "Generation",ylab="Genetic Variance")
matplot(x, (G.var.mat), col = cols, pch = 1:2, type = "o", xlab = "Generation",ylab="Genetic Variance")
xx <- rep.int(x, ncol((G.var.mat)))  ## recycle `x` for each column of `mat`
repcols <- rep(cols, each = nrow((G.var.mat)))  ## recycle `col` for each row of `mat`
## adding error bars using vectorization power of `arrow`
#arrows(xx, low, xx, up, col = repcols, angle = 90, length = 0.03, code = 3)

###Genetic gain####
G.gain.mat <- rbind(apply(AAA.64.g.gain,2,mean),apply(AAA.140.g.gain,2,mean),
                    apply(AGA.64.g.gain,2,mean),apply(AGA.140.g.gain,2,mean))
stdev <- rbind(apply(AAA.64.g.gain,2,sd),apply(AAA.140.g.gain,2,sd),
               apply(AGA.64.g.gain,2,sd),apply(AGA.140.g.gain,2,sd))
G.gain.mat <- t(G.gain.mat)
stdev <- t(stdev)
low <- G.gain.mat - stdev  ##  lower bound
up <- G.gain.mat + stdev  ## upper bound

x <- seq(1:10)  ## x-locations to plot against
## your colour setting; should have `ncol(mat)` colours
## as an example I just use `cols = 1:ncol(mat)`
cols <- c(2,2,4,4)
## plot each column of `mat` one by one (set y-axis limit appropriately)
matplot(x, (G.gain.mat), col = cols, pch = 1:4, type = "o", ylim = c(min(low), max(up)),xlab = "Generation",ylab="Genetic Gain")
xx <- rep.int(x, ncol((G.gain.mat)))  ## recycle `x` for each column of `mat`
repcols <- rep(cols, each = nrow((G.gain.mat)))  ## recycle `col` for each row of `mat`
## adding error bars using vectorization power of `arrow`
#arrows(xx, low, xx, up, col = repcols, angle = 90, length = 0.03, code = 3)

###Delt Allele####
da.mat <- rbind(apply(AAA.64.delt.allele,2,mean),apply(AAA.140.delt.allele,2,mean),
                apply(AGA.64.delt.allele,2,mean),apply(AGA.140.delt.allele,2,mean))
stdev <- rbind(apply(AAA.64.delt.allele,2,sd),apply(AAA.140.delt.allele,2,sd),
               apply(AGA.64.delt.allele,2,sd),apply(AGA.140.delt.allele,2,sd))
da.mat <- t(da.mat)
stdev <- t(stdev)
low <- da.mat - stdev  ##  lower bound
up <- da.mat + stdev  ## upper bound

x <- seq(1:10)  ## x-locations to plot against
## your colour setting; should have `ncol(mat)` colours
## as an example I just use `cols = 1:ncol(mat)`
cols <- c(2,2,4,4)
## plot each column of `mat` one by one (set y-axis limit appropriately)
matplot(x, (da.mat), col = cols, pch = 1:4, type = "o", ylim = c(min(low), max(up)),xlab = "Generation",ylab="Avg Delt alleles")
xx <- rep.int(x, ncol((da.mat)))  ## recycle `x` for each column of `mat`
repcols <- rep(cols, each = nrow((da.mat)))  ## recycle `col` for each row of `mat`
## adding error bars using vectorization power of `arrow`
#arrows(xx, low, xx, up, col = repcols, angle = 90, length = 0.03, code = 3)

matplot(x, (G.var.mat), col = cols, pch = c(1,2,1,2), type = "o", xlab = "Generation",ylab="Genetic Variance")
legend("topright", inset=.05, legend=c("AAA 64", "AAA 140", "AGA 64", "AGA 140"), pch = c(1,2,1,2), col=cols, horiz=T)
matplot(x, (G.gain.mat), col = cols, pch = c(1,2,1,2), type = "o", xlab = "Generation",ylab="Genetic Gain")
legend("topleft", inset=.05, legend=c("AAA 64", "AAA 140", "AGA 64", "AGA 140"), pch = c(1,2,1,2), col=cols, horiz=F)
matplot(x, (da.mat), col = cols, pch = c(1,2,1,2), type = "o", xlab = "Generation",ylab="Avg Delt alleles")
legend("bottomleft", inset=.05, legend=c("AAA 64", "AAA 140", "AGA 64", "AGA 140"), pch = c(1,2,1,2), col=cols, horiz=F)

###CT Allele####
ct.mat <- rbind(apply(AAA.ct,2,mean),apply(AGA.ct,2,mean))
stdev <- rbind(apply(AAA.ct,2,sd),apply(AGA.ct,2,sd))
ct.mat <- t(ct.mat)
stdev <- t(stdev)
low <- ct.mat - stdev  ##  lower bound
up <- ct.mat + stdev  ## upper bound

x <- seq(1:10)  ## x-locations to plot against
## your colour setting; should have `ncol(mat)` colours
## as an example I just use `cols = 1:ncol(mat)`
cols <- 1:ncol((ct.mat))
## plot each column of `mat` one by one (set y-axis limit appropriately)
matplot(x, (ct.mat), col = cols, pch = 1:2, type = "o", ylim = c(min(low), max(up)),xlab = "Generation",ylab="Co-ancestry Threshold")
xx <- rep.int(x, ncol((ct.mat)))  ## recycle `x` for each column of `mat`
repcols <- rep(cols, each = nrow((ct.mat)))  ## recycle `col` for each row of `mat`
## adding error bars using vectorization power of `arrow`
arrows(xx, low, xx, up, col = repcols, angle = 90, length = 0.03, code = 3)

###Genetic variance against gain####
G.var.mat <- rbind(apply(AAA.64.gv.gain,2,mean),apply(AAA.140.gv.gain,2,mean),
                   apply(AGA.64.gv.gain,2,mean),apply(AGA.140.gv.gain,2,mean))

G.var.mat <- t(G.var.mat)


x <- G.gain.mat ## x-locations to plot against

low <- G.gain.mat  ##  lower bound
up <- G.gain.mat  ## upper bound
## your colour setting; should have `ncol(mat)` colours
## as an example I just use `cols = 1:ncol(mat)`
#cols <- 1:ncol((G.var.mat))
cols <- c(2,2,4,4)
## plot each column of `mat` one by one (set y-axis limit appropriately)
matplot(x, (G.var.mat) , col = cols, pch = 1:4, type = "o",xlab = "Generation",ylab="Genetic Variance")
xx <- rep.int(x, ncol((G.var.mat)))  ## recycle `x` for each column of `mat`
repcols <- rep(cols, each = nrow((G.var.mat)))  ## recycle `col` for each row of `mat`
## adding error bars using vectorization power of `arrow`
#arrows(xx, low, xx, up, col = repcols, angle = 90, length = 0.03, code = 3)