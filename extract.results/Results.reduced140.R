###AAA####
#file_names=as.list(dir(path ="/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.140/pop10" ,all.files = T,no.. = T,full.names = T))
#lapply(file_names,load,.GlobalEnv); rm(file_names)
#AAApop10 <- lapply(ls(), function(x) if (class(get(x)) == "list") get(x))
#rm(list=setdiff(ls(), "AAApop10"))
#save(AAApop10,file="/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.140/pop10/AAApop10.RDA")


load("/media/titan/ARF/Simulator.Package/AAA.ST/results/reduced.140/pop10/AAApop10.RDA")

AAA <- c(AAApop10)
rm(list=setdiff(ls(), "AAA"))

####AGA#####
#file_names=as.list(dir(path ="/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.140/pop10" ,all.files = T,no.. = T,full.names = T))
#lapply(file_names,load,.GlobalEnv); rm(file_names)
#AGApop10 <- lapply(ls(), function(x) if (class(get(x)) == "list") get(x))
#rm(list=setdiff(ls(), "AGApop10"))
#save(AGApop10,file="/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.140/pop10/AGApop10.RDA")

load("/media/titan/ARF/Simulator.Package/AGA.ST/results/reduced.140/pop10/AGApop10.RDA")

AGA <- c(AGApop10)

rm(list=setdiff(ls(), c("AAA","AGA")))


####AAA results####
n = length(AAA)
AAA.g.gain <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.g.gain[i,] <- AAA[[i]]$g}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.g.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.g.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.g.gain,2,mean)

AAA.gv.gain <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.gv.gain[i,] <- AAA[[i]]$gv}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.gv.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.gv.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.gv.gain,2,mean)

AAA.delt.allele <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.delt.allele[i,] <- AAA[[i]]$da}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.delt.allele[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.delt.allele[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.delt.allele,2,mean)


AAA.ct  <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AAA.ct[i,] <- AAA[[i]]$ct}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AAA.ct[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AAA.ct[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AAA.ct,2,mean)

##aga RESULTS####
n = length(AGA)
AGA.g.gain <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AGA.g.gain[i,] <- AGA[[i]]$g}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AGA.g.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AGA.g.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AGA.g.gain,2,mean)

AGA.gv.gain <- matrix(,ncol=10,nrow=n)

for(i in 1:n){AGA.gv.gain[i,] <- AGA[[i]]$gv}
upper.conf <- vector(); lower.conf <- vector()
for(i in 1:10){
  upper.conf <- c(upper.conf,t.test(AGA.gv.gain[,i])$conf.int[2])
  lower.conf <- c(lower.conf,t.test(AGA.gv.gain[,i])$conf.int[1])
}
(upper.conf-lower.conf)/apply(AGA.gv.gain,2,mean)

AGA.delt.allele <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AGA.delt.allele[i,] <- AGA[[i]]$da}

AGA.ct <- matrix(,ncol=10,nrow=n)
for(i in 1:n){AGA.ct[i,] <- AGA[[i]]$ct}

###Genetic variance####
G.var.mat <- rbind(apply(AAA.gv.gain,2,mean),apply(AGA.gv.gain,2,mean))
stdev <- rbind(apply(AAA.gv.gain,2,sd),apply(AGA.gv.gain,2,sd))
G.var.mat <- t(G.var.mat)
stdev <- t(stdev)
low <- G.var.mat - stdev  ##  lower bound
up <- G.var.mat + stdev  ## upper bound

x <- seq(1:10)  ## x-locations to plot against
## your colour setting; should have `ncol(mat)` colours
## as an example I just use `cols = 1:ncol(mat)`
cols <- 1:ncol((G.var.mat))
## plot each column of `mat` one by one (set y-axis limit appropriately)
matplot(x, (G.var.mat), col = cols, pch = 1:2, type = "o", ylim = c(min(low), max(up)),xlab = "Generation",ylab="Genetic Variance")
xx <- rep.int(x, ncol((G.var.mat)))  ## recycle `x` for each column of `mat`
repcols <- rep(cols, each = nrow((G.var.mat)))  ## recycle `col` for each row of `mat`
## adding error bars using vectorization power of `arrow`
arrows(xx, low, xx, up, col = repcols, angle = 90, length = 0.03, code = 3)


###Genetic gain####
G.gain.mat <- rbind(apply(AAA.g.gain,2,mean),apply(AGA.g.gain,2,mean))
stdev <- rbind(apply(AAA.g.gain,2,sd),apply(AGA.g.gain,2,sd))
G.gain.mat <- t(G.gain.mat)
stdev <- t(stdev)
low <- G.gain.mat - stdev  ##  lower bound
up <- G.gain.mat + stdev  ## upper bound

x <- seq(1:10)  ## x-locations to plot against
## your colour setting; should have `ncol(mat)` colours
## as an example I just use `cols = 1:ncol(mat)`
cols <- 1:ncol((G.gain.mat))
## plot each column of `mat` one by one (set y-axis limit appropriately)
matplot(x, (G.gain.mat), col = cols, pch = 1:2, type = "o", ylim = c(min(low), max(up)),xlab = "Generation",ylab="Genetic Gain")
xx <- rep.int(x, ncol((G.gain.mat)))  ## recycle `x` for each column of `mat`
repcols <- rep(cols, each = nrow((G.gain.mat)))  ## recycle `col` for each row of `mat`
## adding error bars using vectorization power of `arrow`
arrows(xx, low, xx, up, col = repcols, angle = 90, length = 0.03, code = 3)

###Delt Allele####
da.mat <- rbind(apply(AAA.delt.allele,2,mean),apply(AGA.delt.allele,2,mean))
stdev <- rbind(apply(AAA.delt.allele,2,sd),apply(AGA.delt.allele,2,sd))
da.mat <- t(da.mat)
stdev <- t(stdev)
low <- da.mat - stdev  ##  lower bound
up <- da.mat + stdev  ## upper bound

x <- seq(1:10)  ## x-locations to plot against
## your colour setting; should have `ncol(mat)` colours
## as an example I just use `cols = 1:ncol(mat)`
cols <- 1:ncol((da.mat))
## plot each column of `mat` one by one (set y-axis limit appropriately)
matplot(x, (da.mat), col = cols, pch = 1:2, type = "o", ylim = c(min(low), max(up)),xlab = "Generation",ylab="Avg Delt alleles")
xx <- rep.int(x, ncol((da.mat)))  ## recycle `x` for each column of `mat`
repcols <- rep(cols, each = nrow((da.mat)))  ## recycle `col` for each row of `mat`
## adding error bars using vectorization power of `arrow`
arrows(xx, low, xx, up, col = repcols, angle = 90, length = 0.03, code = 3)

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
