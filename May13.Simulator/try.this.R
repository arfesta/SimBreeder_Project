

tot <- vector()
for(each.row in 1:120){
  tot <- c(tot,length(which(p1[each.row,] %in% p2[each.row,]))) }

sum(tot)/240


g.rel.mat <- matrix(0,nrow=ncol(prog.1),ncol=ncol(prog.1))
for(i in 1:120) {
  these <- unique(c(prog.1[i,],prog.2[i,]))
  for(all in 1:length(these)) {
    posit <- (which(prog.1[i,] %in% these[all] | prog.2[i,] %in% these[all]))
    for(each in 1:length(posit)){ 
      fill <- posit[each]
      others <- posit[-c(1:each)]
      g.rel.mat[fill,others] <- g.rel.mat[fill,others] + 1/120
    }}}


these <- c(prog.1[1,1],prog.2[1,1])
prog.1 <- progmarkers[,,1]
prog.2 <- progmarkers[,,2]


final.mat <- matrix(0,nrow = ncol(prog.1),ncol = ncol(prog.1))
for (this.prog in 1:ncol(prog.1)) 
fin <- mclapply(1:ncol(prog.1),function(x){  
  l <- matrix(0,nrow=120,ncol=ncol(prog.1))
  ty <- prog.1[,x]
  ty2 <- prog.2[,x]
for(each.row in 1:120){
  these <- (c(ty[each.row],ty2[each.row]))
l1 <- c(which(prog.1[each.row,] %in% these[1]))
l2 <- c(which(prog.1[each.row,] %in% these[2]))
l3 <- c(which(prog.2[each.row,] %in% these[1]))
l4 <- c(which(prog.2[each.row,] %in% these[2]))

if(length(unique(these)) < 2){
yo <-  c(table(c(l1,l2,l3,l4))/4) 
l[each.row,as.numeric(names(yo))] <- yo
} else {
yo <- c(table(c(l1,l2,l3,l4))/2)
l[each.row,as.numeric(names(yo))] <- yo}}

 colMeans(l)},mc.cores=24)

final.mat <- matrix(0,nrow = ncol(prog.1),ncol = ncol(prog.1))
for(each.list in 1:length(fin)){
  final.mat[each.list,] <- fin[[each.list]]}