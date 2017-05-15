cor.mat <- matrix(0,nrow=5824,ncol=5824)
for(i in 1:120) {
  these <- unique(allele1[i,])
  for(all in 1:length(these)) {
  posit <- which(allele1[i,] == these[all])
  for(each in 1:length(posit)){ 
    fill <- posit[each]
    others <- posit[-c(1:each)]
    cor.mat[fill,others] <- cor.mat[fill,others] + 1/240
  }}}

cor.mat2 <- matrix(0,nrow=5824,ncol=5824)
for(i in 1:120) {
  these <- unique(allele2[i,])
  for(all in 1:length(these)) {
    posit <- which(allele2[i,] == these[all])
    for(each in 1:length(posit)){ 
      fill <- posit[each]
      others <- posit[-c(1:each)]
      cor.mat2[fill,others] <- cor.mat2[fill,others] + 1/240
    }}}

the.cor <- cor.mat + cor.mat2
rm(cor.mat,cor.mat2); gc()
