the.files <- list.files("./140/AAA",full.names = T)
lapply(the.files,load,.GlobalEnv)
rm(the.files)

AAA <- lapply(ls(), get)
rm(list=setdiff(ls(), "AAA"))

the.files <- list.files("./140/AGA",full.names = T)
lapply(the.files,load,.GlobalEnv)
rm(the.files)

AGA <- lapply(setdiff(ls(),"AAA"), get)
rm(list=setdiff(ls(), c("AAA","AGA")))


 for(x in 1:length(AGA)) {
     this <- AGA[[x]]$g
     if (x==1) { AGA.g <- as.matrix(this)} else { AGA.g <- rbind(AGA.g,this)} 
    
 } ; rm(this,x)
AGA.g = apply(AGA.g,2,mean)


for(x in 1:length(AAA)) {
  this <- (AAA[[x]]$g)
 if (x==1) { AAA.g <- as.matrix(this)} else { AAA.g <- rbind(AAA.g,this)} 
  
} ; rm(this,x)
AAA.g = apply(AAA.g,2,mean)
plot(1:10,ylim = c(mi))
plot(1:10,AGA.g,ylim=c(min(AAA.g),max(AGA.g)))
lines(1:10,AAA.g)
