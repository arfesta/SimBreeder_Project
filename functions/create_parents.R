create_parents <- function(map.info, num.parents, max.delt.allele, QTL.sd = 0.25){
  
  num.rQTL <- length(map.info$rQTL.loci)  # the number of additive qtl (rQTL)
  num.markers <- length(map.info$available.Markers)  # the number of potential markers
  
  #  Create Dimension names for 2 arrays
  par.IDs <- 1:num.parents; allele.IDs <- c("a1","a2")
  
  # Parents array holds both alleles for all parents
  # QTLSNP array holds the alleles for the loci of all parents under dominance control
  parents <- array(rep(NA,map.info$total.loci.num*num.parents*2), dim=c(map.info$total.loci.num,num.parents,2),dimnames=list(map.info$genetic.map$loci,par.IDs,allele.IDs))
  QTLSNP.array <- array(0, dim=c(map.info$total.SNPQTL.num,num.parents,2),dimnames=list(map.info$genetic.map$loci[map.info$QTLSNP.loci],par.IDs,allele.IDs))
  
  #  Sample rQTL values from two normal distributions with mean 0 and std deviation half of that specified by QTLsd variable
  #  Use 'ceiling' function for one distribution, 'floor' function for other, then add them together to get centered distribution
  dist1 <- ceiling(rnorm((num.rQTL*2), 0, (QTL.sd*(sqrt(2)/2)))); dist2 <- floor(rnorm((num.rQTL*2), 0, (QTL.sd*(sqrt(2)/2))))
  QTL.alleles <- dist1 + dist2  # Vector contains the values for the rQTL alleles
  
  # Select character values of dominant and recessive alleles to markers
  if(length(which(map.info$genetic.map$types == "fi")) != 0 ){
    d.fi.marker.genotypes <- vector(); r.fi.marker.genotypes <- vector(); 
    for(i in 1:26){
      d.fi.marker.genotypes <- c(d.fi.marker.genotypes,paste(LETTERS[i],LETTERS, sep=""))
      r.fi.marker.genotypes <-  c(r.fi.marker.genotypes,paste(letters[i],letters, sep=""))
    }
  } 
  if(length(which(map.info$genetic.map$types == "m")) != 0) {
    d.marker.genotypes <- "G"
    r.marker.genotypes <- "g"
  }
  
  #  Create empty vectors to hold either allele1 or allele2  for each parent
  
  
  # For each parent assign them alleles
  while(anyNA(parents)){
    if(max.delt.allele == 0 ){
      for(par in 1:num.parents){
        allele1<-rep(NA, map.info$total.loci.num); allele2<-rep(NA, map.info$total.loci.num)
        maf.test <- runif(n = map.info$total.SNPQTL.num,min = 0,max = 1)
        minor1 <- map.info$QTLSNP.loci[which(map.info$genetic.map[map.info$QTLSNP.loci,"MAF"]  > maf.test)]
        maf.test <- runif(n = map.info$total.SNPQTL.num,min = 0,max = 1)
        minor2 <- map.info$QTLSNP.loci[which(map.info$genetic.map[map.info$QTLSNP.loci,"MAF"]  > maf.test)]
        major1 <- map.info$QTLSNP.loci[which(!map.info$QTLSNP.loci %in% minor1)]
        major2 <- map.info$QTLSNP.loci[which(!map.info$QTLSNP.loci %in% minor2)]
        allele1[minor1]<-"c"            # minor allele is always "c", major allele is always "a"
        allele1[major1] <-"a"
        allele2[minor2]<-"c"
        allele2[major2]<-"a"
        
        rand1<-runif(num.markers, min=0, max=1)
        rand2<-runif(num.markers, min=0, max=1)
        minor1<- map.info$available.Markers[which(map.info$genetic.map[map.info$available.Markers,"MAF"] > rand1)]
        major1 <- map.info$available.Markers[which(!map.info$available.Markers %in% minor1)]
        minor2<- map.info$available.Markers[which(map.info$genetic.map[map.info$available.Markers,"MAF"] > rand2)]
        major2 <- map.info$available.Markers[which(!map.info$available.Markers %in% minor2)]
        
        if(heterozygous.markers){
          all.markers.allele1 <- map.info$available.Markers
          minor1 <- sample(map.info$available.Markers,size = length(map.info$available.Markers)/2,replace=F)
          major1 <- map.info$available.Markers[which(!map.info$available.Markers %in% minor1)]
          minor2<- major1
          major2 <- minor1
        }
        if(inbred.parents){
          allele1[map.info$available.Markers] <- marker.genotypes[par]
          allele2[map.info$available.Markers] <- marker.genotypes[par]
        } else {
          
          allele1[minor1] <- r.marker.genotypes[par]
          allele1[major1] <- d.marker.genotypes[par]
          allele2[minor2] <- r.marker.genotypes[par]
          allele2[major2] <- d.marker.genotypes[par]
        }
        
        parents[,par,1]<-allele1
        parents[,par,2]<-allele2
        
        parents[map.info$rQTL.loci,par,1] <- sample(QTL.alleles,num.rQTL)
        parents[map.info$rQTL.loci,par,2] <- sample(QTL.alleles,num.rQTL)
        ## Write a copy of the SNP alleles at the SNPQTL into an array to use in calculating parental values
        QTLSNP.array[,par,1]<-allele1[map.info$QTLSNP.loci]
        QTLSNP.array[,par,2]<-allele2[map.info$QTLSNP.loci]
      }
      
      
    } else {
      available.delt1 <- map.info$QTLSNP.loci; available.delt2 <- map.info$QTLSNP.loci
      ratio <- (max.delt.allele/2)/map.info$total.SNPQTL.num
      
      for(par in 1:num.parents){
        allele1<-rep(NA, map.info$total.loci.num); allele2<-rep(NA, map.info$total.loci.num)
        num.delt <- sample(0:(map.info$total.SNPQTL.num*ratio),1) # pick a number of deleterious alleles for allele 1 for this parent
        if(num.delt > length(available.delt1)) {break}
        minor1 <- sample(available.delt1,num.delt)  # pick out of the available deleterious alleles for this parent
        available.delt1 <- available.delt1[which(!available.delt1 %in% minor1 )] # this new list makes it so the next parent can't have any of these alleles at allele 1 positions
        major1 <- map.info$QTLSNP.loci[which(!map.info$QTLSNP.loci %in% minor1)]
        
        num.delt <- sample(0:(map.info$total.SNPQTL.num*ratio),1)
        if(num.delt > length(available.delt1[which(!available.delt1 %in% minor1)])) {break}
        minor2 <- sample(available.delt1[which(!available.delt1 %in% minor1)],num.delt)
        available.delt1 <- available.delt1[which(!available.delt1 %in% minor2 )]
        major2 <- map.info$QTLSNP.loci[which(!map.info$QTLSNP.loci %in% minor2)]
        
        allele1[minor1]<-"c"            # minor allele is always "c", major allele is always "a"
        allele1[major1] <-"a"
        allele2[minor2]<-"c"
        allele2[major2]<-"a"
        
        rand1<-runif(num.markers, min=0, max=1)
        rand2<-runif(num.markers, min=0, max=1)
        minor1<-map.info$available.Markers[which(map.info$genetic.map$MAF[map.info$available.Markers] > rand1)]
        major1 <- map.info$available.Markers[which(!map.info$available.Markers %in% minor1)]
        minor2<-map.info$available.Markers[which(map.info$genetic.map$MAF[map.info$available.Markers] > rand2)]
        major2 <- map.info$available.Markers[which(!map.info$available.Markers %in% minor2)]
        
        allele1[minor1] <- r.marker.genotypes
        allele1[major1] <- d.marker.genotypes
        allele2[minor2] <- r.marker.genotypes
        allele2[major2] <- d.marker.genotypes
        
        if(length(which(map.info$genetic.map$types == "fi")) != 0) {
          all.markers.allele1 <- which(map.info$genetic.map$types == "fi")
          minor1 <- sample(all.markers.allele1,size = length(all.markers.allele1)/2,replace=F)
          major1 <- all.markers.allele1[which(!all.markers.allele1 %in% minor1)]
          minor2<- major1
          major2 <- minor1
          allele1[minor1] <- r.fi.marker.genotypes[par]
          allele1[major1] <- d.fi.marker.genotypes[par]
          allele2[minor2] <- r.fi.marker.genotypes[par]
          allele2[major2] <- d.fi.marker.genotypes[par]
        }
        
        allele1[which(is.na(allele1) == T)] <- sample(x = QTL.alleles,size = num.rQTL)
        allele2[which(is.na(allele2) == T)] <- sample(x = QTL.alleles,size = num.rQTL)
        
        parents[,par,1]<-allele1
        parents[,par,2]<-allele2
        
        ## Write a copy of the SNP alleles at the SNPQTL into an array to use in calculating parental values
        QTLSNP.array[,par,1]<-allele1[map.info$QTLSNP.loci]
        QTLSNP.array[,par,2]<-allele2[map.info$QTLSNP.loci]
      }
    }
  }
  
  # Output the snpqtls genotypes to a matrix
  parent.SNP.genos <- matrix(paste(parents[map.info$QTLSNP.loci,,1],parents[map.info$QTLSNP.loci,,2],sep=""),nrow=map.info$total.SNPQTL.num,ncol=num.parents)
  dimnames(parent.SNP.genos) <- list(map.info$genetic.map$loci[map.info$QTLSNP.loci],par.IDs)
  
  # Output the marker gentoypes to a matrix
  parent.markers <- matrix(paste(parents[map.info$available.Markers,,1],parents[map.info$available.Markers,,2],sep=""),nrow=num.markers,ncol=num.parents)
  dimnames(parent.markers) <- list(map.info$genetic.map$loci[map.info$available.Markers],par.IDs)
  
  # Output the rqtl gentoypes to a matrix
  rqtl_mat1 <- apply(parents[map.info$rQTL.loci,,1],2,as.integer)
  rqtl_mat2 <- apply(parents[map.info$rQTL.loci,,2],2,as.integer)
  parent.rqtl <- rqtl_mat1+rqtl_mat2
  dimnames(parent.rqtl) <- list(map.info$genetic.map$loci[map.info$rQTL.loci],par.IDs)
  
  result.parent <- sapply(rep(1:num.parents,1),function(x) length(which(parent.SNP.genos[,x]=="ac" | parent.SNP.genos[,x]=="ca" | parent.SNP.genos[,x]=="cc")))
  result.parent <- unlist(result.parent); names(result.parent) <- colnames(parents)
  # If save is true then then we save the parentSNPgenos object and parent.markers object
  
  parentinfo<-list(genos.3d=parents,delt.allele=result.parent, parent.rQTL.matrix = parent.rqtl, parent.SNPQTL.matrix=parent.SNP.genos, parent.Marker.matrix=parent.markers, parent.IDs=par.IDs,num.parents=num.parents)
  #cat("The returned object is a list containing a 3-D array of parent genotypes\n")
  return(parentinfo)
  # # # # #    End of function    # # # # # # # #
}
