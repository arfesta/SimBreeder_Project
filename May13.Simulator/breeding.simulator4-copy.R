#Create Map####
create.map <- function(num.markers, totalQTL,  numSNPQTL, map.dist = "haldane",
                       numchromos,  maplength, totalloci, chromosomesizerange=.2,
                       signifdigits=2, het.markers=F, seed=F, save=F, rep.num=NULL) {
  
  # Determine sizes of each chromosome##
  avgchromo <- maplength/numchromos
  chromosizes <- floor(runif(numchromos, min=(avgchromo - chromosomesizerange*avgchromo), max=(avgchromo + chromosomesizerange*avgchromo)))
  
  # Determining the number of intervals that will be calculated and the number of loci per chromosome
  numintervals <- totalloci-1 #Specifying the # of intervals that will need to be calcuated
  chromoloci <- rep(totalloci/numchromos,numchromos) #Each chromosome will get an equal number of loci
  #Specify the last loci for each chromosome
  all.loci <- vector("list",12)
  for (i in 1:length(chromoloci)){if (i==1){all.loci[[1]] <- chromoloci[i] } else {all.loci[[i]] <- chromoloci[[i]] + all.loci[[i-1]]}}
  chromo.loci <- unlist(all.loci)
  
  ##Create Map and Recombination Frequencies###
  # Creating empty vectors
  list1<-vector("list", totalloci) # List "chr#" character set for all loci
  list2<-vector("list", totalloci) # List locus# for all chromosomes
  intervals<-rep(NA, totalloci)    # Holds interval distance from each locus to the previous one
  positions<- rep(NA, totalloci)   # Holds the positions of each loci on a chromosome in cM
  locusnames<-rep(NA,totalloci)    # Holds each locus name for all of the loci (chr1_locus#)
  recfreqs<-rep(NA, totalloci)     # Holds Recombination frequency for each of the loci
  
  # Define map function used to create recombination frequencies:
  if(map.dist == "haldane") {mapfun<-function(x){(1-(exp(-2*x)))/2 }} 
  if(map.dist == "kosambi") {mapfun <- function(x){((exp(2*x)-exp(-2*x))/(exp(2*x)+exp(-2*x)))/2 }}
  
  ## Start with chromosome 1 - note that "interval" is distance from locus n to n-1, not n+1
  for(k in 1:numchromos){
    # For each chromosome list 1 will contain the number of loci for that chromosome
    list1[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]<-c(1:chromoloci[k])
    # list2 will hold the chr# that can later be pasted with list1 to make locus names
    list2[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]<-rep(paste("chr",k,sep=""))
    
    # Intervals are initally in cM but divided by 100 to be Morgans.  This is done so that it can be input to haldane function
    # which uses Morgans to calculate Recfreqs
    intervals[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]<-c(round(runif((chromoloci[k]-1),min=round(0.2*chromosizes[k]/chromoloci[k],1),
                                                                              max=round(1.8*chromosizes[k]/chromoloci[k],1)),signifdigits)/100,5)
    
    #Recfreqs for each of the loci are determined by using the interval distance (for each locus) in Morgans and the haldane mapping function:
    #r=(1/2)*(1-e^(2*M))
    recfreqs[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])] <- mapfun(intervals[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])])
    
    #Since intervals are in Morgans, we multiply by 100 to turn it back into cM for positions (cM is unit for chromosome lengths)
    positions[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]<- round(100*cumsum(intervals[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]),digits=8)
  }
  
  ## paste together locus names from list2 and list1
  locusnames<- paste(list2,list1,sep="_")
  
  #Create map that has locusnames, distance between loci, empty vector of types (SNPQTL,rQTL, or Marker), MAFs, position on chromosome, & rec freqs
  map<-data.frame(chr=unlist(list2),loci=locusnames, dist= as.numeric(intervals),types= rep(NA,totalloci), MAF=rep(NA,totalloci), pos=as.numeric(positions),recfreqs=as.numeric(recfreqs))
  
  #Now sample from the map to specify SNPQTL & rQTL, the remainder are potential markers that can be used
  all.loci <- 1:totalloci  # vector that contains 1 through the number of all loci
  
  #Sample 10 markers for each chromosome uniformly distributed
  smp <- num.markers/numchromos
  last.pos.chrs <-  sapply(1:numchromos,function(x){round(map$pos[chromo.loci[x]-1],0)})
  Markers <- vector()
  for(i in 1:length(chromoloci)){
    equal.dist.per.chr <- last.pos.chrs[i]/smp
    map.posit <- map$pos[(chromo.loci[i]-chromo.loci[1] +1):(chromo.loci[i])]
    Marker <- sapply(1:smp,function(x){
      the.pos <- equal.dist.per.chr*x
      which.min(abs(round(map.posit - the.pos,1)))})
    if(i > 1){Marker <- Marker + chromo.loci[i-1]}
    Markers <- c(Markers,Marker)}
  
  map$types[Markers] <- "m"    # Specify in map data frame that all loci which are not qtl or snpqtl are markers
  MAFs <- sample(rbeta(200000,.4,.4),totalloci,replace=F)
  MAFs[which(MAFs > .5)] <- 1- MAFs[which(MAFs > .5)]
  #MAFs <- sample(runif(totalloci,.2,.49),size = totalloci)
  map$MAF <- sample(MAFs,length(map$MAF),replace=F) # Assign minor allele frequencies to markers and qtl
  
  SNPQTL <- sample(all.loci[-Markers],numSNPQTL,replace=FALSE)
  SNPQTLMAFs <- runif(numSNPQTL,min=.01,max=.02)
  map$types[SNPQTL] <- "snpqtl"   # Specify in the map data frame that these loci are snpqtl
  map$MAF[SNPQTL] <- SNPQTLMAFs   # Assign these loci the specificed minor allele frequencies generated by user
  
  
  numQTL <- totalQTL - numSNPQTL  # vector that contains the number of QTL which are not under dominance and have addtive effects
  rQTL <- sample(all.loci[-c(SNPQTL,Markers)],numQTL,replace=F)
  map$types[rQTL] <- "qtl"        # Specify in the map data frame that these loci are qtl
  
  if(het.markers){
    map$MAF[which(map$types == "m")] <- .49}
  
  #Now save map output:
  datevalue <- date()
  datevector <- unlist(strsplit(datevalue,"\\s"))
  timevector <- unlist(strsplit(datevector[4],":"))
  newtime <- paste(timevector[1],timevector[2],sep="h")
  newdate <- paste(datevector[3],datevector[2],datevector[5],sep="")
  namestem <- paste(newdate,newtime,sep="_")
  
  if(save){
    mapname<-paste(rep.num,namestem,"_map.txt")
    write.table(map,file=mapname, quote=F, row.names=F, col.names=T, sep="\t")}
  
  mapinfo<-list(genetic.map=map, total.loci.num=totalloci,total.QTL.num=totalQTL, total.SNPQTL.num=numSNPQTL,
                QTLSNP.loci=sort(SNPQTL),date.time=namestem, rQTL.loci=sort(rQTL), available.Markers=Markers,last.locus.per.chrom=chromo.loci)
  cat("The returned object is a list containing:\n")
  cat("a matrix of 5 columns with locus names, intervals, types (marker, snp, or qtl), MAFs, and recombination fractions;\n")
  cat("a set of scalars with the numbers of total loci, total QTL, SNPs that are also QTL, total SNPs, and map intervals,\n")
  cat("vectors identifying which loci are SNPs, which are SNP-QTL, and which are 'invisible' QTL,\n")
  cat("and the namestem prefix used in saving text output files.\n")
  return(mapinfo)
  # End of function #
}

####Create Founder Parent Population####
create.parents <- function(mapinfo=map, numparents, max.delt.allele, heterozygous.markers = F,
                           QTLsd = 0.5, save = F, rep.num = NULL){
  totalloci <- mapinfo$total.loci.num  # Specifies the total # of loci by pulling from the map object
  locusnames <- mapinfo$genetic.map$loci  # Specifies locus names from map object
  all.MAFs <- mapinfo$genetic.map$MAF     # Assign all.MAFs as the minor allele frequencies in the map object
  totalQTL <- mapinfo$total.QTL.num    # Total number of qtl pulled from map object
  QTLSNPs <- mapinfo$QTLSNP.loci      # A vector of the loci which are snpqtl
  numSNPQTL <- mapinfo$total.SNPQTL.num  # the number of snpqtl pulled from map object
  numQTL <- totalQTL - numSNPQTL  # the number of additive qtl (rQTL)
  QTLoci <- sort(mapinfo$rQTL.loci)    # A vector of the loci which are additive QTL
  all.Markers <- as.numeric(rownames(mapinfo$genetic.map)[-c(QTLSNPs,QTLoci)])
  num.markers <- length(all.Markers)  # the number of potential markers
  
  #  Sample rQTL values from two normal distributions with mean 0 and std deviation half of that specified by QTLsd variable
  #  Use 'ceiling' function for one distribution, 'floor' function for other, then add them together to get centered distribution
  Dist1 <- ceiling(rnorm((numQTL*2), 0, (QTLsd*(sqrt(2)/2)))); Dist2 <- floor(rnorm((numQTL*2), 0, (QTLsd*(sqrt(2)/2))))
  QTLalleles <- Dist1 + Dist2  # Vector contains the values for the rQTL alleles
  
  #  Create Dimension names for 2 arrays
  parIDs<-1:numparents; alleleIDs<-c("a1","a2")
  
  # Parents array holds both alleles for all parents
  # QTLSNP array holds the alleles for the loci of all parents under dominance control
  parents<-array(rep(NA,totalloci*numparents*2), dim=c(totalloci,numparents,2),dimnames=list(locusnames,parIDs,alleleIDs))
  QTLSNParray<-array(0, dim=c(numSNPQTL,numparents,2),dimnames=list(QTLSNPs,parIDs,alleleIDs))
  
  #  Create empty vectors to hold either allele1 or allele2  for each parent
  allele1<-rep(NA, totalloci); allele2<-rep(NA, totalloci)
  
  # Assign deleterious alleles
  available.delt1 <- QTLSNPs; available.delt2 <- QTLSNPs
  D.marker.genotypes <- vector(); R.marker.genotypes <- vector()
  for (i in 1:26){
    D.marker.genotypes <- c(D.marker.genotypes,paste(LETTERS[i],LETTERS, sep=""))
    R.marker.genotypes <-  c(R.marker.genotypes,paste(letters[i],letters, sep=""))}
  
  ratio <- (max.delt.allele/2)/numSNPQTL
  # For each parent assign them alleles
  while(anyNA(parents)){
    available.delt1 <- QTLSNPs
    for(par in 1:numparents){
      num.delt <- sample(0:(numSNPQTL*ratio),1) # pick a number of deleterious alleles for allele 1 for this parent
      if(num.delt > length(available.delt1)) {break}
      minor1 <- sample(available.delt1,num.delt)  # pick out of the avaiable deleterious alleles for this parent
      available.delt1 <- available.delt1[which(!available.delt1 %in% minor1 )] # this new list makes it so the next parent can't have any of these alleles at allele 1 positions
      major1 <- QTLSNPs[which(!QTLSNPs %in% minor1)]
      
      num.delt <- sample(0:(numSNPQTL*ratio),1)
      if(num.delt > length(available.delt1[which(!available.delt1 %in% minor1)])) {break}
      minor2 <- sample(available.delt1[which(!available.delt1 %in% minor1)],num.delt)
      available.delt1 <- available.delt1[which(!available.delt1 %in% minor2 )]
      
      
      major2 <- QTLSNPs[which(!QTLSNPs %in% minor2)]
      allele1[minor1]<-"c"            # minor allele is always "c", major allele is always "a"
      allele1[major1] <-"a"
      allele2[minor2]<-"c"
      allele2[major2]<-"a"
      
      rand1<-runif(num.markers, min=0, max=1)
      rand2<-runif(num.markers, min=0, max=1)
      minor1<-all.Markers[which(all.MAFs[all.Markers] > rand1)]
      major1 <- all.Markers[which(!all.Markers %in% minor1)]
      minor2<-all.Markers[which(all.MAFs[all.Markers] > rand2)]
      major2 <- all.Markers[which(!all.Markers %in% minor2)]
      if(heterozygous.markers){
        all.Markers.allele1 <- all.Markers
        minor1 <- sample(all.Markers,size = length(all.Markers)/2,replace=F)
        major1 <- all.Markers[which(!all.Markers %in% minor1)]
        minor2<- major1
        major2 <- minor1
      }
      
      allele1[minor1] <- R.marker.genotypes[par]
      allele1[major1] <- D.marker.genotypes[par]
      allele2[minor2] <- R.marker.genotypes[par]
      allele2[major2] <- D.marker.genotypes[par]
      
      parents[,par,1]<-allele1
      parents[,par,2]<-allele2
      
      parents[QTLoci,par,1] <- sample(QTLalleles,numQTL)
      parents[QTLoci,par,2] <- sample(QTLalleles,numQTL)
      ## Write a copy of the SNP alleles at the SNPQTL into an array to use in calculating parental values
      QTLSNParray[,par,1]<-allele1[QTLSNPs]
      QTLSNParray[,par,2]<-allele2[QTLSNPs]
    }}
  
  # Output the snpqtls genotypes to a matrix
  parentSNPgenos <- matrix(paste(parents[QTLSNPs,,1],parents[QTLSNPs,,2],sep=""),nrow=numSNPQTL,ncol=numparents)
  dimnames(parentSNPgenos) <- list(locusnames[QTLSNPs],parIDs)
  
  # Output the marker gentoypes to a matrix
  parent.markers <- matrix(paste(parents[all.Markers,,1],parents[all.Markers,,2],sep=""),nrow=num.markers,ncol=numparents)
  dimnames(parent.markers) <- list(locusnames[all.Markers],parIDs)
  
  result.parent <- sapply(rep(1:numparents,1),function(x) length(which(parentSNPgenos[,x]=="ac" | parentSNPgenos[,x]=="ca" | parentSNPgenos[,x]=="cc")))
  result.parent <- unlist(result.parent)
  # If save is true then then we save the parentSNPgenos object and parent.markers object
  if(save){
    parfilename=paste(rep.num,"parent-SNPQTL-genotypes.txt")
    write.table(parentSNPgenos,parfilename, quote = F, row.names = T, col.names = T, sep= "\t")
    
    parfilename=paste(rep.num,"parent-marker-genotypes.txt")
    write.table(parent.markers,parfilename, quote = F, row.names = T, col.names = T, sep= "\t")
  }
  
  ###IMP that genos.3d is first list element and numparents is second becuase that is what createTGV routine expects###
  parentinfo<-list(genos.3d=parents,delt.allele=result.parent,num.parents=numparents, parent.SNPQTL.matrix=parentSNPgenos, parent.Marker.matrix=parent.markers, parent.IDs=parIDs)
  cat("The returned object is a list containing a 3-D array of parent genotypes\n")
  cat("a matrix of parent SNP genotypes with individuals in columns, a scalar of the number of parents,\n")
  cat("matrices of the parental genotypes at all SNPs, parental QTL allele values, genotypes and values at SNP QTL,\n")
  cat("two matrices with variances of QTL alleles & genetic values and SNP QTL alleles & genetic values, and a list \n")
  cat("of the numerical values assigned to major and minor alleles at SNP-QTL in the parents.\n")
  return(parentinfo)
  # # # # #    End of function    # # # # # # # #
}

####Create  TGV Genetic Value####
create.parents.TGV <- function(parents, mapinfo,
                               A, a, dom.coeff,
                               usesnpeffect = F,snpeffects = NULL,
                               save = F, rep.num = NULL, prefix = NULL){
  
  locusnames <- mapinfo$genetic.map$loci # The locus names pulled from the mab object
  QTLSNPs <- mapinfo$QTLSNP.loci     # vector of the loci which are snpqtl
  QTLSNPnum <- parents$genos.3d[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
  markers <- mapinfo$available.Markers# a list of all the markers pulled from map object
  num.markers <- length(markers)
  markers.select <- sort(sample(markers,num.markers,replace=F),decreasing=F)
  num.markers <- length(markers.select) # length of markers that were selected
  marker.select.genos <- parents$genos.3d[markers.select,,] # genotypes of the markers pulled from the current generation
  map.markers <- mapinfo$genetic.map[markers.select,c(1,6)]
  numQTL <- length(mapinfo$rQTL.loci) # the number of additive qtl
  parIDs <- parents$parent.IDs
  
  
  
  # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
  numSNPQTL <- mapinfo$total.SNPQTL.num # the number of loci which are snpqtl
  numparents <- length(parIDs) # the number of parents
  QTLSNPvalues <-matrix(NA,nrow=numSNPQTL,ncol=numparents) # matrix to hold snpqtl values
  marker.values <- matrix(NA,nrow=num.markers,ncol=numparents) # matrix to hold marker values
  Capital.genotypes <- vector()
  Lowercase.genotypes <- vector()
  for (i in 1:26){
    Capital.genotypes <- c(Capital.genotypes,paste(LETTERS[i],LETTERS, sep=""))
    Lowercase.genotypes <-  c(Lowercase.genotypes,paste(letters[i],letters, sep=""))}
  
  for (i in 1:length(parIDs)){
    QTLSNPaa <- which(QTLSNPnum[,i,1]=="a" & QTLSNPnum[,i,2]=="a")
    QTLSNPcc <- which(QTLSNPnum[,i,1]=="c" & QTLSNPnum[,i,2]=="c")
    QTLSNPac <- which(QTLSNPnum[,i,1]=="a" & QTLSNPnum[,i,2]=="c")
    QTLSNPca <- which(QTLSNPnum[,i,1]=="c"  & QTLSNPnum[,i,2]=="a")
    if (dom.coeff==0){
      QTLSNPvalues[QTLSNPaa,i] <- A*2
      QTLSNPvalues[QTLSNPcc,i] <- a*2
      QTLSNPvalues[QTLSNPac,i] <- (A+a+dom.coeff)
      QTLSNPvalues[QTLSNPca,i] <- (A+a+dom.coeff) } else {
        QTLSNPvalues[QTLSNPaa,i] <- A
        if(usesnpeffect){
          QTLSNPvalues[QTLSNPcc,i] <- snpeffects[a]
        } else{
          QTLSNPvalues[QTLSNPcc,i] <- a}
        QTLSNPvalues[QTLSNPac,i] <- (A*dom.coeff)
        QTLSNPvalues[QTLSNPca,i] <- (A*dom.coeff)
      }
    
    markers.aa <- which(marker.select.genos[,i,1] %in% Capital.genotypes & marker.select.genos[,i,2] %in% Capital.genotypes)
    markers.cc <- which(marker.select.genos[,i,1] %in% Lowercase.genotypes & marker.select.genos[,i,2] %in% Lowercase.genotypes)
    markers.ac <-which(marker.select.genos[,i,1] %in% Capital.genotypes & marker.select.genos[,i,2] %in% Lowercase.genotypes)
    markers.ca <- which(marker.select.genos[,i,1] %in% Lowercase.genotypes & marker.select.genos[,i,2] %in% Capital.genotypes)
    
    marker.values[markers.aa,i] <- "0"
    marker.values[markers.cc,i] <- "2"
    marker.values[markers.ac,i] <- "1"
    marker.values[markers.ca,i] <- "1"
  }
  marker.values <- t(marker.values)
  colnames(marker.values) <- markers.select
  rownames(marker.values) <- parIDs
  
  # Convert the 'invisible' rQTL genotypes to numeric matrices, merge the alleles to paired values also
  parQTLallele1 <- matrix(as.integer(parents$genos.3d[mapinfo$rQTL,,1]),nrow=numQTL,ncol=numparents)
  colnames(parQTLallele1) <- c(parIDs)
  parQTLallele2 <- matrix(as.integer(parents$genos.3d[mapinfo$rQTL,,2]),nrow=numQTL,ncol=numparents)
  colnames(parQTLallele2) <- c(parIDs)
  QTLvalues     <- matrix(paste(parQTLallele1,parQTLallele2,sep=","),nrow=numQTL,ncol=numparents)
  dimnames(QTLvalues)<-list(locusnames[mapinfo$rQTL],parIDs)
  
  # Genetic values of progeny
  geneticvals <- colSums(QTLSNPvalues) + colSums(parQTLallele1) + colSums(parQTLallele2)
  
  if(save) {
    pedfilename=paste(rep.num,prefix,".txt", sep="")
    write.table(geneticvals,pedfilename, quote = F, row.names = F, col.names = T, sep="\t")}
  
  TGV <- list(genetic.values=geneticvals, SNP.value.matrix=QTLSNPvalues, markers.matrix=marker.values, marker.loci=markers.select, marker.map=map.markers)
  return(TGV)
}

####Create Unscaled Phenos####
create.parents.phenos <- function(parents.TGV, h2, E.var = NULL, save = F, folder = NULL) {
  geneticvals <-parents.TGV$genetic.values
  totalindiv <-length(geneticvals)
  if(length(E.var) > 0){
    phenos <- geneticvals + rnorm(totalindiv,mean = 0,sd = sqrt(E.var))
  } else {
    phenos <- geneticvals + rnorm(totalindiv,mean = 0,sd = sqrt(var(geneticvals)/h2))}
  true.h2 <- round(var(geneticvals)/var(phenos),2)
  pheno.info<-list(phenos=phenos, genetic.values=geneticvals, true.h2 = true.h2)
  return(pheno.info)
}

#Extract OP genotypic values####
select.OP.Phenotype.Values <- function(mapinfo,parentinfo,parent.phenos,parents.TGV,num.select,crossprog,dom.coeff,A,a,
                                       h2, E.var = NULL, run.parallel=F, n.cores= NULL ) {
  library(abind)
  numpars <- as.numeric(length(colnames(parentinfo$genos.3d)))                         
  Parent1 <-  rep(1:numpars,each=(numpars-1)); Parent2 <- vector()
  for(i in 1:numpars){
    all.parents <- 1:numpars
    Parent2 <- c(Parent2,all.parents[-i])}
  Progeny <- rep(crossprog,length(Parent2))
  OP.crosses <- data.frame(par1=Parent1,par2=Parent2,numprog=Progeny)
  
  if(run.parallel){
    library(parallel)
    g.val <- mclapply(1:nrow(OP.crosses),mc.cores = n.cores,FUN = function(x,crosses=OP.crosses){
      
      y= crosses[x,]
      chromoloci <- mapinfo$last.locus.per.chrom   # The last number of each loci for all chromsomes
      numchromos <- length(chromoloci)
      QTLSNPs      <- mapinfo$QTLSNP.loci # Vector of loci which are snpqtl
      m <- mapinfo$genetic.map
      
      par1<-match(y[1],colnames(parentinfo$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
      par2<-match(y[2],colnames(parentinfo$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
      crossprog<-as.numeric(y[3]) #assigns number of progeny to be the third column for cross "X"
      
      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes1<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
      gametes2<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
      
      par1.alleles <- (parentinfo$genos.3d[,par1,])
      par2.alleles <- (parentinfo$genos.3d[,par2,])
      chr.ind.r <- vector("list")
      for(each in 1:crossprog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromoloci)){
          last.pos <- chromoloci[i]
          t <- sample(seq(0,1,.0001),chromoloci[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:crossprog) {
        allele <- sample(1:2,1)
        end <- chromoloci[numchromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes1[z:recombination.spots[each],i] <- par1.alleles[z:recombination.spots[each],allele]
          } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}
      
      chr.ind.r <- vector("list")
      for(each in 1:crossprog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromoloci)){
          last.pos <- chromoloci[i]
          t <- sample(seq(0,1,.0001),chromoloci[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:crossprog) {
        allele <- sample(1:2,1)
        end <- chromoloci[numchromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes2[z:recombination.spots[each],i] <- par2.alleles[z:recombination.spots[each],allele]
          } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}
      array.out <- abind(gametes1,gametes2,along=3)
      
      locusnames <- mapinfo$genetic.map$loci # The locus names pulled from the mab object
      QTLSNPs <- mapinfo$QTLSNP.loci     # vector of the loci which are snpqtl
      QTLSNPnum <- array.out[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
      numQTL <- length(mapinfo$rQTL.loci) # the number of additive qtl
      
      # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
      numSNPQTL <- mapinfo$total.SNPQTL.num # the number of loci which are snpqtl
      QTLSNPvalues <-matrix(NA,nrow=numSNPQTL,ncol=1) # matrix to hold snpqtl values
      
      QTLSNPaa <- which(QTLSNPnum[,1]=="a" & QTLSNPnum[,2]=="a")
      QTLSNPcc <- which(QTLSNPnum[,1]=="c" & QTLSNPnum[,2]=="c")
      QTLSNPac <- which(QTLSNPnum[,1]=="a" & QTLSNPnum[,2]=="c")
      QTLSNPca <- which(QTLSNPnum[,1]=="c"  & QTLSNPnum[,2]=="a")
      if (dom.coeff==0){
        QTLSNPvalues[QTLSNPaa,1] <- A*2
        QTLSNPvalues[QTLSNPcc,1] <- a*2
        QTLSNPvalues[QTLSNPac,1] <- (A+a+dom.coeff)
        QTLSNPvalues[QTLSNPca,1] <- (A+a+dom.coeff) } else {
          QTLSNPvalues[QTLSNPaa,1] <- A
          QTLSNPvalues[QTLSNPcc,1] <- a
          QTLSNPvalues[QTLSNPac,1] <- (A*dom.coeff)
          QTLSNPvalues[QTLSNPca,1] <- (A*dom.coeff)}
      
      parQTLallele1 <- as.numeric(array.out[mapinfo$rQTL.loci,,1])
      parQTLallele2 <- as.numeric(array.out[mapinfo$rQTL.loci,,2])
      
      # Genetic values of progeny
      geneticvalue <- sum(QTLSNPvalues[,1]) + sum(parQTLallele1) + sum(parQTLallele2)
      geneticvalue
    })
  } else {
    g.val <- lapply(1:nrow(OP.crosses),FUN = function(x,crosses=OP.crosses){
      numchromos <- length(chromoloci)
      y= crosses[x,]
      chromoloci <- mapinfo$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- mapinfo$QTLSNP.loci # Vector of loci which are snpqtl
      m <- mapinfo$genetic.map
      
      par1<-match(y[1],colnames(parentinfo$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
      par2<-match(y[2],colnames(parentinfo$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
      crossprog<-as.numeric(y[3]) #assigns number of progeny to be the third column for cross "X"
      
      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes1<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
      gametes2<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
      
      par1.alleles <- (parentinfo$genos.3d[,par1,])
      par2.alleles <- (parentinfo$genos.3d[,par2,])
      chr.ind.r <- vector("list")
      for(each in 1:crossprog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromoloci)){
          last.pos <- chromoloci[i]
          t <- sample(seq(0,1,.0001),chromoloci[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:crossprog) {
        allele <- sample(1:2,1)
        end <- chromoloci[numchromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes1[z:recombination.spots[each],i] <- par1.alleles[z:recombination.spots[each],allele]
          } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}
      
      chr.ind.r <- vector("list")
      for(each in 1:crossprog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromoloci)){
          last.pos <- chromoloci[i]
          t <- sample(seq(0,1,.0001),chromoloci[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:crossprog) {
        allele <- sample(1:2,1)
        end <- chromoloci[numchromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes2[z:recombination.spots[each],i] <- par2.alleles[z:recombination.spots[each],allele]
          } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}
      array.out <- abind(gametes1,gametes2,along=3)
      
      locusnames <- mapinfo$genetic.map$loci # The locus names pulled from the mab object
      QTLSNPs <- mapinfo$QTLSNP.loci     # vector of the loci which are snpqtl
      QTLSNPnum <- array.out[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
      numQTL <- length(mapinfo$rQTL.loci) # the number of additive qtl
      
      # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
      numSNPQTL <- mapinfo$total.SNPQTL.num # the number of loci which are snpqtl
      QTLSNPvalues <-matrix(NA,nrow=numSNPQTL,ncol=1) # matrix to hold snpqtl values
      
      QTLSNPaa <- which(QTLSNPnum[,1]=="a" & QTLSNPnum[,2]=="a")
      QTLSNPcc <- which(QTLSNPnum[,1]=="c" & QTLSNPnum[,2]=="c")
      QTLSNPac <- which(QTLSNPnum[,1]=="a" & QTLSNPnum[,2]=="c")
      QTLSNPca <- which(QTLSNPnum[,1]=="c"  & QTLSNPnum[,2]=="a")
      if (dom.coeff==0){
        QTLSNPvalues[QTLSNPaa,1] <- A*2
        QTLSNPvalues[QTLSNPcc,1] <- a*2
        QTLSNPvalues[QTLSNPac,1] <- (A+a+dom.coeff)
        QTLSNPvalues[QTLSNPca,1] <- (A+a+dom.coeff) } else {
          QTLSNPvalues[QTLSNPaa,1] <- A
          QTLSNPvalues[QTLSNPcc,1] <- a
          QTLSNPvalues[QTLSNPac,1] <- (A*dom.coeff)
          QTLSNPvalues[QTLSNPca,1] <- (A*dom.coeff)}
      
      parQTLallele1 <- as.numeric(array.out[mapinfo$rQTL.loci,,1])
      parQTLallele2 <- as.numeric(array.out[mapinfo$rQTL.loci,,2])
      
      # Genetic values of progeny
      geneticvalue <- sum(QTLSNPvalues[,1]) + sum(parQTLallele1) + sum(parQTLallele2)
      geneticvalue
    })
  }
  
  g.val <- unlist(g.val)
  totalindiv <- length(g.val)
  if(length(E.var) > 0){ phenos <- g.val + rnorm(totalindiv,mean = 0,sd = E.var) } else {
    phenos <- g.val + rnorm(totalindiv,mean = 0,sd = sqrt(var(g.val)/h2))}
  
  mean.phenos <- vector()
  first <- 1
  last <- numpars-1
  for(i in 1:numpars){
    mean.phenos <- c(mean.phenos,mean(phenos[first:last]))
    first <- last+1
    last <- first + numpars-2}
  
  names(mean.phenos) <- colnames(parentinfo$genos.3d)
  pars <- sort(mean.phenos,decreasing = T)[1:num.select]
  delt.alleles <- parentinfo$delt.allele[as.numeric(names(pars))]
  genetic.values <- parents.phenos$genetic.values[as.numeric(names(pars))]
  phenos <- parents.phenos$phenos[as.numeric(names(pars))]
  genos.3d <- parentinfo$genos.3d[,as.numeric(names(pars)),]
  marker.matrix <- parents.TGV$markers.matrix[as.numeric(names(pars)),]
  
  out <- list(pars=pars,delt.alleles=delt.alleles,genetic.values=genetic.values,phenos=phenos, genos.3d=genos.3d,
              marker.matrix=marker.matrix)
  return(out)
}

####Create Cross Design & Pedigree's####
create.cross.design <- function(mapinfo=map, parentinfo, prog.percross, gen,
                                mating.design="AssortiveMating", coancest.thresh=T, 
                                use.par.marker.effects = F, use.op.par.phenos = F,
                                crossfile= NULL, save=F, parentphenos = NULL){
  library(reshape2)
  coancest=NULL
  NumParents=ncol((parentinfo$genos.3d[,,1]))
  parent.names <- colnames(parentinfo$genos.3d)
  
  if(mating.design == "Self"){
    par1 <- parent.names
    par2 <- parent.names
    prog.par <- rep(prog.percross,length(par1))
    crossdesign <- cbind(par1,par2,prog.par)}
  
  #Must be divisible by 4
  if(mating.design == "AssortiveMating"){
    total <- NumParents/4
    if(gen==1){
      if(use.par.marker.effects | use.op.par.phenos) {sorted.ebvs <- parentinfo$pars
      } else {
        if(length(parentinfo$phenos) > 0) {sorted.ebvs <- parentinfo$phenos
        } else {sorted.ebvs <- sort(parentphenos$phenos,decreasing=T)}}
      
      relationshipmatrix <- as.data.frame(diag(1,nrow=NumParents,ncol=NumParents))
      colnames(relationshipmatrix) <- names(sorted.ebvs); rownames(relationshipmatrix) <- names(sorted.ebvs)
      
      top.25 <- names(sorted.ebvs)[1:total]
      top.50 <- names(sorted.ebvs)[(total+1):(total*2)]
      bottom.50 <- names(sorted.ebvs)[(total*2+1):(total*3)]
      bottom.25 <- names(sorted.ebvs)[(total*3+1):(total*4)]
    } else {
      sorted.ebvs <- sort(parentinfo$select.EBVs,decreasing = T)
      relationshipmatrix <- parentinfo$relmat[which(colnames(parentinfo$relmat) %in% names(sorted.ebvs)),which(colnames(parentinfo$relmat) %in% names(sorted.ebvs))]
      top.25 <- names(sorted.ebvs)[1:total]
      top.50 <- names(sorted.ebvs)[(total+1):(total*2)]
      bottom.50 <- names(sorted.ebvs)[(total*2+1):(total*3)]
      bottom.25 <- names(sorted.ebvs)[(total*3+1):(total*4)]
    }
    
    all <- vector("list")
    for(i in 1:total){all[[i]] <- top.25[-i]}
    
    if(coancest.thresh){
      E <- upper.tri(relationshipmatrix); relationshipmatrix[E] <- NA
      coancestry <- seq(.01,1,.005)
      crossdesign <- matrix(NA,nrow=1,ncol=2)
      
      for(i in 1:length(coancestry)){
        new <- subset(melt(as.matrix(t(relationshipmatrix))), value < coancestry[i])
        matches1 <- match(new[,1],names(sorted.ebvs))
        matches2 <- match(new[,2],names(sorted.ebvs))
        a <- sorted.ebvs[matches1]
        b <- sorted.ebvs[matches2]
        new$mean.bv <- (a+b)/2
        new <- new[order(new$mean.bv,decreasing=T),]
        coancest <- coancestry[i]
        crossdesign <- matrix(NA,nrow=1,ncol=2)
        for(each in 1:length(top.25)){
          par <- top.25[each]
          options1 <- which(new[,1] %in% par)
          options2 <- which(new[,2] %in% par )
          options1.other <- options1[which(as.character(new[options1,2]) %in% all[[each]])]
          options2.other <- options2[which(new[options2,1] %in% all[[each]])]
          options <- sort(c(options1.other,options2.other))
          this.cross <- new[options[1:2],c(1:2)] ; colnames(this.cross) <- c("V1","V2")
          #if(anyNA(this.cross)){break}
          crossdesign <- rbind(crossdesign,this.cross)
          if(each == 1){crossdesign <- crossdesign[-1,]}
          new <- new[-c(which(new[,1] %in% this.cross[,1] & new[,2] %in% this.cross[,2] )),]
        }
        
        top50.samp <-top.50
        for(each in 1:length(top.25)){
          par <- top.25[each]
          options1 <-   which(new[,1] %in% par)
          options2 <-  which( new[,2] %in% par )
          options1.other <- options1[which(as.character(new[options1,2]) %in% top50.samp)]
          options2.other <- options2[which(new[options2,1] %in% top50.samp)]
          options <- sort(c(options1.other,options2.other))
          this.cross <- new[options[1],c(1:2)] ; colnames(this.cross) <- c("V1","V2")
          #if(anyNA(this.cross)){break}
          crossdesign <- rbind(crossdesign,this.cross)
          new <- new[-c(which(new[,1] %in% this.cross[,1] & new[,2] %in% this.cross[,2] )),]
          top50.samp <- top50.samp[-c(which(top50.samp %in% this.cross[,1] | top50.samp %in% this.cross[,2]))]
        }
        
        parent1.2.top25 <- sample(top.25,length(top.25),replace=F)
        parent2.2.bottom50 <- sample(bottom.50,length(top.25),replace=F)
        parent1.3.top50 <- sample(top.50,length(top.25),replace=F)
        parent2.3.bottom25 <- sample(bottom.25,length(top.25),replace=F)
        parent1.4.top50 <- sample(top.50,length(top.25),replace=F)
        parent2.4.bottom50 <- sample(bottom.50,length(top.25),replace=F)
        
        all.other.crosses1 <- c(parent1.2.top25,parent1.3.top50,parent1.4.top50)
        all.other.crosses2 <- c(parent2.2.bottom50,parent2.3.bottom25,parent2.4.bottom50)
        all.other.crosses <- cbind(all.other.crosses1,all.other.crosses2)
        
        colnames(all.other.crosses) <- c("V1","V2")
        crossdesign <- rbind(crossdesign,all.other.crosses)
        if(!anyNA(crossdesign)){break}}
      crossdesign[,3] <- rep(prog.percross,nrow(crossdesign))
    } else {
      parent2 <- vector()
      for(each in 1:length(top.25)){
        par1 <- each
        par2 <- sample(all[[each]],2,replace=F)
        parent2 <- c(parent2,par2)
        for(i in 1:length(par2)){
          all[[par2[i]]] <- all[[par2[i]]][-each]
        }}
      parent1.top25 <- rep(top.25,each=2)
      parent2.top25 <- parent2
      parent1.1.top25 <- sample(top.25,length(top.25),replace=F)
      parent2.1.top50 <- sample(top.50,length(top.25),replace=F)
      parent1.2.top25 <- sample(top.25,length(top.25),replace=F)
      parent2.2.bottom50 <- sample(bottom.50,length(top.25),replace=F)
      parent1.3.top50 <- sample(top.50,length(top.25),replace=F)
      parent2.3.bottom25 <- sample(bottom.25,length(top.25),replace=F)
      parent1.4.top50 <- sample(top.50,length(top.25),replace=F)
      parent2.4.bottom50 <- sample(bottom.50,length(top.25),replace=F)
      
      pars1 <- c(parent1.top25,parent1.1.top25,parent1.2.top25,parent1.3.top50,parent1.4.top50)
      
      pars2 <- c(parent2.top25,parent2.1.top50,parent2.2.bottom50,parent2.3.bottom25,parent2.4.bottom50)
      
      
      crossdesign <- cbind(par1=pars1,
                           par2=pars2,
                           prog=rep(prog.percross,length(pars1)))}}
  
  if(mating.design == "RandomMating"){
    crossdesign= matrix(sample(colnames(parentinfo$genos.3d),replace=FALSE, size=NumParents),nrow = NumParents/2,ncol=2)
    crossdesign <- cbind(crossdesign,as.numeric(rep(prog.percross,length(crossdesign[,1]))))}
  
  if(mating.design == "SP" ){
    if(coancest.thresh){
      par1 <- vector(); par2 <- vector()
      relationshipmatrix <- parentinfo$relmat[which(colnames(parentinfo$relmat) %in% parentinfo$select.ped.ids),which(colnames(parentinfo$relmat) %in% parentinfo$select.ped.ids)]
      E <- upper.tri(relationshipmatrix); relationshipmatrix[E] <- NA
      coancestry <- seq(.01,1,.005)
      for(i in 1:length(coancestry)){
        new <- subset(melt(as.matrix(t(relationshipmatrix))), value < coancestry[i])
        matches1 <- match(new[,1],names(parentinfo$select.EBVs))
        matches2 <- match(new[,2],names(parentinfo$select.EBVs))
        a <- parentinfo$select.EBVs[matches1]
        b <- parentinfo$select.EBVs[matches2]
        new$mean.bv <- (a+b)/2
        new <- new[order(new$mean.bv,decreasing=T),]
        coancest <- coancestry[i]
        crossdesign <- matrix(NA,nrow=(NumParents/2),ncol=3)
        for(i in 1:(NumParents/2)) {
          crossdesign[i,1] <-  new[1,1]
          crossdesign[i,2] <-  new[1,2]
          crossdesign[i,3] <- prog.percross
          new <- new[which(new[,1] != crossdesign[i,1] & new[,1] != crossdesign[i,2]),]
          new <- new[which(new[,2] != crossdesign[i,1] & new[,2] != crossdesign[i,2]),]}
        if (!anyNA(crossdesign)) {break}
      }} else {
        if(gen==1){
          if(use.par.marker.effects | use.op.par.phenos){sorted.ebvs <- parentinfo$pars
          }} else {
            sorted.ebvs <- sort(parentinfo$select.EBVs,decreasing = T)}
        n.s.ebvs <- length(sorted.ebvs)
        parent1 <- names(sorted.ebvs)[c(seq(from = 1,to = n.s.ebvs,by = 2))]
        parent2 <- names(sorted.ebvs)[c(seq(2,n.s.ebvs,2))]
        crossdesign <- cbind(par1=parent1,par2=parent2,prog=rep(prog.percross,length(parent1)))}}
  
  if(mating.design == "MP") {
    if(coancest.thresh){
      par1 <- vector()
      par2 <- vector()
      relationshipmatrix <- parentinfo$relmat[which(colnames(parentinfo$relmat) %in% parentinfo$select.ped.ids),which(colnames(parentinfo$relmat) %in% parentinfo$select.ped.ids)]
      D <- relationshipmatrix
      coancestry <- seq(.01,1,.005)
      for(i in 1:length(coancestry)){
        D[D < coancestry[i]] <- 0
        D[D > coancestry[i]] <- 1
        D <- 1-D
        minimum <- colSums(D, na.rm=T)
        coancest <- coancestry[i]
        D <- relationshipmatrix
        #if(min(minimum) >= 90) break }
        sorted.minimum <- sort(minimum, decreasing=F)
        relationshipmatrix <- relationshipmatrix[,names(sorted.minimum)]
        E <- upper.tri(relationshipmatrix)
        relationshipmatrix[E] <- NA
        new <- subset(melt(as.matrix(t(relationshipmatrix))), value < coancest)
        matches1 <- match(new[,1],names(parentinfo$select.EBVs))
        matches2 <- match(new[,2],names(parentinfo$select.EBVs))
        a <- parentinfo$select.EBVs[matches1]
        b <- parentinfo$select.EBVs[matches2]
        new$mean.bv <- (a+b)/2
        sib1 <- sample(colnames(relationshipmatrix)[seq(1,NumParents,2)],NumParents/2)
        sib2 <- colnames(relationshipmatrix)[seq(2,NumParents,2)]
        test <- new[order(-new[,4]),]
        par1 <- test[1:(length(parentinfo$selections)/2),1]
        par2 <- test[1:(length(parentinfo$selections)/2),2]
        progeny <- rep(prog.percross, length(par1))
        crossdesign=matrix(data=c(par1,par2,progeny),nrow=length(progeny),ncol=3)
        if (!anyNA(crossdesign)) {break}}} else {
          if(gen==1){
            if(use.par.marker.effects | use.op.par.phenos){sorted.ebvs <- parentinfo$pars
            }} else {
              sorted.ebvs <- sort(parentinfo$select.EBVs,decreasing = T)}
          n.s.ebvs <- as.numeric(length(sorted.ebvs))
          parent1 <- rep(names(sorted.ebvs[1]),n.s.ebvs-1)
          parent2 <- names(sorted.ebvs[-1])
          parent3 <- rep(names(sorted.ebvs[2]),(n.s.ebvs/2)+1)
          parent4 <- names(sorted.ebvs[3:((n.s.ebvs/2)+3)])
          crossdesign <- cbind(par1=c(parent1,parent3),par2=c(parent2,parent4),prog=rep(prog.percross,length(c(parent1,parent3))))
          
        }}
  
  if(mating.design == "cross.file.input") {
    colClasses = c("integer", "integer", "integer")
    crossdesign = as.matrix(read.table(crossfile,header=F,sep="\t", colClasses))}
  
  numcrosses <- nrow(crossdesign)
  totalprogeny <- sum(as.numeric(crossdesign[,3]))
  par1list<-vector(); par2list<-vector()
  par1id<-vector(); par2id<-vector()
  generation<- rep(gen,totalprogeny)
  
  if (gen==1){
    totalindiv   <- totalprogeny+length(unique(c(crossdesign[,1],crossdesign[,2])))
    cumul.total <- length(unique(c(crossdesign[,1],crossdesign[,2]))) + totalprogeny
    indivIDs<-c(seq(max(as.numeric(parent.names))+1,max(as.numeric(parent.names))+(totalprogeny)))
    for (m in 1:numcrosses){
      crossprog<-as.numeric(crossdesign[m,3])
      par1<-as.character(crossdesign[m,1])
      par2<-as.character(crossdesign[m,2])
      par1list<-c(par1list,rep(par1,crossprog))
      par2list<-c(par2list,rep(par2,crossprog))}
    
    uni <- length(unique(c(crossdesign[,1],crossdesign[,2])))
    parentped <- cbind("ID"=rep(unique(c(crossdesign[,1],crossdesign[,2])),1),"Par1" = rep(0,uni), "Par2"= rep(0,uni),"gener"=rep(0,uni))
    progped<-cbind("ID"=indivIDs[1:length(indivIDs)], "Par1"=par1list[1:length(par1list)],"Par2"=par2list[1:length(par1list)],"gener"=generation[1:length(generation)])
    full.ped <- rbind(parentped,progped)
    selection.pedigree <- full.ped
    
  } else {
    indivIDs<- c("indivIDs",seq(max(as.numeric(parentinfo$fullped[,1]))+1,max(as.numeric(parentinfo$fullped[,1]))+totalprogeny))
    cumul.total<- parentinfo$cumulative.total + totalprogeny
    vector2 <- parentinfo$select.ped.ids
    for (m in 1:numcrosses){
      crossprog<-as.numeric(crossdesign[m,3])
      par1<-c(crossdesign[m,1])
      par2<-c(crossdesign[m,2])
      #par1id<- (vector2[par1])
      #par2id<- (vector2[par2])
      par1list<-c(par1list,rep(as.vector(par1),crossprog))
      par2list<-c(par2list,rep(as.vector(par2),crossprog))}
    progped<-cbind("ID"=indivIDs[2:length(indivIDs)], "Par1"=par1list,"Par2"=par2list,"gener"=generation)
    full.ped <- rbind(parentinfo$fullped,progped)
    selection.pedigree <- rbind(parentinfo$ped,progped)
  }
  
  #if(save) {
  #  pedfilename=paste(rep.num,mapinfo$date.time,prefix,"-pedigree.txt",sep="")
  #  write.table(progped,pedfilename, quote = F, row.names = F, col.names = T, sep=" ")}
  
  ###Imp that genos.3d is first in list item, and total progeny is next list item bc in createtgv routine it expects genos.3d as first element and total prog as second
  cross.design <-list(crossdesign=crossdesign, total.progeny.number=totalprogeny, progeny.pedigree=progped, full.pedigree= full.ped,
                      numparents=parentinfo$num.parents, coancestry.threshold=coancest, selection.ped=selection.pedigree,
                      numcrosses=numcrosses,crossprog=crossprog,cumul.total=cumul.total, parent.IDs=progped[,1])
  return(cross.design)}

####Create Make Crosses####
make.crosses <- function(parentinfo,mapinfo=map,crossdesign=cross.design, run.parallel = F, num.cores = NULL){
  crossdesign <- crossdesign$crossdesign
  num.crosses <- as.numeric(nrow(crossdesign))
  num.chromos <- length(mapinfo$last.locus.per.chrom)
  library(parallel); library(abind)
  if(run.parallel){
    gametes1 <- mclapply(1:num.crosses,function(x){
      chromoloci <- mapinfo$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- mapinfo$QTLSNP.loci # Vector of loci which are snpqtl
      m <- mapinfo$genetic.map
      
      par1<-match(crossdesign[x,1],colnames(parentinfo$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
      par2<-match(crossdesign[x,2],colnames(parentinfo$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
      crossprog<-as.numeric(crossdesign[x,3]) #assigns number of progeny to be the third column for cross "X"
      
      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes1<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
      
      par1.alleles <- (parentinfo$genos.3d[,par1,])
      par2.alleles <- (parentinfo$genos.3d[,par2,])
      chr.ind.r <- vector("list")
      for(each in 1:crossprog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromoloci)){
          last.pos <- chromoloci[i]
          t <- sample(seq(0,1,.0001),chromoloci[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1}
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:crossprog) {
        allele <- sample(1:2,1)
        end <- chromoloci[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes1[z:recombination.spots[each],i] <- par1.alleles[z:recombination.spots[each],allele]
          } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1}}
      
      gametes1
    },mc.cores=num.cores)
    gametes1 <- matrix(unlist(gametes1), ncol = num.crosses*mean(as.numeric(crossdesign[,3])))
    gametes2 <- mclapply(1:num.crosses,function(x){
      chromoloci <- mapinfo$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- mapinfo$QTLSNP.loci # Vector of loci which are snpqtl
      m <- mapinfo$genetic.map
      
      par1<-match(crossdesign[x,1],colnames(parentinfo$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
      par2<-match(crossdesign[x,2],colnames(parentinfo$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
      crossprog<-as.numeric(crossdesign[x,3]) #assigns number of progeny to be the third column for cross "X"
      
      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes2<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
      
      par1.alleles <- (parentinfo$genos.3d[,par1,])
      par2.alleles <- (parentinfo$genos.3d[,par2,])
      
      chr.ind.r <- vector("list")
      for(each in 1:crossprog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromoloci)){
          last.pos <- chromoloci[i]
          t <- sample(seq(0,1,.0001),chromoloci[1],replace = F)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:crossprog) {
        allele <- sample(1:2,1)
        end <- chromoloci[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes2[z:recombination.spots[each],i] <- par2.alleles[z:recombination.spots[each],allele]
          } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each] + 1
        }}
      
      gametes2
    },mc.cores=num.cores)
    gametes2 <- matrix(unlist(gametes2), ncol = num.crosses*mean(as.numeric(crossdesign[,3])))
    genos.3d <- abind(gametes1,gametes2,along = 3)
    out <- list(genos.3d=genos.3d)
  } else {
    gametes1 <- lapply(1:num.crosses,function(x){
      chromoloci <- mapinfo$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- mapinfo$QTLSNP.loci # Vector of loci which are snpqtl
      m <- mapinfo$genetic.map
      
      par1<-match(crossdesign[x,1],colnames(parentinfo$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
      par2<-match(crossdesign[x,2],colnames(parentinfo$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
      crossprog<-as.numeric(crossdesign[x,3]) #assigns number of progeny to be the third column for cross "X"
      
      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes1<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
      
      par1.alleles <- (parentinfo$genos.3d[,par1,])
      par2.alleles <- (parentinfo$genos.3d[,par2,])
      chr.ind.r <- vector("list")
      for(each in 1:crossprog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromoloci)){
          last.pos <- chromoloci[i]
          t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          if(length(l) > 2){ l <- l[1:2]} else{l <- l}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:crossprog) {
        allele <- sample(1:2,1)
        end <- chromoloci[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes1[z:recombination.spots[each]-1,i] <- par1.alleles[z:recombination.spots[each]-1,allele]
          } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each]
        }}
      
      gametes1
    })
    gametes1 <- matrix(unlist(gametes1), ncol = num.crosses*mean(as.numeric(crossdesign[,3])))
    gametes2 <- lapply(1:num.crosses,function(x){
      chromoloci <- mapinfo$last.locus.per.chrom   # The last number of each loci for all chromsomes
      QTLSNPs      <- mapinfo$QTLSNP.loci # Vector of loci which are snpqtl
      m <- mapinfo$genetic.map
      
      par1<-match(crossdesign[x,1],colnames(parentinfo$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
      par2<-match(crossdesign[x,2],colnames(parentinfo$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
      crossprog<-as.numeric(crossdesign[x,3]) #assigns number of progeny to be the third column for cross "X"
      
      # Create empty matrix to hold gametes
      # dimensions are (total # of loci) x  (# of cross progeny)
      # rownames are the loci names
      gametes2<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
      
      par1.alleles <- (parentinfo$genos.3d[,par1,])
      par2.alleles <- (parentinfo$genos.3d[,par2,])
      
      chr.ind.r <- vector("list")
      for(each in 1:crossprog){
        first.pos <- 1
        ch.r <- vector("list")
        for(i in 1:length(chromoloci)){
          last.pos <- chromoloci[i]
          t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])
          while(length(l) < 1){t <- sample(seq(0,1,.0001),chromoloci[1],replace = T)
          l <- which(t < m[first.pos:last.pos,7])}
          if(length(l) > 2){ l <- l[1:2]} else{l <- l}
          ch.r[[i]] <- seq(first.pos,last.pos,1)[l]
          first.pos <- last.pos+ 1
        }
        chr.ind.r[[each]] <- unlist(ch.r)}
      
      for(i in 1:crossprog) {
        allele <- sample(1:2,1)
        end <- chromoloci[num.chromos]
        z <-1
        recombination.spots <- chr.ind.r[[i]]
        for (each in 1:length(recombination.spots)) {
          if (each < length(recombination.spots)){
            gametes2[z:recombination.spots[each]-1,i] <- par2.alleles[z:recombination.spots[each]-1,allele]
          } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
          if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
          z <- recombination.spots[each]
        }}
      
      gametes2
    })
    gametes2 <- matrix(unlist(gametes2), ncol = num.crosses*mean(as.numeric(crossdesign[,3])))
    genos.3d <- abind(gametes1,gametes2,along = 3)
    out <- list(genos.3d=genos.3d)
  }
  return(out)}

####Create TGV2####
create.progeny.TGV <- function(proginfo, mapinfo=map, crossdesign=cross.design, A=1,a=-100, dom.coeff=Dominance.Coeff,
                               save=F, rep.num=NULL, prefix=NULL){
  
  QTLSNPs <- mapinfo$QTLSNP.loci     # vector of the loci which are snpqtl
  QTLSNPnum <- proginfo$genos.3d[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
  markers <- mapinfo$available.Markers# a list of all the markers pulled from map object
  num.markers <- length(markers) # length of markers that were selected
  markers.select <- sort(sample(markers,num.markers,replace=F),decreasing=F)
  marker.select.genos <- proginfo$genos.3d[markers.select,,] # genotypes of the markers pulled from the current generation
  map.markers <- mapinfo$genetic.map[markers.select,c(1,6)]
  numQTL <- length(mapinfo$rQTL.loci) # the number of additive qtl
  
  parIDs <- crossdesign$parent.IDs
  length.prog <- length(parIDs)
  
  numSNPQTL <- mapinfo$total.SNPQTL.num # the number of loci which are snpqtl
  numparents <- length(parIDs) # the number of parents
  #QTLSNPvalues <-matrix(NA,nrow=numSNPQTL,ncol=numparents) # matrix to hold snpqtl values
  marker.values <- matrix(NA,nrow=num.markers,ncol=numparents) # matrix to hold marker values
  Capital.genotypes <- vector()
  Lowercase.genotypes <- vector()
  for (i in 1:26){
    Capital.genotypes <- c(Capital.genotypes,paste(LETTERS[i],LETTERS, sep=""))
    Lowercase.genotypes <-  c(Lowercase.genotypes,paste(letters[i],letters, sep=""))}
  
  QTLSNPaa <- sapply(1:length.prog,function(x){
    A*length(which(QTLSNPnum[,x,1]=="a" & QTLSNPnum[,x,2]=="a"))},simplify = T)
  QTLSNPcc <- sapply(1:length.prog,function(x){
    a*length(which(QTLSNPnum[,x,1]=="c" & QTLSNPnum[,x,2]=="c"))},simplify = T)
  QTLSNPac <- sapply(1:length.prog,function(x){
    (A * dom.coeff) * length(which(QTLSNPnum[,x,1]=="a" & QTLSNPnum[,x,2]=="c"))},simplify = T)
  QTLSNPca <- sapply(1:length.prog,function(x){
    (A * dom.coeff) * length(which(QTLSNPnum[,x,1]=="c" & QTLSNPnum[,x,2]=="a"))},simplify = T)
  QTLSNPvalues <- QTLSNPaa+QTLSNPcc+QTLSNPac+QTLSNPca
  
  markers.aa <-   sapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% Capital.genotypes & marker.select.genos[,x,2] %in% Capital.genotypes)},simplify = T)
  markers.cc <-   sapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% Lowercase.genotypes & marker.select.genos[,x,2] %in% Lowercase.genotypes)},simplify = T)
  markers.ac <-   sapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% Capital.genotypes & marker.select.genos[,x,2] %in% Lowercase.genotypes)},simplify = T)
  markers.ca <-   sapply(1:length.prog,function(x){which(marker.select.genos[,x,1] %in% Lowercase.genotypes & marker.select.genos[,x,2] %in% Capital.genotypes)},simplify = T)
  marker.values <- matrix(NA,nrow=num.markers,ncol=length.prog) # matrix to hold marker values
  
  for(i in 1:length.prog){
    marker.values[markers.aa[[i]],i] <- "0"
    marker.values[markers.cc[[i]],i] <- "2"
    marker.values[markers.ac[[i]],i] <- "1"
    marker.values[markers.ca[[i]],i] <- "1"}
  
  marker.values <- t(marker.values)
  colnames(marker.values) <- markers.select
  rownames(marker.values) <- parIDs
  #removed marker.vals from return (markers.matrix=marker.values)
  geneticvals <- QTLSNPvalues + colSums( matrix(as.integer(proginfo$genos.3d[mapinfo$rQTL,,1]),nrow=numQTL,ncol=length.prog)) + colSums(matrix(as.integer(proginfo$genos.3d[mapinfo$rQTL,,2]),nrow=numQTL,ncol=length.prog))
  names(geneticvals) <- parIDs
  
  TGV <- list(genetic.values=geneticvals, SNP.value.matrix=QTLSNPvalues,markers.matrix=marker.values, marker.loci=markers.select, marker.map=map.markers)
  return(TGV)}

####Create Unscaled Phenos####
create.progeny.phenos <- function(progeny.TGV, h2, crossdesign=cross.design, E.var = NULL, save = F,
                                  folder = NULL) {
  geneticvals <-progeny.TGV$genetic.values
  totalindiv <-length(geneticvals)
  if(length(E.var) > 0){
    phenos <- geneticvals + rnorm(totalindiv,mean = 0,sd = sqrt(E.var))
  } else {
    phenos <- geneticvals + rnorm(totalindiv,mean = 0,sd = sqrt(var(geneticvals)/h2))}
  true.h2 <- round(var(geneticvals)/var(phenos),2)
  
  num.progeny <- as.numeric(crossdesign$crossdesign[1,3])
  crosses <- nrow(crossdesign$crossdesign)
  fam.mean.g <- vector(); fam.mean.p <- vector()
  first <- 1
  last <- num.progeny
  for(i in 1:crosses){
    fam.mean.g <- c(fam.mean.g,mean(geneticvals[first:last]))
    fam.mean.p <- c(fam.mean.p,mean(phenos[first:last]))
    first <- last +1
    last <- last + num.progeny}
  fam.mean.h2 <-  var(fam.mean.g)/var(fam.mean.p)
  
  pheno.info<-list(phenos=phenos, genetic.values=geneticvals, true.h2 = true.h2,fam.mean.h2 = fam.mean.h2)
  return(pheno.info)
}

####Extract Selections####
ex <- function(num.cores=NULL, mapinfo=map, crossdesign=cross.design, past.tgv, past.phenos, parentinfo,
               progenyinfo, progenyTGV, progenyphenos, selection.strategy, 
               among.family.selection = NULL, within.family.selection = NULL, numSelections.within.family = 1,
               numSelections.among.family = NULL, relmatrix = "pedigree",
               prefix = NULL, rep.num = NULL, reduced = F, weighted=F) {
  
  library(MatrixModels); library(parallel); library(pedigreemm); library(pedigree); library(bigmemory)
  
  # Extracting crossdesign, number of parents, complete(full) & progeny pedigree from progenyinfo object
  cross.design <- crossdesign$crossdesign   #Get cross design file for current generation
  prog.percross <- as.numeric(cross.design[1,3]) #Determine the number of progeny per cross
  prog.pedigree <- crossdesign$progeny.pedigree  #Extract the full progeny pedigree for this generation
  gen <- as.numeric(prog.pedigree[1,4])  #The current generation number is the last column in the prog pedigree
  numparents <- crossdesign$numparents  #Get number of parents from parent generation
  full.ped <- crossdesign$full.pedigree 
  selection.ped <- crossdesign$selection.ped
  
  # Extracting parent markers, progeny tgv & markers, as well as creating marker map from progenyTGV object
  if (gen==1){
    if(reduced){
      parent.markers <- parentinfo$genos.3d[mapinfo$available.Markers,(names(parentinfo$pars)),]
    } else { parent.markers <- past.tgv$markers.matrix}} else {parent.markers <- parentinfo$all.markers}
  prog.markers <- progenyTGV$markers.matrix
  progeny.phenos <- progenyphenos$phenos
  prog.genetic.values <- progenyTGV$genetic.values
  map.markers <- progenyTGV$marker.map
  colnames(map.markers) <- c("chr","pos")
  map.markers$chr <- as.character(map.markers$chr)
  
  # Generate pedigree using synbreed if 1st gen otherwise just use full.ped object
  # Create 2 objects which hold all markers and all phenotypes respectively
  if (gen==1) {
    ped <- pedigree(label = full.ped[,1],sire = full.ped[,2],dam = full.ped[,3])
    #all.markers <- rbind(parent.markers,)
    if(reduced){  all.phenos <- c(parentinfo$phenos,progeny.phenos)
    all.genetic.vals <- c(parentinfo$genetic.values,prog.genetic.values)
    } else { all.phenos <- c(past.phenos$phenos,progeny.phenos)
    all.genetic.vals <- c(past.phenos$genetic.values,prog.genetic.values)}
  } else {
    ped <- pedigree(label = selection.ped[,1],sire = selection.ped[,2],dam = selection.ped[,3])
    #all.markers <- rbind(parent.markers,prog.markers)
    all.phenos <- c(parentinfo$all.phenos,progeny.phenos)
    all.genetic.vals <- c(past.phenos$genetic.values,prog.genetic.values)
  }
  
  
  first.in.family <- 1
  Selections <- vector()
  Capital.genotypes <- vector(); Lowercase.genotypes <- vector()
  for (i in 1:26){
    Capital.genotypes <- c(Capital.genotypes,paste(LETTERS[i],LETTERS,LETTERS[i],LETTERS, sep=""))
    Lowercase.genotypes <-  c(Lowercase.genotypes,paste(letters[i],letters, letters[i], letters, sep=""))}
  
  if (among.family.selection == "ABLUP"){
    AI <- as.big.matrix(as.matrix(getA(ped)))
    P=all.phenos[match(colnames(AI),names(all.phenos))]
    s <-  ncol(AI)/8
    e <- s-1
    first.seq <- seq(1,ncol(AI),ncol(AI)/8)
    {
      my.env <- new.env()
      
      my.env$l1 <-  AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e)),
                       j=c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e)),drop=T]
      my.env$l2 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e)),
                      j=c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e)),drop=T]
      my.env$l3 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e)),
                      j=c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e)),drop=T]
      my.env$l4 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e)),
                      j=c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e)),drop=T]
      my.env$l5 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e)),
                      j=c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
      my.env$l6 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[7]:(first.seq[7]+e)),
                      j=c(first.seq[1]:(first.seq[1]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
      my.env$l7 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[8]:(first.seq[8]+e)),
                      j=c(first.seq[1]:(first.seq[1]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
      
      my.env$l8 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e)),
                      j=c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e)),drop=T]
      my.env$l9 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e)),
                      j=c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e)),drop=T]
      my.env$l10 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e)),
                       j=c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e)),drop=T]
      my.env$l11 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e)),
                       j=c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
      my.env$l12 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[7]:(first.seq[7]+e)),
                       j=c(first.seq[2]:(first.seq[2]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
      my.env$l13 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[8]:(first.seq[8]+e)),
                       j=c(first.seq[2]:(first.seq[2]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
      
      my.env$l14 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e)),
                       j=c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e)),drop=T]
      
      my.env$l15 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e)),
                       j=c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e)),drop=T]
      my.env$l16 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e)),
                       j=c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
      my.env$l17 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[7]:(first.seq[7]+e)),
                       j=c(first.seq[3]:(first.seq[3]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
      my.env$l18 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[8]:(first.seq[8]+e)),
                       j=c(first.seq[3]:(first.seq[3]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
      
      
      my.env$l19 <- AI[i=c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e)),
                       j=c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e)),drop=T]
      my.env$l20 <- AI[i=c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e)),
                       j=c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
      my.env$l21<- AI[i=c(first.seq[4]:(first.seq[4]+e),first.seq[7]:(first.seq[7]+e)),
                      j=c(first.seq[4]:(first.seq[4]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
      my.env$l22<- AI[i=c(first.seq[4]:(first.seq[4]+e),first.seq[8]:(first.seq[8]+e)),
                      j=c(first.seq[4]:(first.seq[4]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
      
      
      my.env$l23 <- AI[i=c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e)),
                       j=c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
      my.env$l24 <- AI[i=c(first.seq[5]:(first.seq[5]+e),first.seq[7]:(first.seq[7]+e)),
                       j=c(first.seq[5]:(first.seq[5]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
      my.env$l25 <- AI[i=c(first.seq[5]:(first.seq[5]+e),first.seq[8]:(first.seq[8]+e)),
                       j=c(first.seq[5]:(first.seq[5]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
      
      my.env$l26 <- AI[i=c(first.seq[6]:(first.seq[6]+e),first.seq[7]:(first.seq[7]+e)),
                       j=c(first.seq[6]:(first.seq[6]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
      my.env$l27 <- AI[i=c(first.seq[6]:(first.seq[6]+e),first.seq[8]:(first.seq[8]+e)),
                       j=c(first.seq[6]:(first.seq[6]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
      
      my.env$l28 <- AI[i=c(first.seq[7]:(first.seq[7]+e),first.seq[8]:(first.seq[8]+e)),
                       j=c(first.seq[7]:(first.seq[7]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
    }
    
    rm(AI); gc()
    
    {
      s.phenos1 <- P[c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e))]
      s.phenos2 <- P[c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e))]
      s.phenos3 <- P[c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos4 <- P[c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos5 <- P[c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos6 <- P[c(first.seq[1]:(first.seq[1]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos7 <- P[c(first.seq[1]:(first.seq[1]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos8 <- P[c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e))]
      s.phenos9 <- P[c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos10 <- P[c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos11 <- P[c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos12 <- P[c(first.seq[2]:(first.seq[2]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos13 <- P[c(first.seq[2]:(first.seq[2]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos14 <- P[c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos15 <- P[c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos16 <- P[c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos17 <- P[c(first.seq[3]:(first.seq[3]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos18 <- P[c(first.seq[3]:(first.seq[3]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos19 <- P[c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos20 <- P[c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos21 <- P[c(first.seq[4]:(first.seq[4]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos22 <- P[c(first.seq[4]:(first.seq[4]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos23 <- P[c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos24 <- P[c(first.seq[5]:(first.seq[5]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos25 <- P[c(first.seq[5]:(first.seq[5]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos26 <- P[c(first.seq[6]:(first.seq[6]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos27 <- P[c(first.seq[6]:(first.seq[6]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos28 <- P[c(first.seq[7]:(first.seq[7]+e),first.seq[8]:(first.seq[8]+e))]
      
      pheno.list <- list(s.phenos1,s.phenos2,s.phenos3,s.phenos4,s.phenos5,s.phenos6,s.phenos7,s.phenos8,s.phenos9,
                         s.phenos10,s.phenos11,s.phenos12,s.phenos13,s.phenos14,s.phenos15,s.phenos16,s.phenos17,
                         s.phenos18,s.phenos19,s.phenos20,s.phenos21,s.phenos22,s.phenos23,s.phenos24,s.phenos25,s.phenos26,s.phenos27,s.phenos28)
    }
    
    applyMyFun <- function(idx, env) {
      eval(parse(text = paste0("result <- env$", ls(env)[idx])))
      the.data <- as.matrix(result)
      the.data <- chol2inv(chol(the.data))
      n.col <- ncol(the.data)
      h.2 <- var(prog.genetic.values)/var(progeny.phenos)
      lambda <- (1-h.2)/h.2 
      I <- diag(n.col)
      s.phenos <- pheno.list[[idx]]
      sol <- solve(rBind(cBind(n.col, t(rep(1,n.col))),
                         cBind(rep(1,n.col), (I + (lambda*the.data)))),
                   matrix(c(sum(s.phenos),
                            c(as.vector(s.phenos)))))
      sol <- sol[-1,1]
      sol
    }
    
    index <- 1:28
    names(index) <- 1:28
    
    the.blups <- mclapply(index,applyMyFun,env = my.env, mc.cores=28)
    rm(my.env);gc()
    
    {
      SG1 <- apply(rbind(the.blups[[1]][1:(length(s.phenos1)/2)],
                         the.blups[[2]][1:(length(s.phenos1)/2)],
                         the.blups[[3]][1:(length(s.phenos1)/2)],
                         the.blups[[4]][1:(length(s.phenos1)/2)],
                         the.blups[[5]][1:(length(s.phenos1)/2)],
                         the.blups[[6]][1:(length(s.phenos1)/2)],
                         the.blups[[7]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG2 <- apply(rbind(the.blups[[1]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[8]][1:(length(s.phenos1)/2)],
                         the.blups[[9]][1:(length(s.phenos1)/2)],
                         the.blups[[10]][1:(length(s.phenos1)/2)],
                         the.blups[[11]][1:(length(s.phenos1)/2)],
                         the.blups[[12]][1:(length(s.phenos1)/2)],
                         the.blups[[13]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG3 <- apply(rbind(the.blups[[2]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[8]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[14]][1:(length(s.phenos1)/2)],
                         the.blups[[15]][1:(length(s.phenos1)/2)],
                         the.blups[[16]][1:(length(s.phenos1)/2)],
                         the.blups[[17]][1:(length(s.phenos1)/2)],
                         the.blups[[18]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG4 <- apply(rbind(the.blups[[3]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[9]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[14]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[19]][1:(length(s.phenos1)/2)],
                         the.blups[[20]][1:(length(s.phenos1)/2)],
                         the.blups[[21]][1:(length(s.phenos1)/2)],
                         the.blups[[22]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG5 <- apply(rbind(the.blups[[4]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[10]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[15]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[19]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[23]][1:(length(s.phenos1)/2)],
                         the.blups[[24]][1:(length(s.phenos1)/2)],
                         the.blups[[25]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG6 <- apply(rbind(the.blups[[5]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[11]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[16]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[20]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[23]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[26]][1:(length(s.phenos1)/2)],
                         the.blups[[27]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG7 <- apply(rbind(the.blups[[6]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[12]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[17]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[21]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[24]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[26]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[28]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG8 <- apply(rbind(the.blups[[7]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[13]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[18]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[22]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[25]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[27]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[28]][((length(s.phenos1)/2)+1):length(s.phenos1)]),2,mean)
    }
    
    progeny.blups <- c(SG1,SG2,SG3,SG4,SG5,SG6,SG7,SG8); names(progeny.blups) <- names(P)
    l <- which(names(progeny.blups) %in% names(progeny.phenos))
    progeny.blups <- progeny.blups[l]
    
  }
  
  
  first <- 1
  last <- prog.percross
  mean.progeny.blups <- vector()
  for(each in 1:nrow(cross.design)){
    mean.progeny.blups <- c(mean.progeny.blups,mean(progeny.blups[first:last]))
    first <- last+1
    last <- first+prog.percross-1
  }
  sorted.mean.progeny.blups <- sort(mean.progeny.blups,decreasing=T)
  top.280.families <- match(sorted.mean.progeny.blups, mean.progeny.blups)
  
  if(selection.strategy == "within.family"){
    if(within.family.selection == "GBLUP"){
      require(parallel)
      m <- mapinfo$available.Markers
      if(gen==1){parentmarkers <- parentinfo$genos.3d[m,,]} else {parentmarkers <- parentinfo$all.markers}
      progmarkers <- progenyinfo$genos.3d[m,,]
      prog.1 <- progmarkers[,,1]; colnames(prog.1) <- names(progeny.phenos)
      prog.2 <- progmarkers[,,2]; colnames(prog.2) <- names(progeny.phenos)
      allele1 <- cbind(parentmarkers[,,1],prog.1)
      allele2 <- cbind(parentmarkers[,,2],prog.2)
      A <- as.big.matrix(as.matrix(getA(ped)))
      #this.order <- match(colnames(allele1),rownames(A))
      #A <- A[c(this.order),c(this.order)]
      #rm(progmarkers,parentmarkers)
      cor.mat <- matrix(nrow=ncol(allele1),ncol=ncol(allele1))
      for(i in 1:120) {
        these <- unique(allele1[i,])
        for(all in 1:length(these)) {
          posit <- which(allele1[i,] == these[all])
          for(each in 1:length(posit)){ 
            fill <- posit[each]
            others <- posit[-c(1:each)]
            cor.mat[fill,others] <- cor.mat[fill,others] + 1/240
          }}
        print(i)}
      
      cor.mat2 <- matrix(0,nrow=ncol(allele1),ncol=ncol(allele1))
      for(i in 1:120) {
        these <- unique(allele2[i,])
        for(all in 1:length(these)) {
          posit <- which(allele2[i,] == these[all])
          for(each in 1:length(posit)){ 
            fill <- posit[each]
            others <- posit[-c(1:each)]
            cor.mat2[fill,others] <- cor.mat2[fill,others] + 1/240
          }}}
      
      g.rel.mat <- cor.mat + cor.mat2
      rm(cor.mat,cor.mat2); gc()
      diag(g.rel.mat) <- 1
      diag(g.rel.mat) <- 1
      g.rel.mat[lower.tri(g.rel.mat)] <- g.rel.mat[upper.tri(g.rel.mat)]
      rm(allele1,allele2,prog.1,prog.2); gc()
      
      AI <- .97*g.rel.mat + .03*A
      rm(g.rel.mat,A); gc()
      
      s <-  ncol(AI)/8
      e <- s-1
      first.seq <- seq(1,ncol(AI),ncol(AI)/8)
      {
        my.env <- new.env()
        
        my.env$l1 <-  AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e)),
                         j=c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e)),drop=T]
        my.env$l2 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e)),
                        j=c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e)),drop=T]
        my.env$l3 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e)),
                        j=c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e)),drop=T]
        my.env$l4 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e)),
                        j=c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e)),drop=T]
        my.env$l5 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e)),
                        j=c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
        my.env$l6 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[7]:(first.seq[7]+e)),
                        j=c(first.seq[1]:(first.seq[1]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
        my.env$l7 <- AI[i=c(first.seq[1]:(first.seq[1]+e),first.seq[8]:(first.seq[8]+e)),
                        j=c(first.seq[1]:(first.seq[1]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
        
        my.env$l8 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e)),
                        j=c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e)),drop=T]
        my.env$l9 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e)),
                        j=c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e)),drop=T]
        my.env$l10 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e)),
                         j=c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e)),drop=T]
        my.env$l11 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e)),
                         j=c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
        my.env$l12 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[7]:(first.seq[7]+e)),
                         j=c(first.seq[2]:(first.seq[2]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
        my.env$l13 <- AI[i=c(first.seq[2]:(first.seq[2]+e),first.seq[8]:(first.seq[8]+e)),
                         j=c(first.seq[2]:(first.seq[2]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
        
        my.env$l14 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e)),
                         j=c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e)),drop=T]
        
        my.env$l15 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e)),
                         j=c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e)),drop=T]
        my.env$l16 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e)),
                         j=c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
        my.env$l17 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[7]:(first.seq[7]+e)),
                         j=c(first.seq[3]:(first.seq[3]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
        my.env$l18 <- AI[i=c(first.seq[3]:(first.seq[3]+e),first.seq[8]:(first.seq[8]+e)),
                         j=c(first.seq[3]:(first.seq[3]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
        
        
        my.env$l19 <- AI[i=c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e)),
                         j=c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e)),drop=T]
        my.env$l20 <- AI[i=c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e)),
                         j=c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
        my.env$l21<- AI[i=c(first.seq[4]:(first.seq[4]+e),first.seq[7]:(first.seq[7]+e)),
                        j=c(first.seq[4]:(first.seq[4]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
        my.env$l22<- AI[i=c(first.seq[4]:(first.seq[4]+e),first.seq[8]:(first.seq[8]+e)),
                        j=c(first.seq[4]:(first.seq[4]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
        
        
        my.env$l23 <- AI[i=c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e)),
                         j=c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e)),drop=T]
        my.env$l24 <- AI[i=c(first.seq[5]:(first.seq[5]+e),first.seq[7]:(first.seq[7]+e)),
                         j=c(first.seq[5]:(first.seq[5]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
        my.env$l25 <- AI[i=c(first.seq[5]:(first.seq[5]+e),first.seq[8]:(first.seq[8]+e)),
                         j=c(first.seq[5]:(first.seq[5]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
        
        my.env$l26 <- AI[i=c(first.seq[6]:(first.seq[6]+e),first.seq[7]:(first.seq[7]+e)),
                         j=c(first.seq[6]:(first.seq[6]+e),first.seq[7]:(first.seq[7]+e)),drop=T]
        my.env$l27 <- AI[i=c(first.seq[6]:(first.seq[6]+e),first.seq[8]:(first.seq[8]+e)),
                         j=c(first.seq[6]:(first.seq[6]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
        
        my.env$l28 <- AI[i=c(first.seq[7]:(first.seq[7]+e),first.seq[8]:(first.seq[8]+e)),
                         j=c(first.seq[7]:(first.seq[7]+e),first.seq[8]:(first.seq[8]+e)),drop=T]
      }
      rm(AI); gc()
      s.phenos1 <- all.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e))]
      s.phenos2 <- all.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e))]
      s.phenos3 <- all.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos4 <- all.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos5 <- all.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos6 <- all.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos7 <- all.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos8 <- all.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e))]
      s.phenos9 <- all.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos10 <- all.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos11 <- all.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos12 <- all.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos13 <- all.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos14 <- all.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos15 <- all.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos16 <- all.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos17 <- all.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos18 <- all.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos19 <- all.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos20 <- all.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos21 <- all.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos22 <- all.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos23 <- all.phenos[c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos24 <- all.phenos[c(first.seq[5]:(first.seq[5]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos25 <- all.phenos[c(first.seq[5]:(first.seq[5]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos26 <- all.phenos[c(first.seq[6]:(first.seq[6]+e),first.seq[7]:(first.seq[7]+e))]
      s.phenos27 <- all.phenos[c(first.seq[6]:(first.seq[6]+e),first.seq[8]:(first.seq[8]+e))]
      s.phenos28 <- all.phenos[c(first.seq[7]:(first.seq[7]+e),first.seq[8]:(first.seq[8]+e))]
      
      pheno.list <- list(s.phenos1,s.phenos2,s.phenos3,s.phenos4,s.phenos5,s.phenos6,s.phenos7,s.phenos8,s.phenos9,
                         s.phenos10,s.phenos11,s.phenos12,s.phenos13,s.phenos14,s.phenos15,s.phenos16,s.phenos17,
                         s.phenos18,s.phenos19,s.phenos20,s.phenos21,s.phenos22,s.phenos23,s.phenos24,s.phenos25,s.phenos26,s.phenos27,s.phenos28)
      applyMyFun <- function(idx, env) {
        eval(parse(text = paste0("result <- env$", ls(env)[idx])))
        the.data <- as.matrix(result)
        the.data <- chol2inv(chol(the.data))
        n.col <- ncol(the.data)
        h.2 <- var(prog.genetic.values)/var(progeny.phenos)
        lambda <- (1-h.2)/h.2 
        I <- diag(n.col)
        s.phenos <- pheno.list[[idx]]
        sol <- solve(rBind(cBind(n.col, t(rep(1,n.col))),
                           cBind(rep(1,n.col), (I + (lambda*the.data)))),
                     matrix(c(sum(s.phenos),
                              c(as.vector(s.phenos)))))
        sol <- sol[-1,1]
        sol
      }
      
      index <- 1:28
      names(index) <- 1:28
      
      the.blups <- mclapply(index,applyMyFun,env = my.env, mc.cores=28)
      rm(my.env);gc()
      
      SG1 <- apply(rbind(the.blups[[1]][1:(length(s.phenos1)/2)],
                         the.blups[[2]][1:(length(s.phenos1)/2)],
                         the.blups[[3]][1:(length(s.phenos1)/2)],
                         the.blups[[4]][1:(length(s.phenos1)/2)],
                         the.blups[[5]][1:(length(s.phenos1)/2)],
                         the.blups[[6]][1:(length(s.phenos1)/2)],
                         the.blups[[7]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG2 <- apply(rbind(the.blups[[1]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[8]][1:(length(s.phenos1)/2)],
                         the.blups[[9]][1:(length(s.phenos1)/2)],
                         the.blups[[10]][1:(length(s.phenos1)/2)],
                         the.blups[[11]][1:(length(s.phenos1)/2)],
                         the.blups[[12]][1:(length(s.phenos1)/2)],
                         the.blups[[13]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG3 <- apply(rbind(the.blups[[2]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[8]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[14]][1:(length(s.phenos1)/2)],
                         the.blups[[15]][1:(length(s.phenos1)/2)],
                         the.blups[[16]][1:(length(s.phenos1)/2)],
                         the.blups[[17]][1:(length(s.phenos1)/2)],
                         the.blups[[18]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG4 <- apply(rbind(the.blups[[3]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[9]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[14]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[19]][1:(length(s.phenos1)/2)],
                         the.blups[[20]][1:(length(s.phenos1)/2)],
                         the.blups[[21]][1:(length(s.phenos1)/2)],
                         the.blups[[22]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG5 <- apply(rbind(the.blups[[4]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[10]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[15]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[19]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[23]][1:(length(s.phenos1)/2)],
                         the.blups[[24]][1:(length(s.phenos1)/2)],
                         the.blups[[25]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG6 <- apply(rbind(the.blups[[5]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[11]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[16]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[20]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[23]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[26]][1:(length(s.phenos1)/2)],
                         the.blups[[27]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG7 <- apply(rbind(the.blups[[6]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[12]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[17]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[21]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[24]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[26]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[28]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG8 <- apply(rbind(the.blups[[7]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[13]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[18]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[22]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[25]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[27]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[28]][((length(s.phenos1)/2)+1):length(s.phenos1)]),2,mean)
      
      gprogeny.blups <- c(SG1,SG2,SG3,SG4,SG5,SG6,SG7,SG8); names(gprogeny.blups) <- names(all.phenos)
      l <- which(names(gprogeny.blups) %in% names(progeny.phenos))
      gprogeny.blups <- gprogeny.blups[l]
      plot(gprogeny.blups,prog.genetic.values)
      Selections <- vector()
      first.in.family <- 1
      for(family in 1:length(cross.design[,3])){
        num.offspring <- as.numeric(cross.design[family,3])
        last.in.family <- num.offspring + first.in.family - 1
        temp <- (gprogeny.blups[first.in.family:last.in.family]) ; sorted <-sort(temp,decreasing=TRUE)
        BestOne <- which(temp==sorted[1:numSelections.within.family]) ; selected <- BestOne + first.in.family - 1 ;Selections <- c(Selections, selected)
        first.in.family <- last.in.family + 1 }
      if(length(numSelections.among.family) > 0) {
        Selections <- Selections[top.280.families[1:numSelections.among.family]]} else {
          Selections <- Selections[top.280.families]
        }
    } else {
      Selections <- vector()
      first.in.family <- 1
      for(family in 1:length(cross.design[,3])){
        num.offspring <- as.numeric(cross.design[family,3])
        last.in.family <- num.offspring + first.in.family - 1
        temp <- (progeny.blups[first.in.family:last.in.family]) ; sorted <-sort(temp,decreasing=TRUE)
        BestOne <- which(temp==sorted[1:numSelections.within.family]) ; selected <- BestOne + first.in.family - 1 ;Selections <- c(Selections, selected)
        first.in.family <- last.in.family + 1 }
      if(length(numSelections.among.family) > 0) {
        Selections <- Selections[top.280.families[1:numSelections.among.family]]} else {
          Selections <- Selections[top.280.families]
        }
    }}
  
  # Extract phenotypes of new selections
  Selection.phenos<- progeny.phenos[Selections]
  new.pars.genval <- prog.genetic.values[Selections]
  Selection.EBvs <- progeny.blups[Selections]
  sorted.top192 <- sort(Selection.phenos,decreasing=T)
  #all.markers <- rbind(parent.markers,prog.markers[Selections,])
  
  if(gen==1){
    if(reduced){  all.phenos <- c(parentinfo$phenos,Selection.phenos)
    all.genetic.vals <- c(parentinfo$genetic.values,new.pars.genval)
    } else {all.phenos <- c(past.phenos$phenos,Selection.phenos)
    all.genetic.vals <- c(past.phenos$genetic.values,new.pars.genval)}
  } else { all.phenos <- c(parentinfo$all.phenos,Selection.phenos)
  all.genetic.vals <- c(parentinfo$all.genetic.vals,new.pars.genval)}
  ped <- full.ped[match(names(all.phenos),full.ped[,1]),]
  
  # Extract genotypes and parent ids of new selections
  #colnames(progenyinfo$genos.3d) <- names(progeny.phenos)
  new.parent.genos <- progenyinfo$genos.3d[,Selections,]
  numselections <- dim(new.parent.genos)[2]
  select.ids <- as.numeric(names(new.pars.genval))-numparents
  select.ped.ids <- as.numeric(names(new.pars.genval))
  if (selection.strategy == "self.test") {colnames(new.parent.genos) <- select.ped.ids} else {colnames(new.parent.genos) <- select.ped.ids}
  if(gen==1){
    all.markers1 <- cbind(parentinfo$genos.3d[mapinfo$available.Markers,,1],new.parent.genos[mapinfo$available.Markers,,1]) 
    all.markers2 <- cbind(parentinfo$genos.3d[mapinfo$available.Markers,,2],new.parent.genos[mapinfo$available.Markers,,2])
    all.markers <- abind(all.markers1,all.markers2,along=3)} else {
      all.markers1 <- cbind(parentinfo$all.markers[,,1],new.parent.genos[mapinfo$available.Markers,,1]) 
      all.markers2 <- cbind(parentinfo$all.markers[,,2],new.parent.genos[mapinfo$available.Markers,,2])
      all.markers <- abind(all.markers1,all.markers2,along=3)
    }
  # Calculate:
  #       Inbreding level of progeny/selections
  #       Genetic Variance of progeny and selections (Bulmer Effect)
  f.ped <- pedigree(label = full.ped[,1],sire = full.ped[,2],dam = full.ped[,3])
  pedigree.inbreeding <- inbreeding(f.ped)
  names(pedigree.inbreeding) <- full.ped[,1]
  progeny.inbreeding <- pedigree.inbreeding[prog.pedigree[,1]]
  selections.inbreeding <- progeny.inbreeding[Selections]
  progeny.gen.var <- var(prog.genetic.values[Selections])
  bulmer.effect <- var(prog.genetic.values[Selections]) - var(prog.genetic.values)
  
  #Calculate the number of deletrious alleles for each selection
  result <- sapply(rep(1:length(Selections),1),function(x) {out <-length(which(new.parent.genos[mapinfo$QTLSNP.loci,,1][,x]=="c")); out2 <- length(which(new.parent.genos[mapinfo$QTLSNP.loci,,2][,x]=="c"))
  outer <- out+out2})
  result <- unlist(result)
  
  if(relmatrix == "pedigree"){
    if(gen == 1){A <- as.matrix(getA(f.ped)) 
    } else {
      selection.ped <- crossdesign$selection.ped
      s.ped <- pedigree(label = selection.ped[,1],sire = selection.ped[,2],dam = selection.ped[,3])
      A <- as.matrix(getA(s.ped))}
    l <- match(names(Selection.phenos), rownames(A))
    rel.mat <- A[l,l]
  }
  
  if(relmatrix == "markers"){
    #m <- mapinfo$available.Markers
    #progeny1.markers <- new.parent.genos[m,,]
    allele1 <- all.markers[,,1]
    allele2 <- all.markers[,,2]
    
    numS <- ncol(allele1)
    g.rel.mat <- matrix(NA,nrow=numS,ncol=numS)
    yo <- lapply(1:(numS-1),FUN = function(x){
      com <- allele1[,x]
      com2 <- allele2[,x]
      ss <- x + 1
      b <- sapply(ss:ncol(allele1),FUN=function(x){
        (length(which(com == allele1[,x])) + length(which(com2 == allele2[,x])))/240
      })
    })
    diag(g.rel.mat) <- 1
    g.rel.mat[lower.tri(g.rel.mat)] <- unlist(yo)
    g.rel.mat <- t(g.rel.mat)
    g.rel.mat[lower.tri(g.rel.mat)] <- unlist(yo)
    
    rownames(g.rel.mat) <- colnames(allele1); 
    colnames(g.rel.mat) <- colnames(allele1)
    l <- which(rownames(g.rel.mat) %in% names(Selection.phenos))
    rel.mat <- g.rel.mat[l,l]
  }
  
  extraction.info<-list(relmat=rel.mat, delt.alleles=result, selections=Selections, bulmer.effect=bulmer.effect,
                        select.EBVs = Selection.EBvs, all.genetic.vals=all.genetic.vals,
                        selection.phenos=Selection.phenos,ped=ped, 
                        prog.inbred.level=progeny.inbreeding, select.inbred.level=selections.inbreeding,
                        genos.3d=new.parent.genos, num.parents=numselections,select.genval=new.pars.genval, fullped=full.ped,
                        par.ids=select.ids,select.ped.ids=select.ped.ids,all.markers=all.markers,all.phenos=all.phenos, cumulative.total=crossdesign$cumul.total)
  cat("The returned object is a list containing a matrix of phenotypic data with\n")
  cat("the specified heritability, a vector of unscaled true genetic values,\n")
  cat("to the same total variance as the phenotypic values, and a vector of\n" )
  cat("selected individuals, four per family, with the highest phenotype value.\n" )
  cat("Both phenotypes and true genetic values include values for parents first followed by progeny.\n")
  
  return(extraction.info)
  # # # # #   End of createPhenos() function # # # # # # # #
}
