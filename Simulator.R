####Load Libraries####
library(abind); library(doSNOW); library(synbreed)
library(pedigree); library(reshape2); library(parallel)

####Create Map#######
create.Map <- function(haldane=Haldane,save=Save,numMarkers=NumMarkers,
                      totalQTL=TotalQTL,numSNPQTL=NumSNPQTL,numchromos=num.chromos,
                      maplength=map.length,lowchromosizerange=0.25, 
                      highchromosizerange =.42, totalloci=total.loci,
                      signifdigits=2,rep.num=Rep) {
  
  # Determine sizes of 12 chromosomes
  # Size is based on the average chromo size multiplied by an upper/lower size range as specified above
  # The last chromosome size is determined by subtracting map length specificed by the summation of previous 11 chromosome sizes
  avgchromo <- maplength/numchromos 
  set.seed(32423) #A seed is set so that chromosome sizes are the exact same every time
  chromosizes <- floor(runif(numchromos, min=(avgchromo - lowchromosizerange*avgchromo), max=(avgchromo + highchromosizerange*avgchromo)))
  chromosizes[numchromos] <- maplength-sum(chromosizes[1:(numchromos-1)])
  
  # Determining the number of intervals that will be calculated and the number of loci per chromosome
  numintervals <- totalloci-1 #Specifying the # of intervals that will need to be calcuated
  #chromoloci <- vector("list",numchromos)
  #for(i in 1:numchromos){if(i == numchromos){chromoloci[i] <- totalloci - sum(unlist(chromoloci))
  #} else {chromoloci[i] <- floor((chromosizes[i]/maplength) * totalloci)}}
  #chromoloci <- unlist(chromoloci)
  chromoloci <- rep(totalloci/numchromos,numchromos)
  
  all.loci <- vector("list",numchromos)
  for (i in 1:length(chromoloci)){
    if (i==1){all.loci[[1]] <- chromoloci[i] } else {
      all.loci[[i]] <- chromoloci[[i]] + all.loci[[i-1]]}}
  chromo.loci <- unlist(all.loci)
  
  # Creating empty vectors 
  list1<-vector("list", totalloci) # List "chr#" character set for all loci
  list2<-vector("list", totalloci) # List locus# for all chromosomes
  intervals<-rep(NA, totalloci)    # Holds interval distance from each locus to the previous one
  positions<- rep(NA, totalloci)   # Holds the positions of each loci on a chromosome in cM
  locusnames<-rep(NA,totalloci)    # Holds each locus name for all of the loci (chr1_locus#)
  recfreqs<-rep(NA, totalloci)     # Holds Recombination frequency for each of the loci
  
  # Define map function used to create recombination frequencies:
  if(haldane) {mapfun<-function(x){(1-(exp(-2*x)))/2 }} else {mapfun <- function(x){ ( (exp(2*x) - exp(-2*x))/(exp(2*x) + exp(-2*x)) )/2}}
  
  ## Start with chromosome 1 - note that "interval" is distance from locus n to n-1, not n+1
  for(k in 1:numchromos){
    # For each chromosome list 1 will contain the number of loci for that chromsome
    list1[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]<-c(1:chromoloci[k])
    # list2 will hold the chr# that can later be pasted with list1 to make locus names
    list2[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]<-rep(paste("chr",k,sep=""))
    
    # Intervals are initally in cM but divided by 100 to be Morgans.  This is done so that it can be input to haldane function 
    # which uses Morgans to calculate Recfreqs
    intervals[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]<-c((round(runif((chromoloci[k]-1),min=(0.2*chromosizes[k]/chromoloci[k]), 
                                                                               max=(1.8*chromosizes[k]/chromoloci[k])),signifdigits))/100,5)
    # The last locus interval is determined by the mean interval of all previous loci
    intervals[sum(chromoloci[1:k])]<- mean(intervals[(1+sum(chromoloci[1:k-1])):(sum(chromoloci[1:k]-1))])
    
    #Recfreqs for each of the loci are determined by using the interval distance (for each locus) in Morgans and the haldane mapping function:
    #r=(1/2)*(1-e^(2*M))
    recfreqs[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])] <- mapfun(intervals[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])])
    
    #Since intervals are in Morgans, we multiply by 100 to turn it back into cM for positions (cM is unit for chromosome lengths)
    positions[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]<- round(100*cumsum(intervals[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]),digits=8)}
  
  ## paste together locus names from list2 and list1
  locusnames<- paste(list2,list1,sep="_")
  
  #Create map that has locusnames, distance between loci, empty vector of types (SNPQTL,rQTL, or Marker), MAFs, position on chromosome, & rec freqs
  map<-data.frame(chr=unlist(list2),loci=locusnames, dist= as.numeric(intervals),types= rep(NA,totalloci), MAF=rep(NA,totalloci), pos=as.numeric(positions),recfreqs=as.numeric(recfreqs))
  
  #Now sample from the map to specify SNPQTL & rQTL, the remainder are potential markers that can be used
  all.loci <- 1:totalloci  # vector that contains 1 through the number of all loci
  set.seed(123135)         # seed is set so that same SNPQTL positions are samples each time from all.loci
  SNPQTL <- sample(all.loci,numSNPQTL,replace=FALSE)
  set.seed(23423)     # Set seed so that snp qtl MAF are reproducible
  SNPQTLMAFs <- runif(numSNPQTL,min=.01,max=.02)
  map$types[SNPQTL] <- "snpqtl"   # Specify in the map data frame that these loci are snpqtl
  map$MAF[SNPQTL] <- SNPQTLMAFs   # Assign these loci the specificed minor allele frequencies generated by user
  
  
  numQTL <- totalQTL - numSNPQTL  # vector that contains the number of QTL which are not under dominance and have addtive effects
  set.seed(234253)                # seed is set so that the same loci are sampled to be rQTLs
  rQTL <- sample(all.loci[-SNPQTL],numQTL,replace=F)
  map$types[rQTL] <- "qtl"        # Specify in the map data frame that these loci are qtl
  
  set.seed(1323523)               # Set seed so that the same set of minor allele frequencies are drawn from beta distribution
  MAFs <- sample(rbeta(200000,.4,.4),totalloci,replace=F)
  map$MAF[-SNPQTL] <- sample(MAFs,length(map$MAF[-SNPQTL]),replace=F) # Assign minor allele frequencies to markers and qtl
  #map$MAF <- sample(MAFs,length(map$MAF),replace=F)
  map$types[-c(SNPQTL,rQTL)] <- "m"    # Specify in map data frame that all loci which are not qtl or snpqtl are markers
  Markers <- all.loci[-c(SNPQTL,rQTL)] # vector that contains all potential markers that can be used
  
  
  #Now save map output:
  datevalue <- date()
  datevector <- unlist(strsplit(datevalue,"\\s"))
  timevector <- unlist(strsplit(datevector[4],":"))
  newtime <- paste(timevector[1],timevector[2],sep="h")
  newdate <- paste(datevector[3],datevector[2],datevector[5],sep="")
  namestem <- paste(newdate,newtime,sep="_")
  
  if(save) {
    mapname<-paste(rep.num,namestem,"_map.txt")
    write.table(map,file=mapname, quote=F, row.names=F, col.names=T, sep="\t")
  }
  mapinfo<-list(genetic.map=map, 
                total.loci.num=totalloci,total.QTL.num=totalQTL, total.SNPQTL.num=numSNPQTL,
                QTLSNP.loci=sort(SNPQTL),rQTL.loci=sort(rQTL), available.Markers=Markers,
                last.locus.per.chrom=chromo.loci, 
                date.time=namestem)
  cat("The returned object is a list containing a matrix of 7 columns with chromosome #, locus names, intervals,\n")
  cat("types(marker,rqtl,snpqtl),MAFs, position on chromosome, and recombination fractions.\n")
  cat("A set of scalars is included with the numbers of total loci, total QTL, SNPs that are QTL, and random QTL\n")
  cat("Vectors identifying which loci are avaiable markers, which are SNP-QTL, and which are 'invisible' random QTL,\n")
  cat("and which locus is the last position per chromosome. Finally the namestem prefix used in saving text output files.\n")
  return(mapinfo)
  # End of function #
}

####Create Founder Parent Population####
create.Parents <- function(mapinfo=map,numparents=NumParents, max.delt.allele=MaxDeltAlleles, QTLsd=.25,save=Save,random=FALSE,rep.num=Rep){
  totalloci <- mapinfo$total.loci.num  # Specifies the total # of loci by pulling from the map object
  locusnames <- mapinfo$genetic.map$loci  # Specifies locus names from map object
  all.MAFs <- mapinfo$genetic.map$MAF     # Assign all.MAFs as the minor allele frequencies in the map object
  totalQTL <- mapinfo$total.QTL.num    # Total number of qtl pulled from map object
  QTLSNPs <- mapinfo$QTLSNP.loci      # A vector of the loci which are snpqtl
  numSNPQTL <- mapinfo$total.SNPQTL.num  # the number of snpqtl pulled from map object
  numQTL <- totalQTL - numSNPQTL  # the number of additive qtl (rQTL)
  all.Markers <- map$available.Markers      # A vector of loci which can potentially be markers
  num.markers <- length(all.Markers)  # the number of potential markers
  QTLoci <- mapinfo$rQTL.loci    # A vector of the loci which are additive QTL
  
  
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
  #D.marker genotypes are going to be assign to uppercase letters 
  #R.marker genotypes are going to be assigned lowercase letters

  ratio <- (max.delt.allele/2)
  # For each parent assign them alleles 
      for(par in 1:numparents){
      num.delt <- sample(0:ratio,1) # pick a number of deleterious alleles for allele 1 for this parent
      minor1 <- sample(available.delt1,num.delt)  # pick out of the avaiable deleterious alleles for this parent
      available.delt1 <- available.delt1[which(!available.delt1 %in% minor1 )] # this new list makes it so the next parent can't have any of these alleles at allele 1 positions
      major1 <- QTLSNPs[which(!QTLSNPs %in% minor1)]

      num.delt <- sample(0:ratio,1)
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
    }
  
  # Output the snpqtls genotypes to a matrix
  parentSNPgenos <- matrix(paste(parents[QTLSNPs,,1],parents[QTLSNPs,,2],sep=""),nrow=numSNPQTL,ncol=numparents)
  dimnames(parentSNPgenos) <- list(locusnames[QTLSNPs],parIDs)
  
  # Output the marker gentoypes to a matrix
  parent.markers <- matrix(paste(parents[all.Markers,,1],parents[all.Markers,,2],sep=""),nrow=num.markers,ncol=numparents)
  dimnames(parent.markers) <- list(locusnames[all.Markers],parIDs)
  

  result.parent <- sapply(rep(1:numparents,1),function(x) length(which(parentSNPgenos[,x]=="ac" | parentSNPgenos[,x]=="ca" | parentSNPgenos[,x]=="cc")))
  result.parent <- unlist(result.parent)
  # If save is true then then we save the parentSNPgenos object and parent.markers object
  if(save) {
    parfilename=paste(rep.num,"parent-SNPQTL-genotypes.txt")
    write.table(parentSNPgenos,parfilename, quote = F, row.names = T, col.names = T, sep= "\t")
    
    parfilename=paste(rep.num,"parent-marker-genotypes.txt")
    write.table(parent.markers,parfilename, quote = F, row.names = T, col.names = T, sep= "\t")    
  }
  
  ###IMP that genos.3d is first list element and numparents is second becuase that is what createTGV routine expects###
  parentinfo<-list(genos.3d=parents,delt.allele=result.parent,num.parents=numparents, parent.SNPQTL.matrix=parentSNPgenos, 
                   parent.Marker.matrix=parent.markers, parent.IDs=parIDs)
  cat("The returned object is a list containing a 3-D array of parent genotypes.\n")
  cat("A vector containing # of deleterious alleles per parent and scalar of total # of parents\n")
  cat("Two matrices, one including parent SNP genotypes and the other with parent marker loci, with indnivdiuals in columns.\n")
  cat("A vector of parent IDs.\n")
  return(parentinfo)
  # # # # #    End of function    # # # # # # # #
}

#####The above steps are done only one time, everything below is repeated for each generation####

####Create  TGV Genetic Value####
create.TGV <- function(data1=parents, mapinfo=map, crossdesign=cross.design, gen=gen,
                               Dom.coeff=dom.coeff, num.markers=NumMarkers,
                               save=Save, rep.num=Rep, prefix=prefix){
  
  dom.coeff <- Dom.coeff         # Dominance coeffecient that is set by user
  locusnames <- mapinfo$genetic.map$loci # The locus names pulled from the mab object
  QTLSNPs <- mapinfo$QTLSNP.loci     # vector of the loci which are snpqtl
  QTLSNPnum <- data1$genos.3d[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
  markers <- map$available.Markers# a list of all the markers pulled from map object
  
  set.seed(45234235) # seed is set making the loci which are selected as markers the same every time
  markers.select <- sort(sample(markers,num.markers,replace=F),decreasing=F)
  num.markers <- length(markers.select) # length of markers that were selected
  marker.select.genos <- data1$genos.3d[markers.select,,] # genotypes of the markers pulled from the current generation
  map.markers <- mapinfo$genetic.map[markers.select,c(1,6)]
  numQTL <- length(mapinfo$rQTL.loci) # the number of additive qtl
  if (gen==0){
    parIDs <- data1$parent.IDs
  } else {
  parIDs <- crossdesign$parent.IDs} # the parent ids pulled from the current generation
  
  
  # Assign each snpqtl a value
  A <- 1    # Major allele is assigned 1  
  a <- -1000  # Minor allele is assigned -20
  
  # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
  numSNPQTL <- mapinfo$total.SNPQTL.num # the number of loci which are snpqtl
  numparents <- length(parIDs) # the number of parents
  QTLSNPvalues <-matrix(NA,nrow=numSNPQTL,ncol=numparents) # matrix to hold snpqtl values
  marker.values <- matrix(NA,nrow=num.markers,ncol=numparents) # matrix to hold marker values
  Capital.genotypes <- vector()
  Lowercase.genotypes <- vector()
  for (i in 1:26){
    Capital.genotypes <- c(Capital.genotypes,paste(LETTERS[i],LETTERS, sep=""))
    Lowercase.genotypes <-  c(Lowercase.genotypes,paste(letters[i],letters, sep=""))
  }

  dead <- NULL
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
        QTLSNPvalues[QTLSNPcc,i] <- a
        QTLSNPvalues[QTLSNPac,i] <- (A*dom.coeff)
        QTLSNPvalues[QTLSNPca,i] <- (A*dom.coeff)}
    markers.aa <- which(marker.select.genos[,i,1] %in% Capital.genotypes & marker.select.genos[,i,2] %in% Capital.genotypes)
    markers.cc <- which(marker.select.genos[,i,1] %in% Lowercase.genotypes & marker.select.genos[,i,2] %in% Lowercase.genotypes)
    markers.ac <-which(marker.select.genos[,i,1] %in% Capital.genotypes & marker.select.genos[,i,2] %in% Lowercase.genotypes)
    markers.ca <- which(marker.select.genos[,i,1] %in% Lowercase.genotypes & marker.select.genos[,i,2] %in% Capital.genotypes)
    
    marker.values[markers.aa,i] <- "0"
    marker.values[markers.cc,i] <- "2"
    marker.values[markers.ac,i] <- "1"
    marker.values[markers.ca,i] <- "1"
  dead <- c(dead,QTLSNPcc)
    }

  marker.values <- t(marker.values)
  colnames(marker.values) <- markers.select
  rownames(marker.values) <- parIDs
  
  # Convert the 'invisible' rQTL genotypes to numeric matrices, merge the alleles to paired values also
  parQTLallele1 <- matrix(as.integer(data1$genos.3d[mapinfo$rQTL,,1]),nrow=numQTL,ncol=numparents)
  colnames(parQTLallele1) <- c(parIDs)
  parQTLallele2 <- matrix(as.integer(data1$genos.3d[mapinfo$rQTL,,2]),nrow=numQTL,ncol=numparents)
  colnames(parQTLallele2) <- c(parIDs)
  QTLvalues     <- matrix(paste(parQTLallele1,parQTLallele2,sep=","),nrow=numQTL,ncol=numparents)
  dimnames(QTLvalues)<-list(locusnames[mapinfo$rQTL],parIDs)
  
  # Genetic values of progeny
  geneticvals <- colSums(QTLSNPvalues) + colSums(parQTLallele1) + colSums(parQTLallele2)
  
  if(save) {
    pedfilename=paste(rep.num,prefix,".txt", sep="")
    write.table(geneticvals,pedfilename, quote = F, row.names = F, col.names = T, sep="\t")
  }
  
  TGV <- list(genetic.values=geneticvals, SNP.value.matrix=QTLSNPvalues, markers.matrix=marker.values,
              marker.loci=markers.select, marker.map=map.markers, dead=dead)
  cat("The returned object is a vector containing genetic values for individuals\n")
  cat("Two matrices, one including SNP genetic values assigned and the other with marker loci, indnivdiuals in columns.\n")
  cat("A vector of loci which were selected as markers and a matrix of the genetic map for those markers.\n")
  return(TGV)
}

####Create Unscaled Phenos####
create.Unscaled.Phenos <- function(h2=Heritability, E.var= Env.var, gen=gen, TGV=parents.TGV, folder=Rep,namestem=map$date.time ) {
  
  geneticvals<-TGV$genetic.values
  totalindiv<-length(geneticvals)
  phenos<-vector(length=totalindiv)
  trueh<-vector()
  
  if(gen==0){
    for (i in 1:10000){
      set.seed(2342342)
      E <- rnorm(totalindiv, 0, (sqrt((i)*(1-h2))))# for U-shape MAF and 0 dominance 192 parents
      phenos <- round(geneticvals+E,4)
      trueh <- round(var(geneticvals)/var(phenos),2)
      print(trueh)
      print(i)
      E.var <- i
      if (trueh==.3) break
    }} else {
      E <- rnorm(totalindiv, 0, (sqrt(E.var*(1-h2))))
      phenos <- round(geneticvals+E,4)
      trueh <- var(geneticvals)/var(phenos)
      print(trueh)
    }
  
  pheno.info<-list(phenos=phenos, genetic.values=geneticvals, E.var = E.var)
  cat("Returned object contains a vector containing phenotypic values for individuals,\n")
  cat("vector with corresponding genetic values, and multiplier used for simulating an initial trait heritability.\n")
  return(pheno.info)
}

####Create Cross Design & Pedigree's####
createCrossDesign <- function(mapinfo=map,parentinfo=parents, num.prog.percross=prog.percross, cross.file.input=crossfileinput, gen=1, saveped=F,
                              crossfile=CrossFile1, RandomMating=random.mating , Self=self, SP.MateSelect=sp.mateselect, MP.MateSelect=mp.mateselect){
  
  coancest=NULL
  NumParents=parentinfo$num.parents # Pull number of parents from either parents (if gen is 1) or last phenos if gen>1
  
  # Make Cross Design File Depending on Strategy
  if ( RandomMating ){
    crossdesign= matrix(sample(names(parentinfo$selections),replace=FALSE, size=NumParents),nrow = NumParents/2,ncol=2)
    crossdesign <- cbind(crossdesign,rep(num.prog.percross,length(crossdesign[,1])))
    }else if (Self){
    crossdesign=matrix(data=c(rep(1:NumParents,times=2),rep(NA,NumParents)),nrow=NumParents,ncol=3)
    crossdesign[,3] = rep(num.prog.percross)
  }  else if (SP.MateSelect ) { 
    par1 <- vector(); par2 <- vector()
    relationshipmatrix <- parentinfo$relmat[which(colnames(parentinfo$relmat) %in% parentinfo$select.ped.ids),which(colnames(parentinfo$relmat) %in% parentinfo$select.ped.ids)]
    E <- upper.tri(relationshipmatrix); relationshipmatrix[E] <- NA
    coancestry <- seq(.01,.995,.005)
    for(j in 1:length(coancestry)){
    new <- subset(melt(as.matrix(t(relationshipmatrix))), value < coancestry[j])
    matches1 <- match(new[,1],names(parentinfo$select.EBVs))
    matches2 <- match(new[,2],names(parentinfo$select.EBVs))
    a <- parentinfo$select.EBVs[matches1]
    b <- parentinfo$select.EBVs[matches2]
    new$mean.bv <- (a+b)/2
    new <- new[order(new$mean.bv,decreasing=T),]
    coancest <- coancestry[j]
    crossdesign <- matrix(,nrow=(NumParents/2),ncol=3)
    sib1 <- sample(colnames(relationshipmatrix)[seq(1,NumParents,2)],NumParents/2)
    sib2 <- colnames(relationshipmatrix)[seq(2,NumParents,2)]
      for (i in 1:(length(parentinfo$selections)/2)) {
        newww <- new[new[,2] == sib1[i],]
        newww.1 <- new[new[,1] == sib1[i],]
        c.new <- rbind(newww,newww.1)
        id <- sib1[i]
        l <- sib2
        cc.new <- c.new[apply(c.new[,1:2], MARGIN=1, function(x) any(x %in% l)),]
        okay <- cc.new[order(-cc.new[,4]),]
        c <- unlist(okay[1,1:2])
        crossdesign[i,1] <- c[1]
        crossdesign[i,2] <- c[2]
        crossdesign[i,3] <- num.prog.percross
        new <- new[new$X1!=c[1] & new$X2!=c[1],]
        new <- new[new$X1!=c[2] & new$X2!=c[2],]
      }
    if (!anyNA(crossdesign[,1]) ) break }} else if  (MP.MateSelect) {
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
        if(min(minimum) >= 90) break }
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
    progeny <- rep(num.prog.percross, length(par1))
    crossdesign=matrix(data=c(par1,par2,progeny),nrow=length(progeny),ncol=3)
    if (!anyNA(crossdesign[,1])) break} else if (cross.file.input) {
    colClasses = c("integer", "integer", "integer")
    crossdesign = as.matrix(read.table(crossfile,header=F,sep="\t", colClasses))
  }
  
  numcrosses <- nrow(crossdesign) # Total number of crosses to be made
  totalprogeny <- sum(as.numeric(crossdesign[,3])) # Total number of progeny to be created
  par1list<-vector(); par2list<-vector() # Vectors to hold the parents of each cross fro all progeny
  par1id<-vector(); par2id<-vector()
  generation<- rep(gen,totalprogeny) # Vector of generation number to be attached to pedigree
  
  if (gen==1){
    totalindiv   <- totalprogeny+length(unique(crossdesign[,1:2])) # Number of total individuals including parents
    cumul.total <- length(unique(crossdesign[,1:2])) + totalprogeny
    #Now create ID names for each progeny which begins with the ID of last parent + 1
    indivIDs<-c(seq((parentinfo$num.parents+1),(parentinfo$num.parents)+(totalprogeny)))
    
    for (m in 1:numcrosses){ #For each cross
      crossprog<-as.numeric(crossdesign[m,3]) #Store the number of progeny to create
      par1<-as.character(crossdesign[m,1])    #Identify parent one for that cross
      par2<-as.character(crossdesign[m,2])    #Identify parent two for that cross
      par1list<-c(par1list,rep(par1,crossprog)) #Rep the ID for parent 1 for however many progeny
      par2list<-c(par2list,rep(par2,crossprog))} #Rep the ID for parent 2 for however many progeny
    
    uni <- length(unique(c(crossdesign[,1],crossdesign[,2]))) # For all crosses what is length of unique parents
    #Create parent pedigree by reping 1 through number of unique parents
    #Since this is first generation Par1 and Par2 are uknown and generation is 0
    parentped <- cbind("ID"=rep(1:uni,1),"Par1" = rep(0,uni), "Par2"= rep(0,uni),"gener"=rep(0,uni))
    
    #Create progeny pedigree by combining indiv IDs, parents identefied from par1list and par2list, and the generation which is 1
    progped<-cbind("ID"=indivIDs[1:length(indivIDs)], "Par1"=par1list[1:length(par1list)],"Par2"=par2list[1:length(par1list)],"gener"=generation[1:length(generation)])
    
    #The full pedigree up to this point contains the parent pedigree and progeny pedigree
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
      par1list<-c(par1list,rep(par1,crossprog))
      par2list<-c(par2list,rep(par2,crossprog))}
    progped<-cbind("ID"=indivIDs[2:length(indivIDs)], "Par1"=par1list,"Par2"=par2list,"gener"=generation)
    full.ped <- rbind(parentinfo$fullped,progped)
    selection.pedigree <- rbind(parentinfo$ped,progped)
  }
  
  if(saveped) {
    pedfilename=paste(rep.num,mapinfo$date.time,prefix,"-pedigree.txt",sep="")
    write.table(progped,pedfilename, quote = F, row.names = F, col.names = T, sep=" ")}
  
  ###Imp that genos.3d is first in list item, and total progeny is next list item bc in createtgv routine it expects genos.3d as first element and total prog as second
  progeny<-list(crossdesign=crossdesign, progeny.pedigree=progped, full.pedigree= full.ped, selection.ped=selection.pedigree,
                total.progeny.number=totalprogeny, numparents=parentinfo$num.parents, numcrosses=numcrosses,crossprog=crossprog,
                parent.IDs=progped[,1], namestem=mapinfo$date.time,cumul.total=cumul.total, coancestry.threshold=coancest)
  cat("Returned object contains matrix of cross design with 3 columns: parent 1, parent 2, and number of progeny.\n")
  cat("Three matrices - progeny, selections, and full pedigree which have 4 columns: Indiv ID, par1, par2, generation #.\n")
  cat("Scalars including total number of progeny, parents, crosses made, and progeny per cross.\n")
  cat("A vector of progeny IDs is returned. Finally, the cumulative total # of prog, name stamp used to save files,\n")
  cat("and the coancestry threshold used to decide which crosses should be made is returned.\n")
  return(progeny)}

####Create Make Crosses####
makecrosses <- function(x,mapinfo=map,parentinfo=parents, cumul.total=parentinfo$cumulative.total, num.prog.percross=prog.percross,
                        crossfile=CrossFile1, gen=gen, savegenos=Save,saveped=Save,prefix=prefix,rep.num=Rep, cross.file=as.matrix(crossfile),
                        sbpm=SBPM, selfl=SELFL,selfp=SELFP, sbpw=SBPW, MATESELECT=mate.select, selftest=SELFTEST){
  
  library(abind)
  NumParents=parentinfo$num.parents # Pull number of parents from either parents (if gen is 1) or last phenos if gen>1
  chromoloci <- mapinfo$last.locus.per.chrom   # The last number of each loci for all chromsomes
  QTLSNPs      <- mapinfo$QTLSNP.loci # Vector of loci which are snpqtl
  

  par1<-match(x[1],colnames(parentinfo$genos.3d)) # assigns par1 to be the first parent in crossdesign matrix
  par2<-match(x[2],colnames(parentinfo$genos.3d)) # assigns par2 to be the second parent in the crossdesign matrix
  crossprog<-as.numeric(x[3]) #assigns number of progeny to be the third column for cross "x"
  
  # Create empty matrix to hold gametes
  # dimensions are (total # of loci) x  (# of cross progeny)
  # rownames are the loci names
  gametes1<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
  gametes2<-matrix(rep(NA,mapinfo$total.loci.num*crossprog),nrow=mapinfo$total.loci.num,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
  
  
  r <- vector("list",crossprog)
  par1.alleles <- (parentinfo$genos.3d[,par1,]) # Store alleles from parent 1
  par2.alleles <- (parentinfo$genos.3d[,par2,]) # Store alleles from parent 2

  # Simulating 1st set of gametes for all progeny of given cross
  {
  #For each chromsome sample loci (1.5 x the number of progeny, with replacement) which will have a recombination event 
  chr1.r <- sort(sample(1:chromoloci[1],crossprog*1.5,replace=T))
  chr2.r <- sort(sample((chromoloci[1]+1):chromoloci[2],crossprog*1.5,replace=T))
  chr3.r <- sort(sample((chromoloci[2]+1):chromoloci[3],crossprog*1.5,replace=T))
  chr4.r <- sort(sample((chromoloci[3]+1):chromoloci[4],crossprog*1.5,replace=T))
  chr5.r <- sort(sample((chromoloci[4]+1):chromoloci[5],crossprog*1.5,replace=T))
  chr6.r <- sort(sample((chromoloci[5]+1):chromoloci[6],crossprog*1.5,replace=T))
  chr7.r <- sort(sample((chromoloci[6]+1):chromoloci[7],crossprog*1.5,replace=T))
  chr8.r <- sort(sample((chromoloci[7]+1):chromoloci[8],crossprog*1.5,replace=T))
  chr9.r <- sort(sample((chromoloci[8]+1):chromoloci[9],crossprog*1.5,replace=T))
  chr10.r <- sort(sample((chromoloci[9]+1):chromoloci[10],crossprog*1.5,replace=T))
  chr11.r <- sort(sample((chromoloci[10]+1):chromoloci[11],crossprog*1.5,replace=T))
  chr12.r <- sort(sample((chromoloci[11]+1):chromoloci[12],crossprog*1.5,replace=T))
  
  #Now for each chromsome sample the progeny in that cross which will be assigned the recombination event
  #Since 1 progeny can have up to 2 recombination events on the same chrosome we sample twice without replacement
  indiv1 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv2 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv3 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv4 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv5 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv6 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv7 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv8 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv9 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv10 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv11 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv12 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  }
  
  # Finally for all progeny in this cross we create a set of switches using the loci and progeny involved in 
  #each recombination event for the 12 chromosomes
  #Ex if there are 20 progeny per cross in an organism that has 12 chromsomes;
  #   Then 317 total recombination events spread among 20 progeny is a recombination rate of 317/240 = 1.32
  for (i in 1:(crossprog*1.5)){
    #This stores the recombination event for a given loci across all progeny
    a <- indiv1[i]; b <- indiv2[i]; c <- indiv3[i]; d <- indiv4[i]; e <- indiv5[i]
    f <- indiv6[i]; g <- indiv7[i]; h <- indiv8[i]; i <- indiv9[i]; j <- indiv10[i]
    k <- indiv11[i]; l <- indiv12[i]
    r[[a]] <- unique(sort(c(r[[a]],chr1.r[i]))); r[[b]]<- unique(sort(c(r[[b]],chr2.r[i])))
    r[[c]]<- unique(sort(c(r[[c]],chr3.r[i]))); r[[d]]<- unique(sort(c(r[[d]],chr4.r[i])))
    r[[e]]<- unique(sort(c(r[[e]],chr5.r[i]))); r[[f]]<- unique(sort(c(r[[f]],chr6.r[i])))
    r[[g]]<- unique(sort(c(r[[g]],chr7.r[i]))); r[[h]]<- unique(sort(c(r[[h]],chr8.r[i])))
    r[[i]]<- unique(sort(c(r[[i]],chr9.r[i]))); r[[j]]<- unique(sort(c(r[[j]],chr10.r[i])))
    r[[k]]<- unique(sort(c(r[[k]],chr11.r[i]))); r[[l]]<- unique(sort(c(r[[l]],chr12.r[i])))
  }
  #Ex.output:  print(r[[20]])
  #[1]  2239  2306  8384 12316 18194 20577 23283 25262 30852 31882 38367 43465 46694 52381 55888
  #These are the loci associated with this progeny which will have a recombination event
  
  #Finally for 1 through all cross prog we populate gametes 1 using either the first or second allele from parent 1
  for(i in 1:crossprog) {
    allele <- sample(1:2,1) #Start out by sampling either allele 1 or 2 for parent
    end <- chromoloci[12]  #Specify that the last loci the end of last chromosome
    z <-1
    recombination.spots <- r[[i]]
    for (each in 1:length(recombination.spots)) {
      if (each < length(recombination.spots)){
        gametes1[z:recombination.spots[each],i] <- par1.alleles[z:recombination.spots[each],allele]
      } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
      if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
      z <- recombination.spots[each] +1
    }}
  
  # Simulating 2nd set of gametes for all progeny of given cross
  {
  chr1.r <- sort(sample(1:chromoloci[1],crossprog*1.5,replace=T))
  chr2.r <- sort(sample((chromoloci[1]+1):chromoloci[2],crossprog*1.5,replace=T))
  chr3.r <- sort(sample((chromoloci[2]+1):chromoloci[3],crossprog*1.5,replace=T))
  chr4.r <- sort(sample((chromoloci[3]+1):chromoloci[4],crossprog*1.5,replace=T))
  chr5.r <- sort(sample((chromoloci[4]+1):chromoloci[5],crossprog*1.5,replace=T))
  chr6.r <- sort(sample((chromoloci[5]+1):chromoloci[6],crossprog*1.5,replace=T))
  chr7.r <- sort(sample((chromoloci[6]+1):chromoloci[7],crossprog*1.5,replace=T))
  chr8.r <- sort(sample((chromoloci[7]+1):chromoloci[8],crossprog*1.5,replace=T))
  chr9.r <- sort(sample((chromoloci[8]+1):chromoloci[9],crossprog*1.5,replace=T))
  chr10.r <- sort(sample((chromoloci[9]+1):chromoloci[10],crossprog*1.5,replace=T))
  chr11.r <- sort(sample((chromoloci[10]+1):chromoloci[11],crossprog*1.5,replace=T))
  chr12.r <- sort(sample((chromoloci[11]+1):chromoloci[12],crossprog*1.5,replace=T))
  indiv1 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv2 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv3 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv4 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv5 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv6 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv7 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv8 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv9 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv10 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv11 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv12 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  }
  r <- vector("list",crossprog)
  for (i in 1:(crossprog*1.5)){
    a <- indiv1[i]; b <- indiv2[i]; c <- indiv3[i]; d <- indiv4[i]; e <- indiv5[i]
    f <- indiv6[i]; g <- indiv7[i]; h <- indiv8[i]; i <- indiv9[i]; j <- indiv10[i]
    k <- indiv11[i]; l <- indiv12[i]
    r[[a]] <- unique(sort(c(r[[a]],chr1.r[i]))); r[[b]]<- unique(sort(c(r[[b]],chr2.r[i])))
    r[[c]]<- unique(sort(c(r[[c]],chr3.r[i]))); r[[d]]<- unique(sort(c(r[[d]],chr4.r[i])))
    r[[e]]<- unique(sort(c(r[[e]],chr5.r[i]))); r[[f]]<- unique(sort(c(r[[f]],chr6.r[i])))
    r[[g]]<- unique(sort(c(r[[g]],chr7.r[i]))); r[[h]]<- unique(sort(c(r[[h]],chr8.r[i])))
    r[[i]]<- unique(sort(c(r[[i]],chr9.r[i]))); r[[j]]<- unique(sort(c(r[[j]],chr10.r[i])))
    r[[k]]<- unique(sort(c(r[[k]],chr11.r[i]))); r[[l]]<- unique(sort(c(r[[l]],chr12.r[i])))}
  for(i in 1:crossprog) {
    allele <- sample(1:2,1)
    end <- chromoloci[12]
    z <-1
    recombination.spots <- r[[i]]
    for (each in 1:length(recombination.spots)) {
      if (each < length(recombination.spots)){
        gametes2[z:recombination.spots[each],i] <- par2.alleles[z:recombination.spots[each],allele]
      } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
      if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
      z <- recombination.spots[each] +1
    }}
  gametes <- list(gametes1,gametes2)
  array.out <- abind(gametes, along=3)
  out <- list(as.array(array.out))
}

####Extract Genos####
create.genos <- function(mapinfo=map,crossdesign=cross.design$crossdesign,progeny=progeny1) {
  numprog <- crossdesign[1,3] * 2 * nrow(crossdesign) #Number of progeny * # of alleles 
  x <- seq(1,numprog,crossdesign[1,3])
  y <- seq(1,(nrow(crossdesign)*2),2)
  z <- seq(2,(nrow(crossdesign)*2),2)
  x.y <- x[y]
  x.z <- x[z]
  d <- NULL; e <- NULL
  for( i in 1:length(x.y)) {
    d <- c(d,seq(x.y[i],(x.y[i]+19)))
    e <- c(e,seq(x.z[i],(x.z[i]+19)))
  }
  yo.1 <- data.frame(progeny[c(1:nrow(crossdesign))])[,d]
  yo.2 <- data.frame(progeny[c(1:nrow(crossdesign))])[,e]
  geno <- list(yo.1,yo.2)
  genos.3d <- abind(geno, along=3)

  out <- list(genos.3d=genos.3d)
  return(out)}

####Extract Selections####
Extract.Selections <- function(parents.tgv=parents.TGV, past.phenos=parents.phenos, gen=1,
                               progenyinfo=progeny1, progenyTGV=progeny1.TGV, progenyphenos=progeny1.phenos,
                               #Selection types and Methods
                               sbpm=F,selfl=F, selfp=F, MATESELECT=mateselect,
                               selftest=SELFTEST, gblup=Gblup,ablup=Ablup, hblup=F, p=P,
                               #Constants
                               mapinfo=map,h2=Heritability, crossdesign.inuse=cross.design,Selfed.progeny=selfed.progeny,
                               numSelections=NumSelections, prefix=prefix,rep.num=Rep) {
  
  # Extracting crossdesign, number of parents, complete(full) & progeny pedigree from progenyinfo object
  crossdesign <- crossdesign.inuse$crossdesign
  prog.pedigree <- crossdesign.inuse$progeny.pedigree[which(crossdesign.inuse$progeny.pedigree[,2] != crossdesign.inuse$progeny.pedigree[,3]),]
  numparents <- crossdesign.inuse$numparents
  full.ped <- crossdesign.inuse$full.pedigree
  selection.ped <- data.frame(crossdesign.inuse$selection.ped)
  
  # Extracting parent markers, progeny tgv & markers, as well as creating marker map from progenyTGV object
  if (gen==1){parent.markers <- parents.tgv$markers.matrix} else {parent.markers <- past.phenos$all.markers}
  prog.markers <- progenyTGV$markers.matrix
  progeny.phenos <- progenyphenos$phenos
  prog.genetic.values <- progenyTGV$genetic.values
  map.markers <- progenyTGV$marker.map
  colnames(map.markers) <- c("chr","pos")
  map.markers$chr <- as.character(map.markers$chr)
  
  # Generate pedigree using synbreed if 1st gen otherwise just use full.ped object
  # Create 2 objects which hold all markers and all phenotypes respectively
  if (gen==1) {
      ped <- create.pedigree(ID=full.ped[,1],full.ped[,2],full.ped[,3],full.ped[,4])
      all.phenos <- c(past.phenos$phenos,progeny.phenos)
      } else{
        ped <- create.pedigree(ID=selection.ped[,1],selection.ped[,2],selection.ped[,3],selection.ped[,4])
        all.phenos <- c(past.phenos$all.phenos,progeny.phenos)
      }
  ped <- ped[order(as.numeric(ped$ID)),]
      all.markers <- rbind(parent.markers,prog.markers)
  
  # Pepare for Selection Extraction
  # Make Selections based on type and method
  # Types: SBPM, SBPW, SELFP, SELFL
  # Methods: Phenotype, ABLUP, GBLUP
  first.in.family <- 1
  Selections <- vector()
  Capital.genotypes <- vector(); Lowercase.genotypes <- vector()
  for (i in 1:26){
    Capital.genotypes <- c(Capital.genotypes,paste(LETTERS[i],LETTERS,LETTERS[i],LETTERS, sep=""))
    Lowercase.genotypes <-  c(Lowercase.genotypes,paste(letters[i],letters, letters[i], letters, sep=""))}
  if (selftest) {
    for (each.cross in 1:length(crossdesign[,1])){
      yo <- vector(mode="list",length=as.numeric(crossdesign[1,3]))
      yoo <- vector(mode="list",length=crossdesign[1,3])
      yo.1 <- vector(mode="list",length=crossdesign[1,3])
      yoo.1 <- vector(mode="list",length=crossdesign[1,3])
      parent1 <- as.numeric(crossdesign[each.cross,1])
      parent2 <- as.numeric(crossdesign[each.cross,2])
      print(parent2)
      num.prog <- seq(first.in.family,first.in.family + crossdesign[each.cross,3] -1)
      if (gen==1){
      parent1.alleles <- progenyinfo$genos.3d[as.numeric(names(Selfed.progeny[[parent1]]$top.indiv)),num.prog,1]
      parent2.alleles <- progenyinfo$genos.3d[as.numeric(names(Selfed.progeny[[parent2]]$top.indiv)),num.prog,2]
      top.parent.1 <- Selfed.progeny[[parent1]]$top.indiv
      top.parent.2 <- Selfed.progeny[[parent2]]$top.indiv
      } else {
      parent1.alleles <-  progenyinfo$genos.3d[as.numeric(names(Selfed.progeny[[which(names(Selfed.progeny) %in% parent1)]]$top.indiv)),num.prog,1]
      parent2.alleles <-  progenyinfo$genos.3d[as.numeric(names(Selfed.progeny[[which(names(Selfed.progeny) %in% parent2)]]$top.indiv)),num.prog,2]  
      top.parent.1 <- Selfed.progeny[[which(names(Selfed.progeny) %in% parent1)]]$top.indiv
      top.parent.2 <- Selfed.progeny[[which(names(Selfed.progeny) %in% parent2)]]$top.indiv
      }
    for(i in 1:20){
      yo[[i]] <- which(parent1.alleles[,i]==top.parent.1)
      yoo[[i]] <- sum(Selfed.progeny[[which(names(Selfed.progeny) %in% parent1)]]$most.imp[yo[[i]]])
      yo.1[[i]] <- which(parent2.alleles[,i]==top.parent.2)
      yoo.1[[i]] <- sum(Selfed.progeny[[which(names(Selfed.progeny) %in% parent2)]]$most.imp[yo.1[[i]]])
    }
      total <- mapply(function(x1,y1) x1[[1]] + y1[[1]], yoo, yoo.1)
     Selections <- c(Selections, num.prog[which(total %in% sort(total,decreasing=F)[1:2])])
    #Selections <- c(Selections,num.prog[which(lengths(mapply(c,yo,yo.1,SIMPLIFY=F), use.names = FALSE) %in% sort(lengths(mapply(c,yo,yo.1,SIMPLIFY=F), use.names = FALSE),decreasing=T)[1:2])][1:2])
    first.in.family <- max(num.prog) + 1
    }
  Selections.blup <- Selections
  }

if (ablup){
  phenos.1 <- as.matrix(all.phenos); colnames(phenos.1) <- "Vol"
  gp <- create.gpData(pheno=phenos.1, pedigree=ped)
  ADD.Rel.Mat <- kin(gp,ret="add") ; rel.mat <- data.frame(ADD.Rel.Mat); colnames(rel.mat) <- ped[,1]
  PBLUP <- gpMod(gp, model="BLUP", trait = "Vol", kin=ADD.Rel.Mat, markerEffects=F,predict=F)
    g <- match(names(progenyphenos$phenos),rownames(PBLUP$fit$predicted)) 
    progeny.blups <- PBLUP$fit$predicted[g] ; names(progeny.blups) <- names(progenyphenos$phenos) 
    }

if (hblup){
  phenos.1 <- as.matrix(all.phenos); colnames(phenos.1) <- "Vol"
  gp <- create.gpData(pheno=phenos.1,geno=all.markers, pedigree=ped); gp <- codeGeno(gp,label.heter = "1")
  ADD.Rel.Mat <- kin(gp,ret="add") ; rel.mat <- data.frame(ADD.Rel.Mat); colnames(rel.mat) <- ped[,1]
  G.Rel.Mat <- kin(gp, ret="realized"); H.Rel.Mat <- G.Rel.Mat + ADD.Rel.Mat
  HBLUP <- gpMod(gp, model="BLUP", trait = "Vol", kin=H.Rel.Mat, markerEffects=F,predict=F)
  g <- match(names(progenyphenos$phenos),rownames(HBLUP$fit$predicted))
  progeny.blups <- HBLUP$fit$predicted[g] ; names(progeny.blups) <- names(progenyphenos$phenos)
}

if (gblup){
  phenos.1 <- as.matrix(all.phenos); colnames(phenos.1) <- "Vol"
  gp <- create.gpData(pheno=phenos.1,geno=all.markers, pedigree=ped); gp <- codeGeno(gp,label.heter = "1")
  G.Rel.Mat <- kin(gp,ret="realized"); A.Rel.Mat <- kin(gp,ret="add"); rel.mat <- data.frame(A.Rel.Mat); colnames(rel.mat) <- ped[,1]
  GBLUP <- gpMod(gp, model="BLUP", trait = "Vol", kin=G.Rel.Mat, predict=T, markerEffects=F)
  g <- match(names(progenyphenos$phenos),rownames(GBLUP$fit$predicted))
  progeny.blups <- GBLUP$fit$predicted[g]; names(progeny.blups) <- names(progenyphenos$phenos)
  }

if(sbpw | selfl | MATESELECT){
  #Selections within each family
  for(family in 1:length(crossdesign[,3])){
    num.offspring <- as.numeric(crossdesign[family,3])
    last.in.family <- num.offspring + first.in.family - 1
    if(p) {temp <- (progeny.phenos[first.in.family:last.in.family]) ; sorted <-sort(temp,decreasing=TRUE)
    phenos.1 <- as.matrix(all.phenos); colnames(phenos.1) <- "Vol"
    gp <- create.gpData(pheno=phenos.1, pedigree=ped); A.Rel.Mat  <- kin(gp,ret="add"); rel.mat <-data.frame(A.Rel.Mat); colnames(rel.mat) <- ped[,1]
    } else {
      temp <- (progeny.blups[first.in.family:last.in.family]) ; sorted <-sort(temp,decreasing=TRUE)}
    for (i in 1:2){BestOne <- which(temp==sorted[i]) ; selected <- BestOne + first.in.family - 1 ;Selections <- c(Selections, selected)}
    first.in.family <- last.in.family + 1 }} else if (sbpm) {
      sorted.EBVs <- sort(progeny.blups, decreasing = T)
      for(i in 1:numSelections){
        selected <- which(progeny.blups==sorted.EBVs[i]); Selections <- c(Selections,selected)}}


# Extract phenotypes of new selections
Selection.phenos<- progeny.phenos[Selections]
if (selftest){
  Selection.EBvs <- progeny.blups[which(as.numeric(names(progeny.blups)) %in% Selections.blup)]
} else if(p){
  Selection.EBvs <- Selection.phenos
} else {
  Selection.EBvs <- progeny.blups[Selections]
  }
#sorted.top192 <- sort(Selection.phenos,decreasing=T)
all.markers <- rbind(parent.markers,prog.markers[Selections,])
if(gen==1){
all.phenos <- c(past.phenos$phenos,Selection.phenos) } else {
  all.phenos <- c(past.phenos$all.phenos,Selection.phenos)
}
ped <- data.frame(full.ped[as.numeric(names(all.phenos)),])

# Extract genotypes and parent ids of new selections
new.parent.genos <- progenyinfo$genos.3d[,Selections,]
#if (selftest) {colnames(new.parent.genos) <- which(as.numeric(names(progeny.blups)) %in% Selections)} else {colnames(new.parent.genos) <- names(Selections)}
numselections <- dim(new.parent.genos)[2]
new.pars.genval <- prog.genetic.values[Selections]
select.ids <- as.numeric(names(new.pars.genval))-numparents
select.ped.ids <- as.numeric(names(new.pars.genval))  
if (selftest) {colnames(new.parent.genos) <- select.ped.ids} else {colnames(new.parent.genos) <- names(Selections)}


# Calculate: 
#       Inbreding level of progeny/selections
#       Genetic Variance of progeny and selections (Bulmer Effect)
pedigree.inbreeding <- calcInbreeding(selection.ped[,1:3])
names(pedigree.inbreeding) <- selection.ped[,1]
progeny.inbreeding <- pedigree.inbreeding[prog.pedigree[,1]]
selections.inbreeding <- progeny.inbreeding[Selections]

dead <- sapply(rep(1:length(prog.genetic.values),1),function(x) {out <-which(progenyinfo$genos.3d[mapinfo$QTLSNP.loci,x,1]=="c" & progenyinfo$genos.3d[mapinfo$QTLSNP.loci,x,2]=="c")})
dead <- unlist(dead)
progeny.gen.var <- var(prog.genetic.values[!(names(prog.genetic.values) %in% dead)])
bulmer.effect <- var(prog.genetic.values[Selections]) - var(prog.genetic.values)

result <- sapply(rep(1:length(Selections),1),function(x) {out <-length(which(new.parent.genos[mapinfo$QTLSNP.loci,x,1]=="c")); out2 <- length(which(new.parent.genos[mapinfo$QTLSNP.loci,x,2]=="c"))
outer <- out+out2})
result <- unlist(result)

extraction.info<-list(relmat=rel.mat, num.dead=length(dead), delt.alleles=result, selections=Selections, bulmer.effect=bulmer.effect, select.EBVs = Selection.EBvs,
                      selection.phenos=Selection.phenos,ped=ped, prog.inbred.level=progeny.inbreeding, select.inbred.level=selections.inbreeding,
                      genos.3d=new.parent.genos, num.parents=numselections,select.genval=new.pars.genval, fullped=full.ped,
                      par.ids=select.ids,select.ped.ids=select.ped.ids,all.markers=all.markers,all.phenos=all.phenos, cumulative.total=crossdesign.inuse$cumul.total)
cat("The returned object is a list containing a matrix of phenotypic data with\n")
cat("the specified heritability, a vector of unscaled true genetic values,\n")
cat("to the same total variance as the phenotypic values, and a vector of\n" )
cat("selected individuals, four per family, with the highest phenotype value.\n" )
cat("Both phenotypes and true genetic values include values for parents first followed by progeny.\n")

return(extraction.info)
# # # # #   End of createPhenos() function # # # # # # # #
}

####Create SELFED for each generation####
createSelfProgeny <- function(x,mapinfo=map, total.loci.num=map$total.loci.num, loci=map$loci,
                              parentinfo=parents,phenos=parents.phenos$phenos, TGV=parents.TGV, 
                              chromoloci=map$last.locus.per.chrom,E.var= Env.var,h2=Heritability,...){
  
  crossprog <- 100
  par1.alleles <- parentinfo$genos.3d[,x,]
  par2.alleles <- parentinfo$genos.3d[,x,]
  
  gametes1<-matrix(rep(NA,total.loci.num*crossprog),nrow=total.loci.num,ncol=crossprog,dimnames=list(loci,NULL))
  gametes2<-matrix(rep(NA,total.loci.num*crossprog),nrow=total.loci.num,ncol=crossprog,dimnames=list(loci,NULL))
  
  {
  chr1.r <- sort(sample(1:chromoloci[1],crossprog*1.5,replace=T))
  chr2.r <- sort(sample((chromoloci[1]+1):chromoloci[2],crossprog*1.5,replace=T))
  chr3.r <- sort(sample((chromoloci[2]+1):chromoloci[3],crossprog*1.5,replace=T))
  chr4.r <- sort(sample((chromoloci[3]+1):chromoloci[4],crossprog*1.5,replace=T))
  chr5.r <- sort(sample((chromoloci[4]+1):chromoloci[5],crossprog*1.5,replace=T))
  chr6.r <- sort(sample((chromoloci[5]+1):chromoloci[6],crossprog*1.5,replace=T))
  chr7.r <- sort(sample((chromoloci[6]+1):chromoloci[7],crossprog*1.5,replace=T))
  chr8.r <- sort(sample((chromoloci[7]+1):chromoloci[8],crossprog*1.5,replace=T))
  chr9.r <- sort(sample((chromoloci[8]+1):chromoloci[9],crossprog*1.5,replace=T))
  chr10.r <- sort(sample((chromoloci[9]+1):chromoloci[10],crossprog*1.5,replace=T))
  chr11.r <- sort(sample((chromoloci[10]+1):chromoloci[11],crossprog*1.5,replace=T))
  chr12.r <- sort(sample((chromoloci[11]+1):chromoloci[12],crossprog*1.5,replace=T))
  indiv1 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv2 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv3 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv4 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv5 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv6 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv7 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv8 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv9 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv10 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv11 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv12 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  }
  r <- vector("list",crossprog)
  for (i in 1:(crossprog*1.5)){
    a <- indiv1[i]; b <- indiv2[i]; c <- indiv3[i]; d <- indiv4[i]; e <- indiv5[i]
    f <- indiv6[i]; g <- indiv7[i]; h <- indiv8[i]; i <- indiv9[i]; j <- indiv10[i]
    k <- indiv11[i]; l <- indiv12[i]
    r[[a]] <- unique(sort(c(r[[a]],chr1.r[i]))); r[[b]]<- unique(sort(c(r[[b]],chr2.r[i])))
    r[[c]]<- unique(sort(c(r[[c]],chr3.r[i]))); r[[d]]<- unique(sort(c(r[[d]],chr4.r[i])))
    r[[e]]<- unique(sort(c(r[[e]],chr5.r[i]))); r[[f]]<- unique(sort(c(r[[f]],chr6.r[i])))
    r[[g]]<- unique(sort(c(r[[g]],chr7.r[i]))); r[[h]]<- unique(sort(c(r[[h]],chr8.r[i])))
    r[[i]]<- unique(sort(c(r[[i]],chr9.r[i]))); r[[j]]<- unique(sort(c(r[[j]],chr10.r[i])))
    r[[k]]<- unique(sort(c(r[[k]],chr11.r[i]))); r[[l]]<- unique(sort(c(r[[l]],chr12.r[i])))
  }
  
  for(i in 1:crossprog) {
    allele <- sample(1:2,1)
    end <- chromoloci[12]
    z <-1
    recombination.spots <- r[[i]]
    for (each in 1:length(recombination.spots)) {
      if (each < length(recombination.spots)){
        gametes1[z:recombination.spots[each],i] <- par1.alleles[z:recombination.spots[each],allele]
      } else { gametes1[z:end,i] <- par1.alleles[z:end,allele] }
      if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
      z <- recombination.spots[each] +1
    }}
  {
  chr1.r <- sort(sample(1:chromoloci[1],crossprog*1.5,replace=T))
  chr2.r <- sort(sample((chromoloci[1]+1):chromoloci[2],crossprog*1.5,replace=T))
  chr3.r <- sort(sample((chromoloci[2]+1):chromoloci[3],crossprog*1.5,replace=T))
  chr4.r <- sort(sample((chromoloci[3]+1):chromoloci[4],crossprog*1.5,replace=T))
  chr5.r <- sort(sample((chromoloci[4]+1):chromoloci[5],crossprog*1.5,replace=T))
  chr6.r <- sort(sample((chromoloci[5]+1):chromoloci[6],crossprog*1.5,replace=T))
  chr7.r <- sort(sample((chromoloci[6]+1):chromoloci[7],crossprog*1.5,replace=T))
  chr8.r <- sort(sample((chromoloci[7]+1):chromoloci[8],crossprog*1.5,replace=T))
  chr9.r <- sort(sample((chromoloci[8]+1):chromoloci[9],crossprog*1.5,replace=T))
  chr10.r <- sort(sample((chromoloci[9]+1):chromoloci[10],crossprog*1.5,replace=T))
  chr11.r <- sort(sample((chromoloci[10]+1):chromoloci[11],crossprog*1.5,replace=T))
  chr12.r <- sort(sample((chromoloci[11]+1):chromoloci[12],crossprog*1.5,replace=T))
  indiv1 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv2 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv3 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv4 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv5 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv6 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv7 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv8 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv9 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv10 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv11 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))
  indiv12 <- c(sample(1:crossprog,crossprog,replace=F),sample(1:crossprog,crossprog/2,replace=F))}
  
  r <- vector("list",crossprog)
  for (i in 1:(crossprog*1.5)){
    a <- indiv1[i]; b <- indiv2[i]; c <- indiv3[i]; d <- indiv4[i]; e <- indiv5[i]
    f <- indiv6[i]; g <- indiv7[i]; h <- indiv8[i]; i <- indiv9[i]; j <- indiv10[i]
    k <- indiv11[i]; l <- indiv12[i]
    r[[a]] <- unique(sort(c(r[[a]],chr1.r[i]))); r[[b]]<- unique(sort(c(r[[b]],chr2.r[i])))
    r[[c]]<- unique(sort(c(r[[c]],chr3.r[i]))); r[[d]]<- unique(sort(c(r[[d]],chr4.r[i])))
    r[[e]]<- unique(sort(c(r[[e]],chr5.r[i]))); r[[f]]<- unique(sort(c(r[[f]],chr6.r[i])))
    r[[g]]<- unique(sort(c(r[[g]],chr7.r[i]))); r[[h]]<- unique(sort(c(r[[h]],chr8.r[i])))
    r[[i]]<- unique(sort(c(r[[i]],chr9.r[i]))); r[[j]]<- unique(sort(c(r[[j]],chr10.r[i])))
    r[[k]]<- unique(sort(c(r[[k]],chr11.r[i]))); r[[l]]<- unique(sort(c(r[[l]],chr12.r[i])))}
  
  for(i in 1:crossprog) {
    allele <- sample(1:2,1)
    end <- chromoloci[12]
    z <-1
    recombination.spots <- r[[i]]
    for (each in 1:length(recombination.spots)) {
      if (each < length(recombination.spots)){
        gametes2[z:recombination.spots[each],i] <- par2.alleles[z:recombination.spots[each],allele]
      } else { gametes2[z:end,i] <- par2.alleles[z:end,allele] }
      if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
      x <- recombination.spots[each] +1
    }}
  genos.3d <- abind(gametes1,gametes2,along=3)
  
  
  
  dom.coeff <- 1         # Dominance coeffecient that is set by user
  QTLSNPs <- mapinfo$QTLSNP.loci     # vector of the loci which are snpqtl
  QTLSNPnum <- genos.3d[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
  markers <- mapinfo$available.Markers# a list of all the markers pulled from map object
  set.seed(45234235) # seed is set making the loci which are selected as markers the same every time
  markers.select <- sort(sample(markers,5000,replace=F),decreasing=F)
  num.markers <- length(markers.select) # length of markers that were selected
  marker.select.genos <- genos.3d[markers.select,,] # genotypes of the markers pulled from the current generation
  map.markers <- mapinfo$genetic.map[markers.select,c(1,6)]
  numQTL <- length(mapinfo$rQTL.loci) # the number of additive qtl
  
  
  # Assign each snpqtl a value
  A <- 1    # Major allele is assigned 1  
  a <- -1000  # Minor allele is assigned -20
  
  # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
  numSNPQTL <- mapinfo$total.SNPQTL.num # the number of loci which are snpqtl
  QTLSNPvalues <-matrix(NA,nrow=numSNPQTL,ncol=crossprog) # matrix to hold snpqtl values
  marker.values <- matrix(NA,nrow=num.markers,ncol=crossprog) # matrix to hold marker values
  Capital.genotypes <- vector()
  Lowercase.genotypes <- vector()
  for (i in 1:26){
    Capital.genotypes <- c(Capital.genotypes,paste(LETTERS[i],LETTERS, sep=""))
    Lowercase.genotypes <-  c(Lowercase.genotypes,paste(letters[i],letters, sep=""))
  }
  
  for (i in 1:crossprog){
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
        QTLSNPvalues[QTLSNPcc,i] <- a
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
  rownames(marker.values) <- 1:crossprog
  
  # Convert the 'invisible' rQTL genotypes to numeric matrices, merge the alleles to paired values also
  parQTLallele1 <- matrix(as.integer(genos.3d[mapinfo$rQTL,,1]),nrow=numQTL,ncol=crossprog)
  colnames(parQTLallele1) <- 1:crossprog
  parQTLallele2 <- matrix(as.integer(genos.3d[mapinfo$rQTL,,2]),nrow=numQTL,ncol=crossprog)
  colnames(parQTLallele2) <- 1:crossprog
  QTLvalues     <- matrix(paste(parQTLallele1,parQTLallele2,sep=","),nrow=numQTL,ncol=crossprog)
  dimnames(QTLvalues)<-list(loci[mapinfo$rQTL],1:crossprog)
  
  # Genetic values of progeny
  geneticvals <- colSums(QTLSNPvalues) + colSums(parQTLallele1) + colSums(parQTLallele2)
  
  E <- rnorm(crossprog, 0, (sqrt(E.var*(1-.3))))
  phenos <- round(geneticvals+E,4)
  trueh <- var(geneticvals)/var(phenos)
  print(trueh)
  
  phenos.1 <- as.matrix(phenos); colnames(phenos.1) <- "Vol"; rownames(phenos.1) <- 1:crossprog
  
  gp <- create.gpData(pheno=phenos.1,geno=marker.values); 
  gp <- codeGeno(gp,label.heter = "1")  
  G.Rel.Mat <- kin(gp,ret="realized"); 
  GBLUP <- gpMod(gp, model="BLUP", trait = "Vol", kin=G.Rel.Mat, predict=T, markerEffects=T)
  
  first.in.family <- 1
  Selections <- vector()
  Capital.genotypes <- vector(); Lowercase.genotypes <- vector()
  for (i in 1:26){
    Capital.genotypes <- c(Capital.genotypes,paste(LETTERS[i],LETTERS,LETTERS[i],LETTERS, sep=""))
    Lowercase.genotypes <-  c(Lowercase.genotypes,paste(letters[i],letters, letters[i], letters, sep=""))}
  all.marker.names <- as.numeric(names(GBLUP$markerEffects))
  selfed.progeny.genos <- matrix(paste(genos.3d[,,1],genos.3d[,,2], sep = ""),ncol=crossprog,nrow=60000)
  c <-1
  least.harmful <- vector()
  best.alleles <- vector()
  S <- vector()
  most.imp.per.parent <- vector()
  top.individuals <- vector()
  nae.last <- 0
  for (chr in 1:12){
    first.chr.loci <- mapinfo$last.locus.per.chrom[chr]-5000 + 1
    last.chr.loci <- mapinfo$last.locus.per.chrom[chr]
    chr.marker.num <- length(all.marker.names[all.marker.names >= first.chr.loci  & all.marker.names <= last.chr.loci]) + c -1
    o <- which.min(GBLUP$markerEffects[c:chr.marker.num])
    D.homozygote <- which(selfed.progeny.genos[as.numeric(names(o)),] %in% Capital.genotypes)
    A <- mean(phenos.1[D.homozygote])
    d.homozygote <- which(selfed.progeny.genos[as.numeric(names(o)),] %in% Lowercase.genotypes)
    B <- mean(phenos.1[d.homozygote])
    R <- c(A=A,B=B)
    least.harmful <- c(least.harmful,which.max(R))
    if (names(least.harmful[chr])=="A"){
      best.alleles <- c(best.alleles,unique(selfed.progeny.genos[as.numeric(names(o)),])[which(unique(selfed.progeny.genos[as.numeric(names(o)),]) %in% Capital.genotypes)]) 
    } else {best.alleles <- c(best.alleles,unique(selfed.progeny.genos[as.numeric(names(o)),])[which(unique(selfed.progeny.genos[as.numeric(names(o)),]) %in% Lowercase.genotypes)])}
    S <- c(S,names(o)) 
    c<- chr.marker.num +1
  }
  names(best.alleles) <- S        
  most.imp.per.parent <- sort(GBLUP$markerEffects[which(all.marker.names %in% S & GBLUP$markerEffects < -1)])
  top.individuals <- best.alleles[match(names(most.imp.per.parent),names(best.alleles))]
  
  for (each in 1:length(top.individuals)){
    sst <- strsplit(top.individuals[each], "")[1]
    top.individuals[each] <- paste0(sst[[1]][1],sst[[1]][2])
  }
  
  yo <- list(top.indiv=top.individuals,most.imp=most.imp.per.parent)

  return(yo)}