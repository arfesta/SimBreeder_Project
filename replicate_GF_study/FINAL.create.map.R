####Load Libraries####
library(abind)
library(doSNOW)
library(synbreed)
library(pedigree)

####Create Map#######
createMap <- function(haldane=Haldane,save=Save,numMarkers=NumMarkers,
                      totalQTL=TotalQTL,numSNPQTL=NumSNPQTL,numchromos=num.chromos,
                      maplength=map.length,lowchromosizerange=0.25, 
                      highchromosizerange =.42, totalloci=total.loci,
                      signifdigits=2,rep.num=Rep, MAFS=SNPQTLMAFs) {
  
  # Determine sizes of 12 chromosomes
  # Size is based on the average chromo size multiplied by an upper/lower size range as specified above
  # The last chromosome size is determined by subtracting map length specificed by the summation of previous 11 chromosome sizes
  avgchromo <- maplength/numchromos 
  set.seed(32423) #A seed is set so that chromosome sizes are the exact same every time
  chromosizes <- floor(runif(numchromos, min=(avgchromo - lowchromosizerange*avgchromo), max=(avgchromo + highchromosizerange*avgchromo)))
  chromosizes[12] <- maplength-sum(chromosizes[1:11])
  
  # Determining the number of intervals that will be calculated and the number of loci per chromosome
  numintervals <- totalloci-1 #Specifying the # of intervals that will need to be calcuated
  chromoloci <- rep(totalloci/numchromos,numchromos) #Each chromosome will get an equal number of loci
  all.loci <- vector("list",12)
  for (i in 1:length(chromoloci)){
    if (i==1){
    all.loci[[1]] <- chromoloci[i] } else {
      all.loci[[i]] <- chromoloci[[i]] + all.loci[[i-1]]
    }
  }
  chromo.loci <- unlist(all.loci)
  
  # Creating empty vectors 
  list1<-vector("list", totalloci) # List "chr#" character set for all loci
  list2<-vector("list", totalloci) # List locus# for all chromosomes
  intervals<-rep(NA, totalloci)    # Holds interval distance from each locus to the previous one
  positions<- rep(NA, totalloci)   # Holds the positions of each loci on a chromosome in cM
  locusnames<-rep(NA,totalloci)    # Holds each locus name for all of the loci (chr1_locus#)
  recfreqs<-rep(NA, totalloci)     # Holds Recombination frequency for each of the loci
  
  # Define map function used to create recombination frequencies:
  if(haldane) {mapfun<-function(x){(1-(exp(-2*x)))/2 }
  } else {mapfun <- function(x){ ( (exp(2*x) - exp(-2*x))/(exp(2*x) + exp(-2*x)) )/2}}
  
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
    positions[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]<- round(100*cumsum(intervals[(1+sum(chromoloci[1:k-1])):sum(chromoloci[1:k])]),digits=8)
  }
  
  ## paste together locus names from list2 and list1
  locusnames<- paste(list2,list1,sep="_")
  
  #Create map that has locusnames, distance between loci, empty vector of types (SNPQTL,rQTL, or Marker), MAFs, position on chromosome, & rec freqs
  map<-data.frame(chr=unlist(list2),loci=locusnames, dist= as.numeric(intervals),types= rep(NA,totalloci), MAF=rep(NA,totalloci), pos=as.numeric(positions),recfreqs=as.numeric(recfreqs))
  
  #Now sample from the map to specify SNPQTL & rQTL, the remainder are potential markers that can be used
  all.loci <- 1:totalloci  # vector that contains 1 through the number of all loci
  set.seed(123135)         # seed is set so that same SNPQTL positions are samples each time from all.loci
  SNPQTL <- sample(all.loci,numSNPQTL,replace=FALSE)
  map$types[SNPQTL] <- "snpqtl"   # Specify in the map data frame that these loci are snpqtl
  map$MAF[SNPQTL] <- MAFS         # Assign these loci the specificed minor allele frequencies generated by user
  
  
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
  
  mapinfo<-list(map=map, totalloci=totalloci,totalQTL=totalQTL,numSNPQTL=numSNPQTL,numintervals=numintervals,
                QTLSNPs=sort(SNPQTL),namestem=namestem, rQTL=rQTL, Markers=Markers,all.loci=chromo.loci)
  cat("The returned object is a list containing:\n")
  cat("a matrix of 5 columns with locus names, intervals, types (marker, snp, or qtl), MAFs, and recombination fractions;\n")
  cat("a set of scalars with the numbers of total loci, total QTL, SNPs that are also QTL, total SNPs, and map intervals,\n")
  cat("vectors identifying which loci are SNPs, which are SNP-QTL, and which are 'invisible' QTL,\n")
  cat("and the namestem prefix used in saving text output files.\n")
  return(mapinfo)
  # End of function #
}

####Create Founder Parent Population####
createParents <- function(mapinfo=map,numparents=NumParents,QTLsd=1,save=Save,random=FALSE,rep.num=Rep){
  totalloci <- mapinfo$totalloci  # Specifies the total # of loci by pulling from the map object
  locusnames <- mapinfo$map$loci  # Specifies locus names from map object
  all.MAFs <- mapinfo$map$MAF     # Assign all.MAFs as the minor allele frequencies in the map object
  totalQTL <- mapinfo$totalQTL    # Total number of qtl pulled from map object
  QTLSNPs <- mapinfo$QTLSNPs      # A vector of the loci which are snpqtl
  numSNPQTL <- mapinfo$numSNPQTL  # the number of snpqtl pulled from map object
  numQTL <- totalQTL - numSNPQTL  # the number of additive qtl (rQTL)
  all.Markers <- map$Markers      # A vector of loci which can potentially be markers
  num.markers <- length(all.Markers)  # the number of potential markers
  QTLoci <- sort(mapinfo$rQTL)    # A vector of the loci which are additive QTL
  
  
  #  Sample rQTL values from two normal distributions with mean 0 and std deviation half of that specified by QTLsd variable
  #  Use 'ceiling' function for one distribution, 'floor' function for other, then add them together to get centered distribution
  #  set.seed(123)
  Dist1 <- ceiling(rnorm((numQTL*2), 0, (QTLsd*(sqrt(2)/2))))
  #  set.seed(745)
  Dist2 <- floor(rnorm((numQTL*2), 0, (QTLsd*(sqrt(2)/2))))
  QTLalleles <- Dist1 + Dist2 # Vector contains the values for the rQTL alleles
  
  
  #  Create Dimension names for 2 arrays
  alleleIDs<-c("a1","a2")
  parIDs<-1:numparents
  
  # Parents array holds both alleles for all parents
  # QTLSNP array holds the alleles for the loci of all parents under dominance control  
  parents<-array(rep(NA,totalloci*numparents*2), dim=c(totalloci,numparents,2),dimnames=list(locusnames,parIDs,alleleIDs))
  QTLSNParray<-array(0, dim=c(numSNPQTL,numparents,2),dimnames=c(QTLSNPs,parIDs,alleleIDs))
  
  #  Create empty vectors to hold either allele1 or allele2  for each parent
  allele1<-rep(NA, totalloci)
  allele2<-rep(NA, totalloci)
  
  # For each parent assign them alleles 
  if (numparents <= 26){
    # If there are less than 26 parents we can assign markers a different letter of the alphabet 
    #    so that we can keep track of where parental haplotypes come from
    for(par in 1:numparents){
      # First we assign parent alleles for all the snpqtl
      rand1<-runif(numSNPQTL, min=0, max=1) # Generate a random set of numbers to test against MAFs to determine first allele
      rand2<-runif(numSNPQTL, min=0, max=1) # Generate a random set of numbers to test against MAFs to determine second allele
      
      # If the MAF of the snpqtl is greater than the random number generated it is a minor allele, otherwise it is the major allele
      minor1 <- QTLSNPs[which(all.MAFs[QTLSNPs] > rand1)] 
      major1 <- QTLSNPs[which(!QTLSNPs %in% minor1)]
      minor2 <- QTLSNPs[which(all.MAFs[QTLSNPs] > rand2)]
      major2 <- QTLSNPs[which(!QTLSNPs %in% minor2)]
      allele1[minor1]<-"c"            # minor allele is always "c", major allele is always "a"
      allele1[major1] <-"a"
      allele2[minor2]<-"c"
      allele2[major2]<-"a"
      
      # Next we repeat the process above for the marker loci
      rand1<-runif(num.markers, min=0, max=1) 
      rand2<-runif(num.markers, min=0, max=1)
      minor1<-all.Markers[which(all.MAFs[all.Markers] > rand1)]
      major1 <- all.Markers[which(!all.Markers %in% minor1)]
      minor2<-all.Markers[which(all.MAFs[all.Markers] > rand2)]
      major2 <- all.Markers[which(!all.Markers %in% minor2)]
      
      #allele1[minor1] <- letters[a]
      #allele1[major1] <- LETTERS[a]
      #allele2[minor2] <- letters[a]
      #allele2[major2] <- LETTERS[a]
      
      allele1[minor1] <- "d"
      allele1[major1] <- "D"
      allele2[minor2] <- "d"
      allele2[major2] <- "D"
      
      # Now we can assign the marker and snpqtl alleles to the parent 
      parents[,par,1]<-allele1
      parents[,par,2]<-allele2
      
      # Finally we sample from the rQTL distribution that was created above to assign each rqtl a numeric value instead of character
      parents[QTLoci,par,1] <- sample(QTLalleles,numQTL)
      parents[QTLoci,par,2] <- sample(QTLalleles,numQTL)
      
      # Write a copy of the SNP alleles at the snpqtl into an array to use in calculating parental values
      QTLSNParray[,par,1]<-allele1[QTLSNPs]
      QTLSNParray[,par,2]<-allele2[QTLSNPs]
    }} else {for(par in 1:numparents){
      rand1<-runif(numSNPQTL, min=0, max=1) 
      rand2<-runif(numSNPQTL, min=0, max=1) 
      minor1 <- QTLSNPs[which(all.MAFs[QTLSNPs] > rand1)]
      major1 <- QTLSNPs[which(!QTLSNPs %in% minor1)]
      minor2 <- QTLSNPs[which(all.MAFs[QTLSNPs] > rand2)]
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
      
      allele1[minor1] <- "d"
      allele1[major1] <- "D"
      allele2[minor2] <- "d"
      allele2[major2] <- "D"
      
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
  
  # If save is true then then we save the parentSNPgenos object and parent.markers object
  if(save) {
    parfilename=paste(rep.num,"parent-SNPQTL-genotypes.txt")
    write.table(parentSNPgenos,parfilename, quote = F, row.names = T, col.names = T, sep= "\t")
    
    parfilename=paste(rep.num,"parent-marker-genotypes.txt")
    write.table(parent.markers,parfilename, quote = F, row.names = T, col.names = T, sep= "\t")    
  }
  
  ###IMP that genos.3d is first list element and numparents is second becuase that is what createTGV routine expects###
  parentinfo<-list(genos.3d=parents,numparents=numparents, parSNPs=parentSNPgenos, par.markers=parent.markers, par.IDs=parIDs)
  cat("The returned object is a list containing a 3-D array of parent genotypes\n")
  cat("a matrix of parent SNP genotypes with individuals in columns, a scalar of the number of parents,\n")
  cat("matrices of the parental genotypes at all SNPs, parental QTL allele values, genotypes and values at SNP QTL,\n")
  cat("two matrices with variances of QTL alleles & genetic values and SNP QTL alleles & genetic values, and a list \n")
  cat("of the numerical values assigned to major and minor alleles at SNP-QTL in the parents.\n")
  return(parentinfo)
  # # # # #    End of function    # # # # # # # #
}

####Create Total Genetic Value####
createTGV <- function(mapinfo=map, data1=parents, save=Save,
                      rep.num=Rep, prefix=prefix, Dom.coeff=dom.coeff, 
                      num.markers=NumMarkers){
  
  dom.coeff <- Dom.coeff         # Dominance coeffecient that is set by user
  locusnames <- mapinfo$map$loci # The locus names pulled from the mab object
  QTLSNPs <- mapinfo$QTLSNPs     # vector of the loci which are snpqtl
  QTLSNPnum <- data1$genos.3d[QTLSNPs,,] # genotypes of both alleles pulled from the current generation
  markers <- mapinfo$Markers # a list of all the markers pulled from map object
  set.seed(45234235) # seed is set making the loci which are selected as markers the same every time
  markers.select <- sort(sample(markers,num.markers,replace=F),decreasing=F)
  num.markers <- length(markers.select) # length of markers that were selected
  marker.select.genos <- data1$genos.3d[markers.select,,] # genotypes of the markers pulled from the current generation
  map.markers <- mapinfo$map[markers.select,c(1,6)]
  numQTL <- length(mapinfo$rQTL) # the number of additive qtl
  parIDs <- data1$par.IDs # the parent ids pulled from the current generation
  
  
  # Assign each snpqtl a value
  A <- 1    # Major allele is assigned 1  
  a <- -.1 # Minor allele is assigned -20
  
  # Create 2 matrices: One to hold snpqtl values for all parents and the other to hold marker marker values for blup analysis
  numSNPQTL <- mapinfo$numSNPQTL # the number of loci which are snpqtl
  numparents <- length(parIDs) # the number of parents
  QTLSNPvalues <-matrix(NA,nrow=numSNPQTL,ncol=numparents) # matrix to hold snpqtl values
  marker.values <- matrix(NA,nrow=num.markers,ncol=numparents) # matrix to hold marker values
  
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
        QTLSNPvalues[QTLSNPca,i] <- (A*dom.coeff)
      }
    markers.aa <- which(marker.select.genos[,i,1]=="D" & marker.select.genos[,i,2]=="D")
    markers.cc <- which(marker.select.genos[,i,1]=="d" & marker.select.genos[,i,2]=="d")
    markers.ac <- which(marker.select.genos[,i,1]=="D" & marker.select.genos[,i,2]=="d")
    markers.ca <- which(marker.select.genos[,i,1]=="d"  & marker.select.genos[,i,2]=="D")
    
    marker.values[markers.aa,i] <- "0"
    marker.values[markers.cc,i] <- "2"
    marker.values[markers.ac,i] <- "1"
    marker.values[markers.ca,i] <- "1"
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
  
  TGV <- list(pargeneticvals=geneticvals, ParSNPval=QTLSNPvalues, gblup.val=marker.values,marker.names=markers.select, marker.map=map.markers)
  return(TGV)
}

####Create Unscaled Parent Phenos####
createUnscaledParentPhenos <- function(ParentTGV=parents.TGV, h2=Heritability, rep.num=Rep) {
 
  prefix="parent.phenos"
  pargeneticvals<-ParentTGV$pargeneticvals
  totalindiv<-length(pargeneticvals)
  phenos<-vector(length=totalindiv)
  trueh<-vector(length=totalindiv)
  
  for(i in 1:5000){
  set.seed(2342342)
  E <- rnorm(totalindiv, 0, (sqrt((i)*(1-h2))))# for U-shape MAF and 0 dominance 192 parents
  P <- round(pargeneticvals+E,4)
  trueh <- round(var(pargeneticvals)/var(P),2)
  print(trueh)
  E.var <- i
  if (trueh==.3) break
  }
  
  pheno.info<-list(phenos=P, pargeneticvals=pargeneticvals, E.var = E.var, E=E)
  
  return(pheno.info)
}

####Create progeny for each generation####
createProgeny <- function(mapinfo=map,parentinfo=parents, cumul.total=progeny1$cumul.total, TGVdata=parents.TGV, 
                          crossfile=CrossFile1, gen=1, savegenos=Save,saveped=Save,prefix=save.prefix,rep.num=Rep,sbpm=SBPM, gblup=GBLUP, ablup=ABLUP,
                          selfl=SELFL,dp=DP,selfp=SELFP, sbpw=SBPW, halfdiallel=HD){
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Read the crossdesign matrix to determine the number of crosses and total number of progeny,
  # and read the identities of parent1, parent2, and number of progeny desired, from those vectors
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  NumParents=parentinfo$numparents
  chromoloci <- mapinfo$all.loci
  if (sbpm | sbpw | gblup |ablup){
    set.seed(gen+23523)
    crossdesign= matrix(data = sample(x=1:NumParents,replace=FALSE, size=NumParents),nrow = NumParents/2,ncol=2)
    crossdesign <- cbind(crossdesign,rep(100,length(crossdesign[,1])))
  }else if (selfl | selfp){
    crossdesign=matrix(data=c(rep(1:NumParents,times=2),rep(NA,NumParents)),nrow=NumParents,ncol=3)
    crossdesign[,3] = rep(50)
  } else if (dp) {
    if(gen==1){
      colClasses = c("integer", "integer", "integer")
      crossdesign = as.matrix(read.table(crossfile,header=F,sep="\t", colClasses))} else {
        crossdesign=matrix(data=c(rep(1:NumParents),sample(1:NumParents,replace=FALSE,size=NumParents),rep(50,NumParents)),nrow=NumParents,ncol = 3)
      }
  } else if(halfdiallel){
    par.1 <- vector()
    par.2 <- vector()
    c <- NumParents
    for (i in 1:NumParents){
      par1 <- rep(i,c)
      if (i < NumParents){
        par2 <- c(rep((i+1):NumParents),i) } else{
          par2 <- i
        }
      par.1 <- c(par.1,par1)
      par.2 <- c(par.2,par2)
      c <- c-1
    }
    progeny <- rep(6, length(par.1))
    crossdesign=matrix(data=c(par.1,par.2,progeny),nrow=length(progeny),ncol=3)
    
  } else {
    colClasses = c("integer", "integer", "integer")
    crossdesign = as.matrix(read.table(crossfile,header=F,sep="\t", colClasses))
  } 
  
  numcrosses   <- length(crossdesign[,1])
  totalprogeny <- sum(as.numeric(crossdesign[,3]))
  totalindiv   <- totalprogeny+parentinfo$numparents
  QTLSNPalleles<- TGVdata$parQTLSNPalleles
  QTLSNPs      <- mapinfo$QTLSNPs
  
  
  for(m in 1:numcrosses){
    par1<-crossdesign[m,1] # assigns par1 to be the first parent in crossdesign matrix
    par2<-crossdesign[m,2] # assigns par2 to be the second parent in the crossdesign matrix
    crossprog<-as.numeric(crossdesign[m,3]) #assigns number of progeny to be the third column for cross "X"
    
    # Create empty matrix to hold gametes
    # dimensions are (total # of loci) x  (# of cross progeny)
    # rownames are the loci names
    gametes1<-matrix(rep(NA,mapinfo$totalloci*crossprog),nrow=mapinfo$totalloci,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
    gametes2<-matrix(rep(NA,mapinfo$totalloci*crossprog),nrow=mapinfo$totalloci,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
    
    
    r <- vector("list",crossprog)
    par1.alleles <- (parentinfo$genos.3d[,par1,])
    par2.alleles <- (parentinfo$genos.3d[,par2,])
    
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
    for (i in 1:(crossprog*1.5)){
      a <- indiv1[i]
      b <- indiv2[i]
      c <- indiv3[i]
      d <- indiv4[i]
      e <- indiv5[i]
      f <- indiv6[i]
      g <- indiv7[i]
      h <- indiv8[i]
      i <- indiv9[i]
      j <- indiv10[i]
      k <- indiv11[i]
      l <- indiv12[i]
      r[[a]] <- unique(sort(c(r[[a]],chr1.r[i])))
      r[[b]]<- unique(sort(c(r[[b]],chr2.r[i])))
      r[[c]]<- unique(sort(c(r[[c]],chr3.r[i])))
      r[[d]]<- unique(sort(c(r[[d]],chr4.r[i])))
      r[[e]]<- unique(sort(c(r[[e]],chr5.r[i])))
      r[[f]]<- unique(sort(c(r[[f]],chr6.r[i])))
      r[[g]]<- unique(sort(c(r[[g]],chr7.r[i])))
      r[[h]]<- unique(sort(c(r[[h]],chr8.r[i])))
      r[[i]]<- unique(sort(c(r[[i]],chr9.r[i])))
      r[[j]]<- unique(sort(c(r[[j]],chr10.r[i])))
      r[[k]]<- unique(sort(c(r[[k]],chr11.r[i])))
      r[[l]]<- unique(sort(c(r[[l]],chr12.r[i])))
    }
    
    for(i in 1:crossprog) {
      allele <- sample(1:2,1)
      end <- chromoloci[12]
      x <-1
      recombination.spots <- r[[i]]
      for (each in 1:length(recombination.spots)) {
        if (each < length(recombination.spots)){
          gametes1[x:recombination.spots[each],i] <- par1.alleles[x:recombination.spots[each],allele]
        } else { gametes1[x:end,i] <- par1.alleles[x:end,allele] }
        if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
        x <- recombination.spots[each] +1
      }}
    
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
    
    r <- vector("list",crossprog)
    for (i in 1:(crossprog*1.5)){
      a <- indiv1[i]
      b <- indiv2[i]
      c <- indiv3[i]
      d <- indiv4[i]
      e <- indiv5[i]
      f <- indiv6[i]
      g <- indiv7[i]
      h <- indiv8[i]
      i <- indiv9[i]
      j <- indiv10[i]
      k <- indiv11[i]
      l <- indiv12[i]
      r[[a]] <- unique(sort(c(r[[a]],chr1.r[i])))
      r[[b]]<- unique(sort(c(r[[b]],chr2.r[i])))
      r[[c]]<- unique(sort(c(r[[c]],chr3.r[i])))
      r[[d]]<- unique(sort(c(r[[d]],chr4.r[i])))
      r[[e]]<- unique(sort(c(r[[e]],chr5.r[i])))
      r[[f]]<- unique(sort(c(r[[f]],chr6.r[i])))
      r[[g]]<- unique(sort(c(r[[g]],chr7.r[i])))
      r[[h]]<- unique(sort(c(r[[h]],chr8.r[i])))
      r[[i]]<- unique(sort(c(r[[i]],chr9.r[i])))
      r[[j]]<- unique(sort(c(r[[j]],chr10.r[i])))
      r[[k]]<- unique(sort(c(r[[k]],chr11.r[i])))
      r[[l]]<- unique(sort(c(r[[l]],chr12.r[i])))
    }
    
    for(i in 1:crossprog) {
      allele <- sample(1:2,1)
      end <- chromoloci[12]
      x <-1
      recombination.spots <- r[[i]]
      for (each in 1:length(recombination.spots)) {
        if (each < length(recombination.spots)){
          gametes2[x:recombination.spots[each],i] <- par2.alleles[x:recombination.spots[each],allele]
        } else { gametes2[x:end,i] <- par2.alleles[x:end,allele] }
        if (allele==1) { allele <- allele +1 } else { allele <- allele -1}
        x <- recombination.spots[each] +1
      }}
    
    if(m==1){
      output<-matrix(paste(gametes1,gametes2,sep=","),nrow=mapinfo$totalloci,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL))
      array.out <- abind(gametes1,gametes2,along=3)
    } else {
      output<-cbind(output,matrix(paste(gametes1,gametes2,sep=","),nrow=mapinfo$totalloci,ncol=crossprog,dimnames=list(mapinfo$map$loci,NULL)))
      array.out <- abind(array.out,abind(gametes1,gametes2,along=3),along=2)
    }
  }
  
  # end of m in 1:numcrosses loop, start a new cross and reinitialize parent IDs & randnum matrices
  colnames(output) <- 1:totalprogeny
  
  # Write diploid genotypes in output matrix to a text file as backup. Matrix and text file have loci as rows, individuals as columns
  # NB: output matrix has SNPs in rows, individuals in columns
  if(savegenos){
    outtablename<-paste(rep.num,prefix,"-genotypes.txt",sep="")
    write.table(output,outtablename,quote=F,row.names=T,col.names=T, sep="\t")
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #     Create a pedigree file from the mating design information in the crossdesign matrix         #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  par1list<-"par1"
  par2list<-"par2"
  par1id<-"par1"
  par2id<-"par2"
  generation<- c("generation",rep(gen,totalprogeny))
  
  
  if (gen==1){
    totalindiv   <- totalprogeny+length(unique(crossdesign[,1:2]))
    cumul.total <- length(unique(crossdesign[,1:2])) + totalprogeny
    indivIDs<-c("indivIDs",seq((length(unique(crossdesign[,1:2]))+1),(totalindiv)))
    for (m in 1:numcrosses){
      crossprog<-as.numeric(crossdesign[m,3])
      par1<-as.character(crossdesign[m,1])
      par2<-as.character(crossdesign[m,2])
      par1list<-c(par1list,rep(par1,crossprog))
      par2list<-c(par2list,rep(par2,crossprog))}
    
    parentped <- cbind("ID"=rep(1:length(unique(crossdesign[,1:2])),1),"Par1" = rep(0,length(unique(crossdesign[,1:2]))), "Par2"= rep(0,length(unique(crossdesign[,1:2]))),"gener"=rep(0,length(unique(crossdesign[,1:2]))))
    progped<-cbind("ID"=indivIDs[2:length(indivIDs)], "Par1"=par1list[2:length(par1list)],"Par2"=par2list[2:length(par1list)],"gener"=generation[2:length(generation)])
    full.ped <- rbind(parentped,progped)
  } else {
    indivIDs<- c("indivIDs",seq(cumul.total+1,cumul.total+totalprogeny))
    cumul.total<- cumul.total + totalprogeny
    vector2 <- parentinfo$select.ped.ids
    for (m in 1:numcrosses){
      crossprog<-as.numeric(crossdesign[m,3])
      par1<-c(crossdesign[m,1])
      par2<-c(crossdesign[m,2])
      par1id<- (vector2[par1])
      par2id<- (vector2[par2])
      par1list<-c(par1list,rep(par1id,crossprog))
      par2list<-c(par2list,rep(par2id,crossprog))}
    progped<-cbind("ID"=indivIDs[2:length(indivIDs)], "Par1"=par1list[2:length(par1list)],"Par2"=par2list[2:length(par1list)],"gener"=generation[2:length(generation)])
    full.ped <- rbind(parentinfo$ped,progped)
  }
  
  
  if(saveped) {
    pedfilename=paste(rep.num,mapinfo$namestem,prefix,"-pedigree.txt",sep="")
    write.table(progped,pedfilename, quote = F, row.names = F, col.names = T, sep=" ")
  }
  ###Imp that genos.3d is first in list item, and total progeny is next list item bc in createtgv routine it expects genos.3d as first element and total prog as second
  progeny<-list(genos.3d=array.out, totalprogeny=totalprogeny, pedigree=progped, fullped= full.ped, 
                crossdesign=crossdesign,progeny.all=output,numparents=parentinfo$numparents, prog.ped=progped,
                numcrosses=numcrosses,crossprog=crossprog,namestem=mapinfo$namestem,cumul.total=cumul.total, par.IDs=progped[,1])
  cat("The returned object is a list containing a matrix of progeny SNP\n")
  cat("genotypes with individuals in rows and SNPs in columns,\n")
  cat("a matrix with the full pedigree (parents first, followed by progeny),\n")
  cat("a matrix with the mating design from the cross design file),\n")
  cat("a vector of progeny true genetic values,\n")
  cat("a matrix with progeny genotypes at both SNP and QTL loci,\n")
  cat("with loci in rows and individuals in columns, and scalars with the \n")
  cat("numbers of parents & crosses in the mating design, the number of\n")
  cat("progeny per cross, the total number of progeny, and the namestem.\n")
  return(progeny)
  # # # # # #   End of function # # # # # # #
}

####Create Unscaled Progeny Phenos & Make Selections####
createUnscaledPhenos <- function(progenyinfo=progeny1, h2=Heritability,past.phenos=NULL,
                                 prefix=prefix,rep.num=Rep,progenyTGV=progeny1.TGV,sbpm=SBPM,selfl=SELFL, dp=DP, E.var = parents.phenos$E.var,
                                 selfp=SELFP, sbpw=SBPW, halfdiallel=HD,mapinfo=map,parentphenos=parents.phenos,gen=1, 
                                 selftest=SELFTEST, gBLUP=GBLUP,aBLUP=ABLUP, parents.tgv=parents.TGV, parentinfo=parents,
                                 numSelections=NumSelections) {
  
  crossdesign <- progenyinfo$crossdesign
  pedigree <- progenyinfo$prog.ped
  numparents <- progenyinfo$numparents
  geneticvalues<-progenyTGV$pargeneticvals
  totalindiv<-length(geneticvalues)
  phenos<-vector(length=totalindiv)
  trueh<-vector(length=totalindiv)
  MAFdist = mapinfo$MAFdist
  prog.markers <- progenyTGV$gblup.val
  map.markers <- progenyTGV$marker.map
  colnames(map.markers) <- c("chr","pos")
  map.markers$chr <- as.character(map.markers$chr)
  
  parent.markers <- parents.tgv$gblup.val

  #####Create Unsclaed Phenotypes
  set.seed(2342342)
  E <- rnorm(totalindiv, 0, (sqrt(E.var*(1-h2))))# for U-shape MAF and 0 dominance 192 parents
  phenos <- round(geneticvalues+E,4)
  trueh <- var(geneticvalues)/var(phenos)
  print(trueh)
  
  ####Generate new pedigree 
  if (gen==1) {
    ancestor <- data.frame(ID=1:numparents,Par1=rep(NA,numparents),Par2=rep(NA,numparents),gener=rep(gen-1,numparents))
    pedigr <- rbind(ancestor,pedigree)
    ped <- create.pedigree(ID=pedigr[,1], pedigr[,2],pedigr[,3]) 
    all.markers <- rbind(parent.markers,prog.markers)
    all.phenos <-c(parentphenos$phenos,phenos)
  } else {
    ped <- progenyinfo$fullped
    all.markers <- rbind(past.phenos$all.markers,prog.markers)
    all.phenos <- c(past.phenos$all.phenos,phenos)
  }
  #ped <- create.pedigree(ID=pedigr[,1], pedigr[,2],pedigr[,3])}
  ped <- ped[order(as.numeric(ped$ID)),] 
  
  
  NUMSelections <- c(1:numSelections)
  first.in.family <- 1
  Selections <- vector()
  
  if(aBLUP) {
    phenos.1 <- as.matrix(all.phenos)
    colnames(phenos.1) <- "Vol"
    gp <- create.gpData(pheno=phenos.1,geno=all.markers, pedigree=ped)
    gp <- codeGeno(gp,label.heter = "1")
    U <- kin(gp,ret="kin")
    a <- gpMod(gp, model="BLUP", trait = "Vol", kin=U, markerEffects=F,predict=T)
    g <- match(names(phenos),rownames(a$fit$predicted))
    yo <- a$g[g]
    names(yo) <- names(phenos)
    
    for(family in 1:length(crossdesign[,3])){
      num.offspring <- crossdesign[family,3]
      last.in.family <- num.offspring + first.in.family - 1 
      temp <- (yo[first.in.family:last.in.family])
      sorted <-sort(temp,decreasing=TRUE)
      
      BestOne <- which(temp==sorted[1])
      SecondOne <- which(temp==sorted[2])
      selected1 <- BestOne + first.in.family - 1
      selected2 <- SecondOne + first.in.family - 1
      Selections <- c(Selections, selected1, selected2)
      first.in.family <- last.in.family + 1
    }
    Selection.phenos<- phenos[Selections]
    sorted.top192 <- sort(Selection.phenos,decreasing=T)
    PP.phenos <- sorted.top192[seq(1,192,8)]
    PP.24 <- which(phenos %in% PP.phenos)
    
    new.parent.genos <- progenyinfo$genos.3d[,Selections,]
    numselections <- dim(new.parent.genos)[2]
    new.pars.genval <- geneticvalues[Selections]
    pp.genval <- geneticvalues[PP.24]
    select.ids <- as.numeric(names(new.pars.genval))-numparents
    select.ped.ids <- as.numeric(names(new.pars.genval)) 
    
  } else if (gBLUP) {
    phenos.1 <- as.matrix(all.phenos)
    colnames(phenos.1) <- "Vol"
    gp <- create.gpData(pheno=phenos.1,geno=all.markers, pedigree=ped)
    gp <- codeGeno(gp,label.heter = "1")
    U <- kin(gp,ret="realized")
    a <- gpMod(gp, model="BLUP", trait = "Vol", kin=U, predict=T, markerEffects=F)
    g <- match(names(phenos),rownames(a$fit$predicted))
    yo <- a$fit$predicted[g]
    names(yo) <- names(phenos)
    
    for(family in 1:length(crossdesign[,3])){
      num.offspring <- crossdesign[family,3]
      last.in.family <- num.offspring + first.in.family - 1 
      temp <- (yo[first.in.family:last.in.family])
      sorted <-sort(temp,decreasing=TRUE)
      
      BestOne <- which(temp==sorted[1])
      SecondOne <- which(temp==sorted[2])
      selected1 <- BestOne + first.in.family - 1
      selected2 <- SecondOne + first.in.family - 1
      Selections <- c(Selections, selected1, selected2)
      first.in.family <- last.in.family + 1
    }
    Selection.phenos<- phenos[Selections]
    sorted.top192 <- sort(Selection.phenos,decreasing=T)
    PP.phenos <- sorted.top192[seq(1,192,8)]
    PP.24 <- which(phenos %in% PP.phenos)
    
    new.parent.genos <- progenyinfo$genos.3d[,Selections,]
    numselections <- dim(new.parent.genos)[2]
    new.pars.genval <- geneticvalues[Selections]
    pp.genval <- geneticvalues[PP.24]
    select.ids <- as.numeric(names(new.pars.genval))-numparents
    select.ped.ids <- as.numeric(names(new.pars.genval)) 
  } else if (selftest) {
    q <- which(crossdesign[,1]==crossdesign[,2])
    lov <- matrix(paste(progenyinfo$genos.3d[,,1],progenyinfo$genos.3d[,,2], sep = ""),ncol=totalindiv,nrow=700)
    p <-1
    a <-1
    lm <-1
    SELECTIONS <- vector()
    ##For each selfing cross, first idnentify the marker effects estimated from a GLBUP
    for (each in 1:length(q)){
      i <- q[each]  #row 15
      first.prog <- sum(crossdesign[1:(i-1),3]) +1 #first.prog is the sum of all previous progeny +1
      num.progeny <- crossdesign[i,3] # number of progeny in the self
      last.prog <- first.prog+num.progeny -1 #last progeny of the selfs
      phenos.2 <- phenos[first.prog:last.prog] # phenos for the self
      phenos.2 <- as.matrix(phenos.2)
      colnames(phenos.2) <- "Vol"
      prog.markers.2 <- prog.markers[first.prog:last.prog,] # markers for the self
      gp <- create.gpData(pheno=phenos.2,geno=prog.markers.2)
      gp <- codeGeno(gp,label.heter = "1")
      U <- kin(gp,ret="realized")
      a <- gpMod(gp, model="BLUP", trait = "Vol", kin=U, markerEffects=T) #returns marker effects for the self
      
      ####Marker effects are now estimated for the selfed progeny
      ###We now estimate which marker for each linkage group has the largest negative effect
      c <-1
      least <- vector()
      S <- vector()
      for (p in 1:12){
        last <- c+map$marker.chrom[p]-1
        o <- which(a$markerEffects[c:last]==(min(a$markerEffects[c:last])))
        RE <- which(progenyinfo$genos.3d[as.numeric(names(o)[1]),first.prog:last.prog,1]=="A" & progenyinfo$genos.3d[as.numeric(names(o)[1]),first.prog:last.prog,2]=="A")
        A <- mean(phenos.2[RE])
        RE.1 <- which(progenyinfo$genos.3d[as.numeric(names(o)[1]),first.prog:last.prog,1]=="a" & progenyinfo$genos.3d[as.numeric(names(o)[1]),first.prog:last.prog,2]=="a")
        B <- mean(phenos.2[RE.3])
        R <- c(A=A,B=B)
        least <- c(least,which.max(R))
        S <- c(S,names(o)[1])
        if (p ==12) break
        c <- map$marker.chrom[p]+c
        p <- p +1
      }
      ###We have identified the allele combination for the marker which is associated with each linakge group to be the most informative
      ##We want to sort those of the marker effects by the 12 most impactful markers which explain each linkage group
      worst <- sort(a$markerEffects[S])
      f <- which(least=="1")
      s <- which(least=="2")
      top.par <- matrix(,nrow=12,ncol=1)
      top.par[f,] <- c("AA")
      top.par[s,] <- c("aa")
      rownames(top.par) <- S
      z <- match(rownames(top.par),names(worst))
      top.par <- top.par[z,]
      #####top.par conatains the most iomportant to
      
      first.in.family <-1
      selections <- vector()
      for(family in lm:i){
        num.offspring <- crossdesign[family,3]
        last.in.family <- num.offspring + first.in.family - 1 
        cs <-1
        tab <- vector()
        for (cs in 1:12){
          a <- which(lov[as.numeric(names(top.par))[cs],first.in.family:last.in.family] == top.par[cs])
          tab <- c(tab,a)
        }
        m <- sort(table(tab),decreasing=T)
        if(family==1){
          y <- as.numeric(names(m))[1:2]}
        else {
          y <- as.numeric(names(m))[1:2] +first.in.family -1 }
        selections <- c(selections,y)
        first.in.family <- last.in.family + 1
      }
      if (each > 1){
        selections <- selections +last.prog - last.in.family
      }
      SELECTIONS <- c(SELECTIONS,selections)
      lm<-i +1
    }
    Selection.phenos<- phenos[SELECTIONS]
    sorted.top192 <- sort(Selection.phenos,decreasing=T)
    PP.phenos <- sorted.top192[seq(1,192,8)]
    PP.24 <- which(phenos %in% PP.phenos)
    
    new.parent.genos <- progenyinfo$genos.3d[,SELECTIONS,]
    numselections <- dim(new.parent.genos)[2]
    new.pars.genval <- geneticvalues[SELECTIONS]
    pp.genval <- geneticvalues[PP.24]
    select.ids <- as.numeric(names(new.pars.genval))-192
    select.ped.ids <- as.numeric(names(new.pars.genval))
    
  }else if(sbpm | selfp){
    sorted.phenos <- sort(phenos,decreasing=TRUE) #all phenos sorted
    Selections <- which(phenos %in% sorted.phenos[1:192]) # id of top 192 sorted phenos
    PP.24 <- which(phenos %in% sorted.phenos[seq(1,192,8)]) # id of equally spread out 24 out of top 192
    
    PP.phenos <- phenos[PP.24] #phenos of production population
    Selection.phenos <- phenos[Selections] #phenos of breeding population
    
    new.parent.genos <- progenyinfo$genos.3d[,Selections,]
    numselections <- dim(new.parent.genos)[2]
    new.pars.genval <- geneticvalues[Selections]
    pp.genval <- geneticvalues[PP.24]
    #selection id within progeny test
    select.ids <- as.numeric(names(new.pars.genval))-192
    #selection id within pedigree
    select.ped.ids <- as.numeric(names(new.pars.genval)) 
    
  } else if(sbpw | halfdiallel){
    #Selections based on unscaled phenotypes
    for(family in 1:length(crossdesign[,3])){
      num.offspring <- crossdesign[family,3]
      last.in.family <- num.offspring + first.in.family - 1 
      temp <- (phenos[first.in.family:last.in.family])
      sorted <-sort(temp,decreasing=TRUE)
      for (i in 1:numSelections){
        BestOne <- which(temp==sorted[i])
        #SecondOne <- which(temp==sorted[2])
        selected1 <- BestOne + first.in.family - 1
        #selected2 <- SecondOne + first.in.family - 1
        Selections <- c(Selections, selected1)}
      first.in.family <- last.in.family + 1
    }
    Selection.phenos<- phenos[Selections]
    sorted.top192 <- sort(Selection.phenos,decreasing=T)
    PP.phenos <- sorted.top192[seq(1,192,8)]
    PP.24 <- which(phenos %in% PP.phenos)
    
    new.parent.genos <- progenyinfo$genos.3d[,Selections,]
    numselections <- dim(new.parent.genos)[2]
    new.pars.genval <- geneticvalues[Selections]
    pp.genval <- geneticvalues[PP.24]
    select.ids <- as.numeric(names(new.pars.genval))-192
    select.ped.ids <- as.numeric(names(new.pars.genval)) 
    
  } else if (selfl){
    for(family in 1:length(crossdesign[,3])){
      num.offspring <- crossdesign[family,3]
      last.in.family <- num.offspring + first.in.family - 1 
      temp <- (phenos[first.in.family:last.in.family])
      sorted <-sort(temp,decreasing=TRUE)
      
      BestOne <- which(temp==sorted[1])
      SecondOne <- which(temp==sorted[2])
      selected1 <- BestOne + first.in.family - 1
      Selections <- c(Selections, selected1)
      first.in.family <- last.in.family + 1
    }
    Selection.phenos<- phenos[Selections]
    sorted.top192 <- sort(Selection.phenos,decreasing=T)
    PP.phenos <- sorted.top192[seq(1,192,8)]
    PP.24 <- which(phenos %in% PP.phenos)
    
    new.parent.genos <- progenyinfo$genos.3d[,Selections,]
    numselections <- dim(new.parent.genos)[2]
    new.pars.genval <- geneticvalues[Selections]
    pp.genval <- geneticvalues[PP.24]
    select.ids <- as.numeric(names(new.pars.genval))-192
    select.ped.ids <- as.numeric(names(new.pars.genval))
    
  } else if(dp){
    for(family in 1:length(crossdesign[,3])){
      num.offspring <- crossdesign[family,3]
      last.in.family <- num.offspring + first.in.family - 1 
      temp <- (phenos[first.in.family:last.in.family])
      sorted <-sort(temp,decreasing=TRUE)
      
      BestOne <- which(temp==sorted[1])
      selected1 <- BestOne + first.in.family - 1
      Selections <- c(Selections, selected1)
      first.in.family <- last.in.family + 1
    }
    Selection.phenos<- phenos[Selections]
    sorted.top192 <- sort(Selection.phenos,decreasing=T)
    PP.phenos <- sorted.top192[seq(1,20,1)]
    PP.24 <- which(phenos %in% PP.phenos)
    
    new.parent.genos <- progenyinfo$genos.3d[,Selections,]
    numselections <- dim(new.parent.genos)[2]
    new.pars.genval <- geneticvalues[Selections]
    pp.genval <- geneticvalues[PP.24]
    select.ids <- as.numeric(names(new.pars.genval))-#00
      select.ped.ids <- as.numeric(names(new.pars.genval)) 
  }
  
  prog.inbreed <- calcInbreeding(ped[,1:3])
  prog.inbreeding <- prog.inbreed[as.numeric(names(Selections))]
  
  
  savename <- paste(rep.num,progenyinfo$namestem,prefix,"-phenos.txt",sep="")
  write.table(phenos,file=savename,col.names=F,row.names=T,quote=F)
  
  pheno.info<-list(phenos=phenos, E=E, TGV=geneticvalues, selections=Selections, 
                   selection.phenos=Selection.phenos,ped=ped, prog.inbred.level=prog.inbreeding,
                   genos.3d=new.parent.genos, numparents=numselections,select.genval=new.pars.genval, pp.genval=pp.genval,
                   par.ids=select.ids,select.ped.ids=select.ped.ids,all.markers=all.markers,all.phenos=all.phenos)
  cat("The returned object is a list containing a matrix of phenotypic data with\n")
  cat("the specified heritability, a vector of unscaled true genetic values,\n")
  cat("to the same total variance as the phenotypic values, and a vector of\n" )
  cat("selected individuals, four per family, with the highest phenotype value.\n" )
  cat("Both phenotypes and true genetic values include values for parents first followed by progeny.\n")
  
  return(pheno.info)
  # # # # #   End of createPhenos() function # # # # # # # #
}