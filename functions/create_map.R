create_genetic_map <- function(num.chromos, map.length, 
                               fullyinformed_markers=NULl, snp.markers, 
                               total.QTL, total.SNPQTL,
                               snp.marker.maf=NULL, snp.qtl.maf=NULL,shape_param = c(1,5),
                               distribute.loci = "even", snp.marker.distribution ="even", snp.qtl.distribution ="random",
                               map.dist = "haldane", chromosome.size.range=.2,
                               signif.digits=2, write.map=F) {
  
  total.markers = fullyinformed_markers + snp.markers
  total.loci = sum(total.markers,total.QTL)
  # If the length of the map is 1, then simluate chromsome sizes
  if(length(map.length) == 1){
    #Determine average length of a chromosome
    avg.chromo.size <- map.length / num.chromos
    
    #Use mean chromosome size length and the chromosomse size range variable to specify chromosome sizes
    all.chromo.sizes <- floor(runif(num.chromos, min=(avg.chromo.size - chromosome.size.range * avg.chromo.size),
                                    max=(avg.chromo.size + chromosome.size.range * avg.chromo.size)))
  }
  
  # If the length of the map size is a vector; ensure it is the same length as number of chromosomes and assign sizes from vector
  if(length(map.length) == num.chromos){
    all.chromo.sizes = map.length
  }
  #Determine the number of intervals that will be calculated between loci and the number of loci per chromosome
  #if(!is.null(total.loci)) {total.loci <- total.loci} else {total.loci <- sum(total.markers,total.QTL,total.SNPQTL) }
  
  num.intervals <- total.loci - 1 #Specifying the number of intervals that will need to be calcuated
  if(distribute.loci == "even") { 
    loci.per.chromo <- round(rep(total.loci / num.chromos, num.chromos),digits = 0)#Each chromosome will get an equal number of loci
  }
  if(is.numeric(distribute.loci)) {
    loci.per.chromo <- distribute.loci
  } 
  
  #Create Map and Recombination Frequencies####
  # Define map function used to create recombination frequencies:
  calc_rec_freq <- function(x,mapdist=map.dist,interval=intervals) {
    loci.intervals <- interval[[x]]
    if(map.dist == "haldane") { out <- (1 - (exp(-2 * loci.intervals))) / 2 }
    if(map.dist == "kosambi") { out <- ((exp(2 * loci.intervals) - exp(-2 * loci.intervals)) / (exp(2 * loci.intervals) + exp(-2 * loci.intervals))) / 2 }
    out
  }
  
  # Set up genetic map data frame
  list1 <- lapply(1:num.chromos,function(x) rep(paste("chr",x,sep=""),loci.per.chromo[x])) # List "chr#" character set for all loci
  list2 <- lapply(1:num.chromos,function(x) 1:loci.per.chromo[x]) # List locus# for all chromosomes
  intervals <- lapply(1:num.chromos,function(x) c(round(runif((loci.per.chromo[x]-1),min=round(0.2*all.chromo.sizes[x]/loci.per.chromo[x],1),
                                                              max=round(1.8*all.chromo.sizes[x]/loci.per.chromo[x],1)),signif.digits)/100,5))
  rec.freqs <- lapply(1:num.chromos,calc_rec_freq)
  #Recfreqs for each of the loci are determined by using the interval distance (for each locus) in Morgans and the haldane mapping function:
  #r=(1/2)*(1-e^(2*M))
  
  #Since intervals are in Morgans, we multiply by 100 to turn it back into cM for positions (cM is unit for chromosome lengths)
  positions <- lapply(1:num.chromos,function(x) {
    all.except.last <- c(round(100*cumsum(intervals[[x]][-length(intervals[[x]])]),digits=8))
    last <- abs(all.chromo.sizes[x]/100 + all.except.last[length(all.except.last)])
    c(all.except.last,last)})
  
  ## paste together locus names from list2 and list1
  locus.names<- paste(unlist(list1),unlist(list2),sep="_")
  
  #Create map that has locusnames, distance between loci, empty vector of types (SNPQTL,rQTL, or Marker), MAFs, position on chromosome, & rec freqs
  map <- data.frame(chr=unlist(list1),loci=locus.names, dist= (unlist(intervals)),types= (rep("NA",total.loci)), MAF=rep("NA",total.loci),
                    pos=as.numeric(unlist(positions)),recfreqs=as.numeric(unlist(rec.freqs)),stringsAsFactors = F)
  
  #Now sample from the map to specify SNPQTL & rQTL, the remainder are potential markers that can be used
  all.loci <- 1:total.loci  # vector that contains 1 through the number of all loci
  
  #Specify the last loci for each chromosome
  chromo.loci.index <- vector("list")
  for (i in 1:length(loci.per.chromo)) {
    if (i==1) {
      chromo.loci.index[[1]] <- loci.per.chromo[i]
    } else {
      chromo.loci.index[[i]] <- loci.per.chromo[[i]] + chromo.loci.index[[i-1]]
    }
  }
  chromo.loci.index <- unlist(chromo.loci.index)
  last.pos.chrs <-  sapply(1:num.chromos,function(x){round(map$pos[chromo.loci.index[x]-1],0)})
  
  f_inf_markers <- vector()
  if(length(fullyinformed_markers) != 0){
    
    markers_per_chr <- fullyinformed_markers/num.chromos
    
    start <- 1
    
    for(i in 1:num.chromos){
      dist.per.marker <- last.pos.chrs[i]/(markers_per_chr+1)
      map.posit <- map$pos[start:chromo.loci.index[i]]
      marker <- sapply(1:markers_per_chr,function(x){
        the.pos <- dist.per.marker*x
        which.min(abs(round(map.posit - the.pos,1)))
      })
      if(i > 1){
        marker <- marker + chromo.loci.index[i-1]
      }
      f_inf_markers <- c(f_inf_markers,marker)
      start <- chromo.loci.index[i] + 1
    } 
    map$types[f_inf_markers] <- "fi"
    map$MAF[f_inf_markers] <- .49
  }
  
  snp_markers <- c()
  if(snp.marker.distribution == "even") { 
    for(i in 1:num.chromos){
      snp_markers <- c(snp_markers,sample(which(map$chr == paste0("chr",i) & map$types != "fi"),snp.markers/num.chromos))
    }
  }
  map$types[snp_markers] <- "m"
  
  # Assign marker minor allele frequencies
  if(length(snp.marker.maf) == snp.markers){
    map$MAF[which(map$types == "m")] <- marker.maf
  }
  if(length(snp.marker.maf) == 1){
    map$MAF[which(map$types == "m")] <- snp.marker.maf
  }
  if(length(snp.marker.maf) == 2){
    y_beta <- rbeta(total.loci*5, shape1 = shape_param[1], shape2 = shape_param[2])
    temp_mafs <- y_beta[intersect(which(y_beta>=snp.marker.maf[1]),which(y_beta<snp.marker.maf[2]))]
    map$MAF[which(map$types=="m")] <- sample(temp_mafs, size = snp.markers, replace = TRUE)
  }
  
  # Assign SNP qtl locations 
  snpqtl <- vector()
  if(is.numeric(snp.qtl.distribution)) { 
    for(i in 1:num.chromos){
      snpqtl <- c(snpqtl,which(map$pos[-map$types] %in% sample(map$pos[which(map$chr == paste0("chr",i))],snp.qtl.distribution[i]) & map$chr == paste0("chr",i)))
    }
  }
  if(snp.qtl.distribution == "random"){
    snpqtl <- sample(all.loci[-c(f_inf_markers,snp_markers)],total.SNPQTL,replace=FALSE)
  }
  
  # Assign SNP qtl MAF
  if(is.null(snp.qtl.maf)){
    SNPQTL.MAFs <- runif(n=total.SNPQTL,.5,.5)
  } else {
    y_beta <- rbeta(total.loci*10, shape1 = 1, shape2 = 5)
    SNPQTL.MAFs <- y_beta[intersect(which(y_beta>=snp.qtl.maf[1]),which(y_beta<snp.qtl.maf[2]))][1:length(snpqtl)]
  }
  
  map$types[snpqtl] <- "snpqtl"   # Specify in the map data frame that these loci are snpqtl
  map$MAF[snpqtl] <- SNPQTL.MAFs   # Assign these loci the specificed minor allele frequencies generated by user
  
  # Assign rQTL and rQTL MAFs
  rQTL <- which(map$types == "NA")
  map$types[rQTL] <- "qtl"
  MAFs <- sample(rbeta(200000,.4,.4),total.loci,replace=F)
  MAFs[which(MAFs > .5)] <- 1- MAFs[which(MAFs > .5)]
  map$MAF[rQTL] <-    sample(MAFs,size = length(rQTL),replace=F) 
  
  #Now save map output:
  datevalue <- date()
  datevector <- unlist(strsplit(datevalue,"\\s"))
  timevector <- unlist(strsplit(datevector[4],":"))
  newtime <- paste(timevector[1],timevector[2],sep="h")
  newdate <- paste(datevector[3],datevector[2],datevector[5],sep="")
  namestem <- paste(newdate,newtime,sep="_")
  
  if(write.map){
    mapname<-paste(namestem,"_map.txt")
    write.table(map,file=mapname, quote=F, row.names=F, col.names=T, sep="\t")}
  map$MAF <- as.numeric(map$MAF)
  mapinfo<-list(genetic.map=map, total.loci.num=total.loci,total.QTL.num=total.QTL, total.SNPQTL.num=total.SNPQTL,
                QTLSNP.loci=sort(snpqtl),date.time=namestem, rQTL.loci=sort(rQTL), available.Markers=snp_markers,last.locus.per.chrom=chromo.loci.index)
  #cat("The returned object is a list containing:\n")
  #cat("Genetic map = a data.frame of 7 columns: chr, locus names, intervals, types (marker, snp, or qtl), MAFs, position, and recombination fractions;\n")
  return(mapinfo)
  # End of function #
}
