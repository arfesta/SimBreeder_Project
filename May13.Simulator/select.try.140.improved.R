####Extract Selections####
extract.selections <- function(num.cores=1, mapinfo=map, crossdesign=cross.design, past.tgv, past.phenos, parentinfo,
                               progenyinfo, progenyTGV, progenyphenos, selection.strategy, 
                               among.family.selection = NULL, within.family.selection = NULL, numSelections.within.family = 1,
                               numSelections.among.family = NULL, relmatrix = "pedigree",
                               prefix = NULL, rep.num = NULL, reduced = F, weighted=F) {
  
 library(MatrixModels); library(parallel); library(pedigreemm); library(pedigree)
  
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
      parent.markers <- parentinfo$markers.matrix[as.numeric(names(parentinfo$pars)),]
    } else { parent.markers <- past.tgv$markers.matrix}} else {parent.markers <- past.phenos$all.markers}
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
    all.markers <- rbind(parent.markers,prog.markers)
    if(reduced){  all.phenos <- c(parentinfo$phenos,progeny.phenos)
    all.genetic.vals <- c(parentinfo$genetic.values,prog.genetic.values)
    } else { all.phenos <- c(past.phenos$phenos,progeny.phenos)
    all.genetic.vals <- c(past.phenos$genetic.values,prog.genetic.values)}
  } else {
    ped <- pedigree(label = selection.ped[,1],sire = selection.ped[,2],dam = selection.ped[,3])
    all.markers <- rbind(parent.markers,prog.markers)
    all.phenos <- c(parentinfo$all.phenos,progeny.phenos)
    all.genetic.vals <- c(past.phenos$genetic.values,prog.genetic.values)
  }
  #ped <- ped[order(as.numeric(ped$ID)),]
  
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
  
  if (selection.strategy == "self.test") {
    class(all.markers) <- "numeric"
    if(withininfamilygblup){
      progeny.blups <- lapply(1:length(cross.design[,1]),function(x){
        par1 <- cross.design[x,1]
        par2 <- cross.design[x,2]
        prog <- which(ped[,2] %in% par1 & ped[,3] %in% par2)
        fam.ped <- ped[prog,]
        fam.phenos <- all.phenos[match(fam.ped[,1],names(all.phenos))]
        fam.markers <- all.markers[match(fam.ped[,1],rownames(all.markers)),]
        c.ped <- cbind(fam.ped[,1:3],y=fam.phenos)
        colnames(c.ped) <- c("ID","SIRE","DAM","y")
        c.ped[,1] <- as.integer(c.ped[,1])
        c.ped[,2] <- as.numeric(c.ped[,2])
        c.ped[,3] <- as.numeric(c.ped[,3])
        c.ped[,4] <- as.numeric(c.ped[,4])
        c.ped <- data.frame(c.ped)
        sol <- gblup(y ~ 1,data=c.ped[,c(1,4)],M=fam.markers,lambda = 9)
        g <- which( rownames(sol) %in% names(progenyphenos$phenos))
        sol[g,1]
      })
      progeny.blups <- unlist(progeny.blups)
      names(progeny.blups) <- names(progenyphenos$phenos)
    } else {
      
      c.ped <- cbind(ped[,1:3],y=all.phenos)
      c.ped <- data.frame(c.ped)
      colnames(c.ped) <- c("ID","SIRE","DAM","y")
      c.ped[,1] <- as.integer(c.ped[,1])
      c.ped[,2] <- as.numeric(c.ped[,2])
      c.ped[,3] <- as.numeric(c.ped[,3])
      c.ped[,4] <- as.numeric(c.ped[,4])
      sol <- gblup(y ~ 1,data=c.ped[,c(1,4)],M=fam.markers,lambda = 9)
      g <- match(names(progenyphenos$phenos),rownames(sol))
      progeny.blups <- sol[g,1] ; names(progeny.blups) <- names(progenyphenos$phenos)}
    
    
    for (each.cross in 1:length(top.280.families)){
      yo <- vector(mode="list",length=as.numeric(cross.design[1,3]))
      yoo <- vector(mode="list",length=cross.design[1,3])
      yo.1 <- vector(mode="list",length=cross.design[1,3])
      yoo.1 <- vector(mode="list",length=cross.design[1,3])
      cross <- top.280.families[each.cross]
      parent1 <- as.numeric(cross.design[cross,1])
      parent2 <- as.numeric(cross.design[cross,2])
      print(parent2)
      if(reduced){first.in.family <- (as.numeric(cross.design[cross,3])*cross) - prog.percross + 1
      num.prog <- seq(first.in.family,first.in.family + prog.percross -1)} else {
        num.prog <- seq(first.in.family,first.in.family + prog.percross -1)}
      parent1.alleles <-  progenyinfo$genos.3d[as.numeric(names(Selfed.progeny[[which(names(Selfed.progeny) %in% parent1)]]$top.indiv)),num.prog,1]
      parent2.alleles <-  progenyinfo$genos.3d[as.numeric(names(Selfed.progeny[[which(names(Selfed.progeny) %in% parent2)]]$top.indiv)),num.prog,2]
      top.parent.1 <- Selfed.progeny[[which(names(Selfed.progeny) %in% parent1)]]$top.indiv
      top.parent.2 <- Selfed.progeny[[which(names(Selfed.progeny) %in% parent2)]]$top.indiv
      
      for(i in 1:prog.percross){
        yo[[i]] <- which(parent1.alleles[,i]==top.parent.1)
        yoo[[i]] <- sum(Selfed.progeny[[which(names(Selfed.progeny) %in% as.character(parent1))]]$most.imp[yo[[i]]])
        yo.1[[i]] <- which(parent2.alleles[,i]==top.parent.2)
        yoo.1[[i]] <- sum(Selfed.progeny[[which(names(Selfed.progeny) %in% as.character(parent2))]]$most.imp[yo.1[[i]]])
      }
      total <- mapply(function(x1,y1) x1[[1]] + y1[[1]], yoo, yoo.1)
      select.from.fam <- names(sort(rank(total)*0 + rank(-progeny.blups[num.prog]))[1])
      select.from.fam.if.g <- names(sort(rank(total)*0 + rank(-gprogeny.blups[num.prog]))[1])
      Selections <- c(Selections,select.from.fam)
      if.within.g.Selections <- c(if.within.g.Selections,select.from.fam.if.g)
      first.in.family <- max(num.prog) + 1
    }
    g <- which( names(progenyphenos$phenos) %in% Selections)
    Selections <- g
  }
  
  if (among.family.selection == "ABLUP"){
    #pedd <- data.frame(label = full.ped[,1],sire = full.ped[,2],dam = full.ped[,3],P=all.phenos[match(full.ped[,1],names(all.phenos))])
    #h.2 <- var(prog.genetic.values)/var(progeny.phenos)
    #lambda <- (1-h.2)/h.2 
    #the.blups <- blup(P ~ 1 ,ped = pedd,alpha = lambda)
    Ai <- as.matrix(getAInv(ped))
    l <- which(rownames(Ai) %in% names(progeny.phenos))
    Ai <- Ai[l,l]
    #the.blups <- the.blups[-1]
    #names(the.blups) <- names(all.phenos)[order(names(all.phenos))]
    #progeny.blups <- the.blups[match(names(all.phenos),names(the.blups))]
    #progeny.blups <- progeny.blups[-c(1:140)]
    s <-  ncol(Ai)/6
    e <- s-1
    first.seq <- seq(1,ncol(Ai),ncol(Ai)/6)
    SG1.2 <- Ai[c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e)),
                c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e))]
    SG1.3 <- Ai[c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e)),
                c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e))]
    SG1.4 <- Ai[c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e)),
                c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e))]
    SG1.5 <- Ai[c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e)),
                c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e))]
    
    SG1.6 <- Ai[c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e)),
                c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e))]
    
    SG2.3 <- Ai[c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e)),
                c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e))]
    SG2.4 <- Ai[c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e)),
                c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e))]
    SG2.5 <- Ai[c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e)),
                c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e))]
    SG2.6 <- Ai[c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e)),
                c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e))]
    
    SG3.4 <- Ai[c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e)),
                c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e))]
    
    SG3.5 <- Ai[c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e)),
                c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e))]
    SG3.6 <- Ai[c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e)),
                c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e))]
    SG4.5 <- Ai[c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e)),
                c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e))]
    SG4.6 <- Ai[c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e)),
                c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e))]
    SG5.6 <- Ai[c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e)),
                c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e))]
    rm(Ai); gc()
    s.phenos1 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e))]
    s.phenos2 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e))]
    s.phenos3 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e))]
    s.phenos4 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e))]
    s.phenos5 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e))]
    s.phenos6 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e))]
    s.phenos7 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e))]
    s.phenos8 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e))]
    s.phenos9 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e))]
    s.phenos10 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e))]
    s.phenos11 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e))]
    s.phenos12 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e))]
    s.phenos13 <- progeny.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e))]
    s.phenos14 <- progeny.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e))]
    s.phenos15 <- progeny.phenos[c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e))]
    
    pheno.list <- list(s.phenos1,s.phenos2,s.phenos3,s.phenos4,s.phenos5,s.phenos6,s.phenos7,s.phenos8,s.phenos9,
                       s.phenos10,s.phenos11,s.phenos12,s.phenos13,s.phenos14,s.phenos15)
    this.list <- list(SG1.2,SG1.3,SG1.4,SG1.5,SG1.6,SG2.3,SG2.4,SG2.5,SG2.6,SG3.4,SG3.5,SG3.6,SG4.5,SG4.6,SG5.6)
    pt <- proc.time()
    the.blups <- mclapply(X=1:length(this.list),function(x){
      the.data <- this.list[[x]]
      n.col <- ncol(the.data)
      h.2 <- var(prog.genetic.values)/var(progeny.phenos)
      lambda <- (1-h.2)/h.2 
      I <- diag(n.col)
      s.phenos <- pheno.list[[x]]
      
      sol <- solve(rBind(cBind(n.col, t(rep(1,n.col))),
                         cBind(rep(1,n.col), (I + (lambda*the.data)))),
                   matrix(c(sum(s.phenos),
                            c(as.vector(s.phenos)))))
      sol <- sol[-1,1]
      sol
    },mc.cores=length(this.list))
    proc.time()-pt
    SG1 <- apply(rbind(the.blups[[1]][1:(length(s.phenos1)/2)],
                       the.blups[[2]][1:(length(s.phenos1)/2)],
                       the.blups[[3]][1:(length(s.phenos1)/2)],
                       the.blups[[4]][1:(length(s.phenos1)/2)],
                       the.blups[[5]][1:(length(s.phenos1)/2)]),2,mean)
    
    SG2 <- apply(rbind(the.blups[[1]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[6]][1:(length(s.phenos1)/2)],
                       the.blups[[7]][1:(length(s.phenos1)/2)],
                       the.blups[[8]][1:(length(s.phenos1)/2)],
                       the.blups[[9]][1:(length(s.phenos1)/2)]),2,mean)
    
    SG3 <- apply(rbind(the.blups[[2]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[6]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[10]][1:(length(s.phenos1)/2)],
                       the.blups[[11]][1:(length(s.phenos1)/2)],
                       the.blups[[12]][1:(length(s.phenos1)/2)]),2,mean)
    
    SG4 <- apply(rbind(the.blups[[3]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[7]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[10]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[13]][1:(length(s.phenos1)/2)],
                       the.blups[[14]][1:(length(s.phenos1)/2)]),2,mean)
    
    SG5 <- apply(rbind(the.blups[[4]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[8]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[11]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[13]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[15]][1:(length(s.phenos1)/2)]),2,mean)
    
    SG6 <- apply(rbind(the.blups[[5]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[9]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[12]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[14]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                       the.blups[[15]][((length(s.phenos1)/2)+1):length(s.phenos1)]),2,mean)
    
    progeny.blups <- c(SG1,SG2,SG3,SG4,SG5,SG6); names(progeny.blups) <- names(progeny.phenos)
    
  }
  
  if (among.family.selection == "GBLUP"){
    A <- as.matrix(getA(ped))
    l <- which(rownames(A) %in% names(progeny.phenos))
    A <- A[l,l]
    
    m <- mapinfo$available.Markers
    parentmarkers <- parentinfo$genos.3d[m,,]
    progmarkers <- progenyinfo$genos.3d[m,,]
    allele1 <- cbind(parentmarkers[,,1],progmarkers[,,1])
    allele2 <- cbind(parentmarkers[,,2],progmarkers[,,2])
    #rm(progmarkers,parentmarkers)
    
    n.col <- ncol(allele1)
    g.rel.mat <- matrix(NA,nrow=n.col,ncol=n.col)
    yo <- mclapply(1:(n.col-1),FUN = function(x,cores){
      com <- allele1[,x]
      com2 <- allele2[,x]
      ss <- x + 1
      b <- sapply(ss:ncol(allele1),FUN=function(x){
        (length(which(com == allele1[,x])) + length(which(com2 == allele2[,x])))/240
      })
    },mc.cores = 60)
    
    diag(g.rel.mat) <- 1
    g.rel.mat[lower.tri(g.rel.mat)] <- unlist(yo)
    g.rel.mat <- t(g.rel.mat)
    g.rel.mat[lower.tri(g.rel.mat)] <- unlist(yo)
    rm(allele1,allele2); gc()
    end <- as.numeric(dim(g.rel.mat)[1]-dim(progenyinfo$genos.3d)[2])
    G.Rel.Mat <- g.rel.mat[-c(1:end),-c(1:end)]
    # get the decomposition: G = UDU'
    G.Rel.Mat1 <- .97*G.Rel.Mat + .03*A
    rm(G.Rel.Mat,A)
    
    out <- chol2inv(chol(G.Rel.Mat1))
    rm(G.Rel.Mat1); gc()
    h.2 <- var(prog.genetic.values)/var(progeny.phenos)
    lambda <- (1-h.2)/h.2 
    n.col <- ncol(out)
    I <- diag(n.col)
    sol <- solve(rBind(cBind(n.col, t(rep(1,n.col))),
                       cBind(rep(1,n.col), (I + (lambda*out)))),
                 matrix(c(sum(progeny.phenos),
                          c(as.vector(progeny.phenos)))))
    sol <- sol[-1,1]
    names(sol) <- names(progeny.phenos)
    progeny.blups <- sol
    gc()
  }
  
  if (among.family.selection == "Phenotype"){
    progeny.blups <- progenyphenos$phenos}
  
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
      A <- as.matrix(getA(ped))
      l <- which(rownames(A) %in% names(progeny.phenos))
      A <- A[l,l]
      
      m <- mapinfo$available.Markers
      if(gen==1){parentmarkers <- parentinfo$genos.3d[m,,]} else {parentmarkers <- parentinfo$all.markers}
      progmarkers <- progenyinfo$genos.3d[m,,]
      allele1 <- cbind(parentmarkers[,,1],progmarkers[,,1])
      allele2 <- cbind(parentmarkers[,,2],progmarkers[,,2])
      #rm(progmarkers,parentmarkers)
      cor.mat <- matrix(0,nrow=ncol(allele1),ncol=ncol(allele1))
      for(i in 1:120) {
        these <- unique(allele1[i,])
        for(all in 1:length(these)) {
          posit <- which(allele1[i,] == these[all])
          for(each in 1:length(posit)){ 
            fill <- posit[each]
            others <- posit[-c(1:each)]
            cor.mat[fill,others] <- cor.mat[fill,others] + 1/240
          }}}
      
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
      rm(allele1,allele2); gc()
      end <- as.numeric(dim(g.rel.mat)[1]-dim(progenyinfo$genos.3d)[2])
      G.Rel.Mat <- g.rel.mat[-c(1:end),-c(1:end)]
      # get the decomposition: G = UDU'
      G.Rel.Mat1 <- .97*G.Rel.Mat + .03*A
      rm(G.Rel.Mat,A); gc()
      
      
      s <-  ncol(G.Rel.Mat1)/6
      e <- s-1
      first.seq <- seq(1,ncol(G.Rel.Mat1),ncol(G.Rel.Mat1)/6)
      SG1.2 <- G.Rel.Mat1[c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e)),
                          c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e))]
      SG1.3 <- G.Rel.Mat1[c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e)),
                          c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e))]
      SG1.4 <- G.Rel.Mat1[c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e)),
                          c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e))]
      SG1.5 <- G.Rel.Mat1[c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e)),
                          c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e))]
      
      SG1.6 <- G.Rel.Mat1[c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e)),
                          c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e))]
      
      SG2.3 <- G.Rel.Mat1[c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e)),
                          c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e))]
      SG2.4 <- G.Rel.Mat1[c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e)),
                          c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e))]
      SG2.5 <- G.Rel.Mat1[c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e)),
                          c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e))]
      SG2.6 <- G.Rel.Mat1[c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e)),
                          c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e))]
      
      SG3.4 <- G.Rel.Mat1[c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e)),
                          c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e))]
      
      SG3.5 <- G.Rel.Mat1[c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e)),
                          c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e))]
      SG3.6 <- G.Rel.Mat1[c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e)),
                          c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e))]
      SG4.5 <- G.Rel.Mat1[c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e)),
                          c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e))]
      SG4.6 <- G.Rel.Mat1[c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e)),
                          c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e))]
      SG5.6 <- G.Rel.Mat1[c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e)),
                          c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e))]
      rm(G.Rel.Mat1); gc()
      s.phenos1 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[2]:(first.seq[2]+e))]
      s.phenos2 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[3]:(first.seq[3]+e))]
      s.phenos3 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos4 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos5 <- progeny.phenos[c(first.seq[1]:(first.seq[1]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos6 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[3]:(first.seq[3]+e))]
      s.phenos7 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos8 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos9 <- progeny.phenos[c(first.seq[2]:(first.seq[2]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos10 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[4]:(first.seq[4]+e))]
      s.phenos11 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos12 <- progeny.phenos[c(first.seq[3]:(first.seq[3]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos13 <- progeny.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[5]:(first.seq[5]+e))]
      s.phenos14 <- progeny.phenos[c(first.seq[4]:(first.seq[4]+e),first.seq[6]:(first.seq[6]+e))]
      s.phenos15 <- progeny.phenos[c(first.seq[5]:(first.seq[5]+e),first.seq[6]:(first.seq[6]+e))]
      
      pheno.list <- list(s.phenos1,s.phenos2,s.phenos3,s.phenos4,s.phenos5,s.phenos6,s.phenos7,s.phenos8,s.phenos9,
                         s.phenos10,s.phenos11,s.phenos12,s.phenos13,s.phenos14,s.phenos15)
      this.list <- list(SG1.2,SG1.3,SG1.4,SG1.5,SG1.6,SG2.3,SG2.4,SG2.5,SG2.6,SG3.4,SG3.5,SG3.6,SG4.5,SG4.6,SG5.6)
      pt <- proc.time()
      the.blups <- mclapply(X=1:length(this.list),function(x){
        the.data <- this.list[[x]]
        TRY <- the.data
        out <- chol2inv(chol(TRY))
        h.2 <- var(prog.genetic.values)/var(progeny.phenos)
        lambda <- (1-h.2)/h.2 
        n.col <- ncol(out)
        I <- diag(n.col)
        s.phenos <- pheno.list[[x]]
        
        sol <- solve(rBind(cBind(n.col, t(rep(1,n.col))),
                           cBind(rep(1,n.col), (I + (lambda*out)))),
                     matrix(c(sum(s.phenos),
                              c(as.vector(s.phenos)))))
        sol <- sol[-1,1]
        sol
      },mc.cores=length(this.list))
      proc.time()-pt
      SG1 <- apply(rbind(the.blups[[1]][1:(length(s.phenos1)/2)],
                         the.blups[[2]][1:(length(s.phenos1)/2)],
                         the.blups[[3]][1:(length(s.phenos1)/2)],
                         the.blups[[4]][1:(length(s.phenos1)/2)],
                         the.blups[[5]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG2 <- apply(rbind(the.blups[[1]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[6]][1:(length(s.phenos1)/2)],
                         the.blups[[7]][1:(length(s.phenos1)/2)],
                         the.blups[[8]][1:(length(s.phenos1)/2)],
                         the.blups[[9]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG3 <- apply(rbind(the.blups[[2]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[6]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[10]][1:(length(s.phenos1)/2)],
                         the.blups[[11]][1:(length(s.phenos1)/2)],
                         the.blups[[12]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG4 <- apply(rbind(the.blups[[3]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[7]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[10]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[13]][1:(length(s.phenos1)/2)],
                         the.blups[[14]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG5 <- apply(rbind(the.blups[[4]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[8]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[11]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[13]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[15]][1:(length(s.phenos1)/2)]),2,mean)
      
      SG6 <- apply(rbind(the.blups[[5]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[9]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[12]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[14]][((length(s.phenos1)/2)+1):length(s.phenos1)],
                         the.blups[[15]][((length(s.phenos1)/2)+1):length(s.phenos1)]),2,mean)
      
      gprogeny.blups <- c(SG1,SG2,SG3,SG4,SG5,SG6); names(gprogeny.blups) <- names(progeny.phenos)
      
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
  
  if(selection.strategy == "mass.selection") {
    sorted.EBVs <- sort(progeny.blups,decreasing = T)
    Selections <- match(names(sorted.EBVs[1:numSelections.among.family]),names(progeny.blups))
  }
  
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
  all.markers1 <- cbind(parentinfo$genos.3d[mapinfo$available.Markers,,1],new.parent.genos[mapinfo$available.Markers,,1]) 
  all.markers2 <- cbind(parentinfo$genos.3d[mapinfo$available.Markers,,2],new.parent.genos[mapinfo$available.Markers,,2])
  all.markers <- abind(all.markers1,all.markers2,along=3)
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
