extract_selections_ibd <- function(map.info, cross.design, past.tgv, past.phenos, parent.info, progeny.info, progeny.TGV, progeny.phenos,
                               selection.strategy, among.family.selection, within.family.selection = NULL, useSNPs="subset", useSNPfxn=0.25,
                               num.selections.within.family = 1, num.selections.among.family = NULL, 
                               relationship.matrix.type = rel.mat.cross, prefix = NULL, rep.num = NULL, 
                               reduced = F, weighted = F,num.cores = NULL) {
  
  library(MatrixModels)
  library(parallel)
  library(pedigreemm)
  library(data.table)
  library(tibble)
  library(tidyr)
  library(Rfast)
  
  {
    marker.loci <- which(map.info$types == "m")
    snpqtl.loci <- which(map.info$types == "snpqtl")
    snpqtl.map <- map.info[snpqtl.loci,]
    current.cross.design <- cross.design$cross.design
    prog.percross <- as.numeric(current.cross.design[1,3])
    prog.pedigree <- cross.design$progeny.pedigree
    generation <- as.numeric(prog.pedigree[1, 4])
    full.ped <- cross.design$full.pedigree
    selection.ped <- cross.design$selection.ped
    selection.snp.genos <- list() # Define an empty list to hold snpqtl genotypes of selected progeny
    
    if (generation == 1) {
      if (reduced) { 
        parent.markers <- parent.info$genos.3d[marker.loci, (names(parent.info$mean.parent.phenos)), ]
      }  else {
        parent.markers <- past.tgv$markers.matrix
      }
    } else { parent.markers <- parent.info$all.markers }
    prog.markers <- progeny.TGV$markers.matrix
    prog.phenos <- progeny.phenos$phenos
    prog.genetic.values <- progeny.TGV$genetic.values
    map.markers <- progeny.TGV$marker.map
    colnames(map.markers) <- c("chr", "pos")
    map.markers$chr <- as.character(map.markers$chr)
    if (generation == 1) {
      ped <- pedigree(label = full.ped[, 1], sire = full.ped[,2], dam = full.ped[, 3])
      if (reduced) {
        all.phenos <- c(parent.info$phenos, prog.phenos)
        all.genetic.vals <- c(parent.info$genetic.values, 
                              prog.genetic.values)
      }
      else {
        all.phenos <- c(past.phenos$phenos, progeny.phenos)
        all.genetic.vals <- c(past.phenos$genetic.values, 
                              prog.genetic.values)
      }
    } else {
      ped <- pedigree(label = selection.ped[, 1], sire = selection.ped[,2], dam = selection.ped[, 3])
      all.phenos <- c(parent.info$all.phenos, prog.phenos)
      all.genetic.vals <- c(past.phenos$genetic.values, prog.genetic.values)
    }
  }
  ### Selection Method ####
  if (among.family.selection == "ABLUP") {
    
    the.data <- as.matrix(getAInv(ped))
    the.data = the.data[i = match(names(prog.phenos), 
                                  colnames(the.data)), j = match(names(prog.phenos), 
                                                                 colnames(the.data))]
    
    n.col <- ncol(the.data)
    h.2 <- var(prog.genetic.values)/var(prog.phenos)
    lambda <- (1 - h.2)/h.2
    I <- diag(n.col)
    DD <- rbind(cbind(n.col, t(rep(1, n.col))), cbind(rep(1, n.col), (I + (lambda * the.data))))
    CC <- matrix(c(sum(prog.phenos),  c(as.vector(prog.phenos))))
    pt <- proc.time()
    sol <- solve(DD,CC)
    proc.time() - pt
    rm(I,the.data,DD,CC); gc(full=T,reset=T)
    sol <- sol[-1, 1]
    progeny.blups <- unlist(sol)
    names(progeny.blups) <- names(prog.phenos)
    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for (each in 1:nrow(current.cross.design)) {
      mean.progeny.blups <- c(mean.progeny.blups, mean(progeny.blups[first:last]))
      first <- last + 1
      last <- first + prog.percross - 1
    }
    sorted.mean.progeny.blups <- sort(mean.progeny.blups, 
                                      decreasing = T)
    top.families <- match(sorted.mean.progeny.blups, 
                              mean.progeny.blups)
    the.selections <- vector()
    first.in.family <- 1
    for (family in 1:length(current.cross.design[, 3])) {
      num.offspring <- as.numeric(current.cross.design[family, 
                                                       3])
      last.in.family <- num.offspring + first.in.family - 
        1
      temp <- (progeny.blups[first.in.family:last.in.family])
      sorted <- sort(temp, decreasing = TRUE)
      best.one <- which(temp == sorted[1:num.selections.within.family])
      selected <- best.one + first.in.family - 1
      the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1
    }
    if (length(num.selections.among.family) > 0) {
      the.selections <- the.selections[top.families[1:num.selections.among.family]]
    } else {
      the.selections <- the.selections[top.families]
    }
  }
  if (among.family.selection == "GBLUP") {
    prog.markers <- progeny.info$genos.3d[marker.loci, , ]
    gc(full=T,verbose = F,reset = T)
    
    all.m <- rep(0,ncol(prog.markers)*ncol(prog.markers))
    for(all.markers in 1:nrow(prog.markers)){
      unique.markers <-  unique(c(prog.markers[all.markers,,1],prog.markers[all.markers,,2]))
      
      for(each.marker in 1:length(unique.markers)){
        these.in.both <- which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] 
                               & prog.markers[all.markers,,2] %in% unique.markers[each.marker])
        if(length(these.in.both) > 1){
          test_in_both <- comb_n(n = these.in.both,k = 2)
          idx <- ((test_in_both[1,]-1)*ncol(prog.markers)) + test_in_both[2,]
          all.m[idx] <- all.m[idx] +2
          these.in.one <-  which(prog.markers[all.markers,,1] %in% unique.markers[each.marker]  | prog.markers[all.markers,,2] %in% unique.markers[each.marker] )
          test_in_one <- comb_n(n = these.in.one,k = 2)
          
          rm1 <- which(test_in_one[1,] %in% test_in_both[1,] & test_in_one[2,] %in% test_in_both[2,])
          rm2 <- which(test_in_one[2,] %in% test_in_both[1,] & test_in_one[1,] %in% test_in_both[2,])
          rmu <- unique(c(rm1,rm2))
          test_in_one <- test_in_one[,-c(rmu)]
          idx <- ((test_in_one[1,]-1)*ncol(prog.markers)) + test_in_one[2,]
          all.m[idx] <- all.m[idx] +1
        } else{
          these.in.one <-  which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] | prog.markers[all.markers,,2] %in% unique.markers[each.marker] )
          if(length(these.in.one) > 1){
            test_in_one <- comb_n(n = these.in.one,k = 2)
            idx <- ((test_in_one[1,]-1)*ncol(prog.markers)) + test_in_one[2,]
            all.m[idx] <- all.m[idx] +1
          } else {}
        }
      }
    }
    
    all.m <- all.m/(nrow(prog.markers)*2)
    g.mat <- (matrix(all.m, nrow = ncol(prog.markers), ncol = ncol(prog.markers),byrow = F))
    rm(all.m); gc(full=T,reset=T)
    diag(g.mat) <- 1
    cc <-  Rfast::transpose(g.mat)[upper.tri( Rfast::transpose(g.mat), diag = F)]
    g.mat <- Rfast::upper_tri.assign(x= g.mat,v = cc, diag = F); rm(cc)
    
    A <- as_tibble(as.matrix(getA(ped)))
    A = A[match(names(prog.phenos), colnames(A)), match(names(prog.phenos),  colnames(A))]
    the.data <- as.matrix(g.mat)*.99 + as.matrix(A) * 0.01
    rm(g.mat,A); gc(full=T,reset = T,verbose = F)
    the.data <- solve(the.data)
    n.col <- ncol(the.data)
    h.2 <- var(prog.genetic.values)/var(prog.phenos)
    lambda <- (1 - h.2)/h.2
    I <- diag(n.col)
    DD <- rbind(cbind(n.col, t(rep(1, n.col))), cbind(rep(1, n.col), (I + (lambda * the.data))))
    CC <- matrix(c(sum(prog.phenos),  c(as.vector(prog.phenos))))
    
    sol <- solve(DD,CC)
    rm(I,the.data,DD,CC); gc(full=T,reset=T,verbose = F)
    
    sol <- sol[-1, 1]
    progeny.blups <- unlist(sol)
    names(progeny.blups) <- names(prog.phenos)
    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for (each in 1:nrow(current.cross.design)) {
      mean.progeny.blups <- c(mean.progeny.blups, mean(progeny.blups[first:last]))
      first <- last + 1
      last <- first + prog.percross - 1
    }
    
    first <- 1
    last <- prog.percross
    mean.progeny.phenos <- vector()
    for (each in 1:nrow(current.cross.design)) {
      mean.progeny.phenos <- c(mean.progeny.phenos, mean(prog.phenos[first:last]))
      first <- last + 1
      last <- first + prog.percross - 1
    }
    
    sorted.mean.progeny.blups <- sort(mean.progeny.blups, decreasing = T)
    top.families <- match(sorted.mean.progeny.blups,  mean.progeny.blups)
    the.selections <- vector()
    first.in.family <- 1
    for (family in 1:length(current.cross.design[, 3])) {
      num.offspring <- as.numeric(current.cross.design[family,  3])
      last.in.family <- num.offspring + first.in.family -  1
      temp <- (progeny.blups[first.in.family:last.in.family])
      sorted <- sort(temp, decreasing = TRUE)
      best.one <- which(temp == sorted[1:num.selections.within.family])
      selected <- best.one + first.in.family - 1
      the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1
    }
    if (length(num.selections.among.family) > 0) {
      the.selections <- the.selections[top.families[1:num.selections.among.family]]
    } else {the.selections <- the.selections[top.families] }
  }
  if (among.family.selection == "snpibd") {
    	  # Subset the genetic map to include only the SNPQTL loci for use in calculating relationships, and remove the "chr" prefix from chr name column
		
		# Test if useSNPs=="subset"; if so, multiply useSNPfxn (default value 0.25) by length(snpqtl.loci) and take floor value as target # of loci
		# Calculate MAF of all SNPQTL loci (count of both "a" and "c" alleles) and filter for subset of loci with MAF>0.05 in progeny
		# Sample the target number of SNPs from the high-MAF subset of SNPs
		# SNPs from snpqtl.loci vector and use those for the genotype matrix. If useSNPs!="subset", use all SNPs in snpqtl.loci vector for genotype matrix.
		
		if(useSNPs=="subset") {
				SNP.subset.num <- floor(length(snpqtl.loci)*useSNPfxn)
				SNP.sample <- sample(snpqtl.loci,SNP.subset.num,replace=FALSE)
				snpqtl.map <- map.info[SNP.sample,]
				snpqtl.map$chr <- gsub("chr","",snpqtl.map$chr)
				snpqtl.loci <- SNP.sample # Assign snpqtl.loci to just the subset of SNPQTL loci used for relationship matrix
			} else {
				snpqtl.map <- map.info[snpqtl.loci,]
				snpqtl.map$chr <- gsub("chr","",snpqtl.map$chr)
			}

		
		# First make a matrix of father and mother IDs of the unique families
		families <- unique(paste(cross.design$progeny.pedigree[,2],cross.design$progeny.pedigree[,3],sep="_"))
		uniq.fam.list <- data.table::tstrsplit(families,"_")
		uniq.fam.mtx <- data.frame(male=uniq.fam.list[[1]],female=uniq.fam.list[[2]])
		
		# Define empty matrix with 6 columns and 0 rows of NA values; rbind lines of family pedigree/phenotype data to this matrix
		fam <- matrix(NA,nrow=0,ncol=6) 
		# Define empty matrix with 0 rows and enough cols of NA values for all snpqtl loci; cbind rows of space-separated allele calls to this matrix
		X <- matrix(NA,nrow=0,ncol=length(snpqtl.loci))

		# Recover genotypes for just the SNPQTL loci for progeny from transposed progeny1$genos.3d array of all loci; 
		# paste the two alleles together with space. Loci now in columns of allele calls; samples in rows. Assign column
		# names of progeny IDs from progeny1.phenos$phenos object
		progeny.snpqtl.genos <- matrix(paste(t(progeny.info$genos.3d[snpqtl.loci,,1]),t(progeny.info$genos.3d[snpqtl.loci,,2]),sep=" "),
				nrow(prog.pedigree),length(snpqtl.loci))
		dimnames(progeny.snpqtl.genos) <- list(names(progeny.phenos$phenos),snpqtl.map$loci)
		
		# Use a loop to reformat the genotypes of SNPQTL loci for parents and progeny into .ped format, cbind to pedigree and phenotype information
		# Code has to be different for generation 1 vs later generations because structure of input data objects is different
		if(generation==1){
			for(row in 1:nrow(uniq.fam.mtx)){
				# First extract the correct rows of cross.design$progeny.pedigree for the 60 progeny of each cross
				# Define the family number as the row number of uniq.fam.mtx
				# Append the male and female parent information and progeny data to bottom of fam matrix with rbind
				# Append the genotypes for 2 parents and 60 progeny to right side of X matrix of genotypes with cbind
				famID <- paste0("fam",row)
				par1 <- c(famID,uniq.fam.mtx[row,1],0,0,1,round(parent.info$phenos[uniq.fam.mtx[row,1]],0))
				par2 <- c(famID,uniq.fam.mtx[row,2],0,0,2,round(parent.info$phenos[uniq.fam.mtx[row,2]],0))
				
				tmp.fam <- cross.design$progeny.pedigree[intersect(which(cross.design$progeny.pedigree[,2]==uniq.fam.mtx[row,1]),
															which(cross.design$progeny.pedigree[,3]==uniq.fam.mtx[row,2])),1:3]
				indivIDs <- tmp.fam[,1]; maleIDs <- tmp.fam[,2]; femaleIDs <- tmp.fam[,3]
				sex <- sample(c(1,2),60,replace=TRUE); prog_phenos <- round(progeny1.phenos$phenos[indivIDs],0)
				progeny.fam <- cbind(rep(famID,60),indivIDs,maleIDs,femaleIDs,sex,prog_phenos)
				fam <- rbind(fam,par1,par2,progeny.fam)
				# Parent genotype array parent.info$genos.3d has loci in rows, parents in columns, and haplotypes in layers (3rd dimension)
				# Need to transpose each layer to have loci in columns and samples in rows, then paste layers for each sample together with space
				# and add that vector to the X matrix of genotypes 
				par1.genos <- paste(t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,1],1]),t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,1],2]),sep=" ")
				par2.genos <- paste(t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,2],1]),t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,2],2]),sep=" ")
				prog.genos <- progeny.snpqtl.genos[indivIDs,]
				X <- rbind(X,par1.genos,par2.genos,prog.genos)
			}
		} else {
			for(row in 1:nrow(uniq.fam.mtx)){
				# First extract the correct rows of cross.design$progeny.pedigree for the 60 progeny of each cross
				# Define the family number as the row number of uniq.fam.mtx
				# Append the male and female parent information and progeny data to bottom of fam matrix with rbind
				# Append the genotypes for 2 parents and 60 progeny to right side of X matrix of genotypes with cbind
				famID <- paste0("fam",row)
				par1 <- c(famID,uniq.fam.mtx[row,1],0,0,1,round(all.phenos[uniq.fam.mtx[row,1]],0))
				par2 <- c(famID,uniq.fam.mtx[row,2],0,0,2,round(all.phenos[uniq.fam.mtx[row,2]],0))
				
				tmp.fam <- cross.design$progeny.pedigree[intersect(which(cross.design$progeny.pedigree[,2]==uniq.fam.mtx[row,1]),
															which(cross.design$progeny.pedigree[,3]==uniq.fam.mtx[row,2])),1:3]
				indivIDs <- tmp.fam[,1]; maleIDs <- tmp.fam[,2]; femaleIDs <- tmp.fam[,3]
				sex <- sample(c(1,2),60,replace=TRUE); prog_phenos <- round(progeny.phenos$phenos[indivIDs],0)
				progeny.fam <- cbind(rep(famID,60),indivIDs,maleIDs,femaleIDs,sex,prog_phenos)
				fam <- rbind(fam,par1,par2,progeny.fam)
				# Parent genotype array parent.info$genos.3d has loci in rows, parents in columns, and haplotypes in layers (3rd dimension)
				# Need to transpose each layer to have loci in columns and samples in rows, then paste layers for each sample together with space
				# and add that vector to the X matrix of genotypes 
				par1.genos <- paste(t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,1],1]),t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,1],2]),sep=" ")
				par2.genos <- paste(t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,2],1]),t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,2],2]),sep=" ")
				prog.genos <- progeny.snpqtl.genos[indivIDs,]
				X <- rbind(X,par1.genos,par2.genos,prog.genos)
			}
		}
		# Join fam and X objects with cbind, write to space-delimited file with no dimnames.
		# Convert snpqtl.map object to PLINK format (chr, SNPid, cM_pos, nt_pos) by multiplying cM_pos x 1e6, save as space-delimited file.
		current.ped <- cbind(fam,X)
		write.table(current.ped,"plink_files/tmp.ped",row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
		plink.map <- data.frame(chr=snpqtl.map$chr,loci=snpqtl.map$loci,cMpos=snpqtl.map$pos,ntpos=(snpqtl.map$pos*1000000))
		write.table(plink.map,"plink_files/tmp.map",row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
		
		# Use system() to invoke plink v1.9 to convert text files to binary format
		system("/usr/local/bin/plink --file /media/hd02/Rsimulations/SimBreeder_Project/plink_files/tmp --out /media/hd02/Rsimulations/SimBreeder_Project/plink_files/tmp_bin",
				intern=FALSE)
				
		# load SNPRelate library, convert binary files to GDS format, then define GDS file as genofile and calculate KING-robust kinship 
		library(SNPRelate)
		snpgdsBED2GDS(bed.fn="plink_files/tmp_bin.bed",bim.fn="plink_files/tmp_bin.bim",
						fam.fn="plink_files/tmp_bin.fam", out.gdsfn="plink_files/tmp.gds")
		genofile <- snpgdsOpen("plink_files/tmp.gds")
		family.mtx <- matrix(unlist(strsplit((read.gdsn(index.gdsn(genofile,"sample.id"))),"-")),byrow=TRUE,ncol=2)
		king.ibd <- snpgdsIBDKING(genofile,type="KING-robust",family.id=family.mtx[,1])
		snpgdsClose(genofile)
				
		# Dimensions of king.ibd$kinship matrix are "num.parents + num.progeny" square
		# Indiv ids are in col 2 of fam object used to create .ped file, showing the 2 parents and 60 progeny IDs for each of 96 families.
		# Parent IDs are also in uniq.fam.mtx, in two columns - can make a vector of unique parent IDs and use that to remove parent rows
		# and columns from kinship matrix
		unique.parents <- unique(c(uniq.fam.mtx[[1]],uniq.fam.mtx[[2]]))
		unique.parent.rows <- which(fam[,2] %in% unique.parents)
		
		# Remove the parent rows from the king.ibd$kinship matrix and multiply by 2 to put on the same scale with g.mat
		king.mtx <- 2*(king.ibd$kinship[-unique.parent.rows,-unique.parent.rows])
		dimnames(king.mtx) <- list(rownames(fam)[-unique.parent.rows],rownames(fam)[-unique.parent.rows])

		l <- match(names(progeny.phenos$phenos), rownames(king.mtx))
		rel.mat <- king.mtx[l, l]
		rm(king.mtx);gc(full=T)

    A <- as_tibble(as.matrix(getA(ped)))
    A = A[match(names(prog.phenos), colnames(A)), match(names(prog.phenos),  colnames(A))]
    the.data <- as.matrix(rel.mat)*.99 + as.matrix(A) * 0.01
    rm(rel.mat,A); gc(full=T,reset = T,verbose = F)
    the.data <- solve(the.data)
    n.col <- ncol(the.data)
    h.2 <- var(prog.genetic.values)/var(prog.phenos)
    lambda <- (1 - h.2)/h.2
    I <- diag(n.col)
    DD <- rbind(cbind(n.col, t(rep(1, n.col))), cbind(rep(1, n.col), (I + (lambda * the.data))))
    CC <- matrix(c(sum(prog.phenos),  c(as.vector(prog.phenos))))
    
    sol <- solve(DD,CC)
    rm(I,the.data,DD,CC); gc(full=T,reset=T,verbose = F)
    
    sol <- sol[-1, 1]
    progeny.blups <- unlist(sol)
    names(progeny.blups) <- names(prog.phenos)
    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for (each in 1:nrow(current.cross.design)) {
      mean.progeny.blups <- c(mean.progeny.blups, mean(progeny.blups[first:last]))
      first <- last + 1
      last <- first + prog.percross - 1
    }
    
    first <- 1
    last <- prog.percross
    mean.progeny.phenos <- vector()
    for (each in 1:nrow(current.cross.design)) {
      mean.progeny.phenos <- c(mean.progeny.phenos, mean(prog.phenos[first:last]))
      first <- last + 1
      last <- first + prog.percross - 1
    }
    
    sorted.mean.progeny.blups <- sort(mean.progeny.blups, decreasing = T)
    top.families <- match(sorted.mean.progeny.blups,  mean.progeny.blups)
    the.selections <- vector()
    first.in.family <- 1
    for (family in 1:length(current.cross.design[, 3])) {
      num.offspring <- as.numeric(current.cross.design[family,  3])
      last.in.family <- num.offspring + first.in.family -  1
      temp <- (progeny.blups[first.in.family:last.in.family])
      sorted <- sort(temp, decreasing = TRUE)
      best.one <- which(temp == sorted[1:num.selections.within.family])
      selected <- best.one + first.in.family - 1
      the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1
    }
    if (length(num.selections.among.family) > 0) {
      the.selections <- the.selections[top.families[1:num.selections.among.family]]
    } else {the.selections <- the.selections[top.families] }
 }  
  if (among.family.selection == "Phenotype") {
    # If using Phenotypes as the among family selection, the progeny blups will be the phenos
    progeny.blups <- prog.phenos
    names(progeny.blups) <- names(prog.phenos) # Names of the progeny.blups will be the names of the progeny phenos
    
    first <- 1
    last <- prog.percross
    mean.progeny.blups <- vector()
    for (each in 1:nrow(current.cross.design)) {
      mean.progeny.blups <- c(mean.progeny.blups, mean(progeny.blups[first:last]))
      first <- last + 1
      last <- first + prog.percross - 1
    }
    sorted.mean.progeny.blups <- sort(mean.progeny.blups,decreasing = T)
    top.families <- match(sorted.mean.progeny.blups,mean.progeny.blups)
    the.selections <- vector()
    first.in.family <- 1
    for (family in 1:length(current.cross.design[, 3])) {
      num.offspring <- as.numeric(current.cross.design[family,3])
      last.in.family <- num.offspring + first.in.family -  1
      temp <- (progeny.blups[first.in.family:last.in.family])
      sorted <- sort(temp, decreasing = TRUE)
      best.one <- which(names(temp) == names(sorted[1:num.selections.within.family]))
      if(length(best.one) == 1){names(best.one) <- names(temp)[best.one]}
      selected <- best.one + first.in.family - 1
      the.selections <- c(the.selections, selected)
      first.in.family <- last.in.family + 1
    }
    if (!is.null(num.selections.among.family)) {
      the.selections <- the.selections[top.families[1:num.selections.among.family]]
    } else {
      the.selections <- the.selections[top.families]
    }    
  }
  
  {
    selection.phenos <- prog.phenos[names(the.selections)]
    new.pars.genval <- prog.genetic.values[names(the.selections)]
    selection.EBVs <- progeny.blups[names(the.selections)]
		# Recover genotypes for just the SNPQTL loci for progeny from transposed progeny.info$genos.3d array of all loci; 
		# paste the two alleles together with space. Loci now in columns of allele calls; samples in rows. Assign column
		# names of progeny IDs from progeny.phenos$phenos object
		progeny.snpqtl.genos <- matrix(paste(t(progeny.info$genos.3d[snpqtl.loci,,1]),t(progeny.info$genos.3d[snpqtl.loci,,2]),sep=" "),
				nrow(prog.pedigree),length(snpqtl.loci))
		dimnames(progeny.snpqtl.genos) <- list(names(progeny.phenos$phenos),snpqtl.map$loci)
		
		# From that matrix, extract just the rows that correspond to selected progeny - don't need the others. Save matrix to selection.snp.genos list
		current <- paste0("gen",generation)
		subset.snp.genos <- progeny.snpqtl.genos[the.selections,]
		selection.snp.genos[[current]] <- subset.snp.genos
    
    if(is.null(parent.info$phenos)){
      parent.info$phenos <- past.phenos$phenos
      parent.info$genetic.values <- past.phenos$genetic.values
    }
    if (generation == 1) {
      all.phenos <- c(parent.info$phenos, selection.phenos)
      all.genetic.vals <- c(parent.info$genetic.values, 
                            new.pars.genval)
    } else {
      all.phenos <- c(parent.info$all.phenos, selection.phenos)
      all.genetic.vals <- c(parent.info$all.genetic.vals, new.pars.genval)
    }
    ped <- full.ped[match(names(all.phenos), full.ped[, 1]), ]
    new.parent.genos <- progeny.info$genos.3d[, names(the.selections),]
    numselections <- dim(new.parent.genos)[2]
    select.ped.ids <- as.numeric(names(new.pars.genval))
    
    if(length(dim(new.parent.genos)) < 3){
      all.markers1 <- cbind(parent.info$genos.3d[marker.loci,,1], new.parent.genos[marker.loci, 1])
      colnames(all.markers1) <- c(colnames(parent.info$genos.3d[marker.loci,,1]),select.ped.ids)
      all.markers2 <- cbind(parent.info$genos.3d[marker.loci,,2], new.parent.genos[marker.loci, 2])
      colnames(all.markers2) <- c(colnames(parent.info$genos.3d[marker.loci,,2]),select.ped.ids)
      all.markers <- abind(all.markers1, all.markers2, along = 3)
    } else {
      all.markers1 <- cbind(parent.info$genos.3d[marker.loci,,1], new.parent.genos[marker.loci, ,1])
      all.markers2 <- cbind(parent.info$genos.3d[marker.loci,, 2], new.parent.genos[marker.loci, ,2])
    }
    all.markers <- abind(all.markers1, all.markers2, along = 3)
    f.ped <- pedigree(label = full.ped[, 1], sire = full.ped[,2], dam = full.ped[, 3])
    pedigree.inbreeding <- inbreeding(f.ped)
    names(pedigree.inbreeding) <- full.ped[, 1]
    progeny.inbreeding <- pedigree.inbreeding[prog.pedigree[,1]]
    selections.inbreeding <- progeny.inbreeding[names(the.selections)]
    progeny.gen.var <- var(prog.genetic.values[names(the.selections)])
    bulmer.effect <- var(prog.genetic.values[names(the.selections)]) -  var(prog.genetic.values)
    if(length(the.selections) == 1){
      out <- length(which(new.parent.genos[snpqtl.loci, 1] == "c"))
      out2 <- length(which(new.parent.genos[snpqtl.loci, 2] == "c"))
      result <- out + out2
    } else { 
      result <- sapply(rep(1:length(the.selections), 1), function(x) {
        out <- length(which(new.parent.genos[snpqtl.loci, , 1][, x] == "c"))
        out2 <- length(which(new.parent.genos[snpqtl.loci, , 2][, x] == "c"))
        outer <- out + out2
      })
      result <- unlist(result) }
  }
  
  
  if (relationship.matrix.type == "pedigree") {
    if (generation == 1) {
      ff <- getA(f.ped)
      A <- data.frame(as.matrix(ff),stringsAsFactors = F)
      colnames(A) <- rownames(A)
    }  else {
      selection.ped <- cross.design$selection.ped
      f.ped2 <- pedigree(label=selection.ped[,1],sire=selection.ped[,2],dam=selection.ped[,3])
      A <- getA(f.ped2)
      A <- data.frame(as.matrix(A),stringsAsFactors = F)
      colnames(A) <- rownames(A)
    }
    l <- match(names(selection.phenos), rownames(A))
    rel.mat <- A[l, l]
  }
  if (relationship.matrix.type == "markers") {
    allele1 <- all.markers[, , 1]
    allele2 <- all.markers[, , 2]
    allele1 <- allele1[, which(colnames(allele1) %in% names(selection.phenos))]
    allele2 <- allele2[, which(colnames(allele2) %in% names(selection.phenos))]
    all.alleles <- abind(allele1, allele2, along = 3)
    
    
    
    pt <- proc.time()
    all.m <- rep(0,ncol(prog.markers)*ncol(prog.markers))
    for(all.markers in 1:nrow(prog.markers)){
      unique.markers <-  unique(c(prog.markers[all.markers,,1],prog.markers[all.markers,,2]))
      
      for(each.marker in 1:length(unique.markers)){
        these.in.both <- which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] 
                               & prog.markers[all.markers,,2] %in% unique.markers[each.marker])
        if(length(these.in.both) > 1){
          test_in_both <- comb_n(n = these.in.both,k = 2)
          idx <- ((test_in_both[1,]-1)*ncol(prog.markers)) + test_in_both[2,]
          all.m[idx] <- all.m[idx] +2
          these.in.one <-  which(prog.markers[all.markers,,1] %in% unique.markers[each.marker]  | prog.markers[all.markers,,2] %in% unique.markers[each.marker] )
          test_in_one <- comb_n(n = these.in.one,k = 2)
          
          rm1 <- which(test_in_one[1,] %in% test_in_both[1,] & test_in_one[2,] %in% test_in_both[2,])
          rm2 <- which(test_in_one[2,] %in% test_in_both[1,] & test_in_one[1,] %in% test_in_both[2,])
          rmu <- unique(c(rm1,rm2))
          test_in_one <- test_in_one[,-c(rmu)]
          idx <- ((test_in_one[1,]-1)*ncol(prog.markers)) + test_in_one[2,]
          all.m[idx] <- all.m[idx] +1
        } else{
          these.in.one <-  which(prog.markers[all.markers,,1] %in% unique.markers[each.marker] | prog.markers[all.markers,,2] %in% unique.markers[each.marker] )
          if(length(these.in.one) > 1){
            test_in_one <- comb_n(n = these.in.one,k = 2)
            idx <- ((test_in_one[1,]-1)*ncol(prog.markers)) + test_in_one[2,]
            all.m[idx] <- all.m[idx] +1
          } else {}
        }
      }
    }
    proc.time() - pt
    all.m <- all.m/240
    g.mat <- matrix(all.m, nrow = ncol(all.alleles), ncol = ncol(all.alleles),byrow = F)
    rm(all.m);gc(full=T,reset=T)
    
    diag(g.mat) <- 1
    # Duplicate the comparisons that have been made to create the complete square matrix
    cc <-  Rfast::transpose(g.mat)[upper.tri( Rfast::transpose(g.mat), diag = F)]
    g.mat <- Rfast::upper_tri.assign(x= g.mat,v = cc, diag = F)
    rownames(g.mat) <- colnames(allele1)
    colnames(g.mat) <- colnames(allele1)
    l <- match(names(selection.phenos), rownames(g.mat))
    rel.mat <- g.mat[l, l]
    rm(g.mat);gc(full=T)
    
  }
  if (relationship.matrix.type == "snpibd") {
	  # Subset the genetic map to include only the SNPQTL loci for use in calculating relationships, and remove the "chr" prefix from chr name column
		
		# Test if useSNPs=="subset"; if so, multiply useSNPfxn (default value 0.25) by length(snpqtl.loci) and take floor value, sample that number of
		# SNPs from snpqtl.loci vector and use those for the genotype matrix. If useSNPs!="subset", use all SNPs in snpqtl.loci vector for genotype matrix.
		
		if(useSNPs=="subset") {
				SNP.subset.num <- floor(length(snpqtl.loci)*useSNPfxn)
				SNP.sample <- sample(snpqtl.loci,SNP.subset.num,replace=FALSE)
				snpqtl.map <- map.info[SNP.sample,]
				snpqtl.map$chr <- gsub("chr","",snpqtl.map$chr)
				snpqtl.loci <- SNP.sample # Assign snpqtl.loci to just the subset of SNPQTL loci used for relationship matrix
			} else {
				snpqtl.map <- map.info[snpqtl.loci,]
				snpqtl.map$chr <- gsub("chr","",snpqtl.map$chr)
			}

		
		# First make a matrix of father and mother IDs of the unique families
		families <- unique(paste(cross.design$progeny.pedigree[,2],cross.design$progeny.pedigree[,3],sep="_"))
		uniq.fam.list <- data.table::tstrsplit(families,"_")
		uniq.fam.mtx <- data.frame(male=uniq.fam.list[[1]],female=uniq.fam.list[[2]])
		
		# Define empty matrix with 6 columns and 0 rows of NA values; rbind lines of family pedigree/phenotype data to this matrix
		fam <- matrix(NA,nrow=0,ncol=6) 
		# Define empty matrix with 0 rows and enough cols of NA values for all snpqtl loci; cbind rows of space-separated allele calls to this matrix
		X <- matrix(NA,nrow=0,ncol=length(snpqtl.loci))
		# Recover genotypes for just the SNPQTL loci for progeny from transposed progeny1$genos.3d array of all loci; 
		# paste the two alleles together with space. Loci now in columns of allele calls; samples in rows. Assign column
		# names of progeny IDs from progeny1.phenos$phenos object
		progeny.snpqtl.genos <- matrix(paste(t(progeny.info$genos.3d[snpqtl.loci,,1]),t(progeny.info$genos.3d[snpqtl.loci,,2]),sep=" "),
				nrow(prog.pedigree),length(snpqtl.loci))
		dimnames(progeny.snpqtl.genos) <- list(names(progeny.phenos$phenos),snpqtl.map$loci)
		
		# Use a loop to reformat the genotypes of SNPQTL loci for parents and progeny into .ped format, cbind to pedigree and phenotype information
		# Code has to be different for generation 1 vs later generations because structure of input data objects is different
		if(generation==1){
			for(row in 1:nrow(uniq.fam.mtx)){
				# First extract the correct rows of cross.design$progeny.pedigree for the 60 progeny of each cross
				# Define the family number as the row number of uniq.fam.mtx
				# Append the male and female parent information and progeny data to bottom of fam matrix with rbind
				# Append the genotypes for 2 parents and 60 progeny to right side of X matrix of genotypes with cbind
				famID <- paste0("fam",row)
				par1 <- c(famID,uniq.fam.mtx[row,1],0,0,1,round(parent.info$phenos[uniq.fam.mtx[row,1]],0))
				par2 <- c(famID,uniq.fam.mtx[row,2],0,0,2,round(parent.info$phenos[uniq.fam.mtx[row,2]],0))
				
				tmp.fam <- cross.design$progeny.pedigree[intersect(which(cross.design$progeny.pedigree[,2]==uniq.fam.mtx[row,1]),
															which(cross.design$progeny.pedigree[,3]==uniq.fam.mtx[row,2])),1:3]
				indivIDs <- tmp.fam[,1]; maleIDs <- tmp.fam[,2]; femaleIDs <- tmp.fam[,3]
				sex <- sample(c(1,2),60,replace=TRUE); prog_phenos <- round(progeny1.phenos$phenos[indivIDs],0)
				progeny.fam <- cbind(rep(famID,60),indivIDs,maleIDs,femaleIDs,sex,prog_phenos)
				fam <- rbind(fam,par1,par2,progeny.fam)
				# Parent genotype array parent.info$genos.3d has loci in rows, parents in columns, and haplotypes in layers (3rd dimension)
				# Need to transpose each layer to have loci in columns and samples in rows, then paste layers for each sample together with space
				# and add that vector to the X matrix of genotypes 
				par1.genos <- paste(t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,1],1]),t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,1],2]),sep=" ")
				par2.genos <- paste(t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,2],1]),t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,2],2]),sep=" ")
				prog.genos <- progeny.snpqtl.genos[indivIDs,]
				X <- rbind(X,par1.genos,par2.genos,prog.genos)
			}
		} else {
			for(row in 1:nrow(uniq.fam.mtx)){
				# First extract the correct rows of cross.design$progeny.pedigree for the 60 progeny of each cross
				# Define the family number as the row number of uniq.fam.mtx
				# Append the male and female parent information and progeny data to bottom of fam matrix with rbind
				# Append the genotypes for 2 parents and 60 progeny to right side of X matrix of genotypes with cbind
				famID <- paste0("fam",row)
				par1 <- c(famID,uniq.fam.mtx[row,1],0,0,1,round(all.phenos[uniq.fam.mtx[row,1]],0))
				par2 <- c(famID,uniq.fam.mtx[row,2],0,0,2,round(all.phenos[uniq.fam.mtx[row,2]],0))
				
				tmp.fam <- cross.design$progeny.pedigree[intersect(which(cross.design$progeny.pedigree[,2]==uniq.fam.mtx[row,1]),
															which(cross.design$progeny.pedigree[,3]==uniq.fam.mtx[row,2])),1:3]
				indivIDs <- tmp.fam[,1]; maleIDs <- tmp.fam[,2]; femaleIDs <- tmp.fam[,3]
				sex <- sample(c(1,2),60,replace=TRUE); prog_phenos <- round(progeny.phenos$phenos[indivIDs],0)
				progeny.fam <- cbind(rep(famID,60),indivIDs,maleIDs,femaleIDs,sex,prog_phenos)
				fam <- rbind(fam,par1,par2,progeny.fam)
				# Parent genotype array parent.info$genos.3d has loci in rows, parents in columns, and haplotypes in layers (3rd dimension)
				# Need to transpose each layer to have loci in columns and samples in rows, then paste layers for each sample together with space
				# and add that vector to the X matrix of genotypes 
				par1.genos <- paste(t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,1],1]),t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,1],2]),sep=" ")
				par2.genos <- paste(t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,2],1]),t(parent.info$genos.3d[snpqtl.loci,uniq.fam.mtx[row,2],2]),sep=" ")
				prog.genos <- progeny.snpqtl.genos[indivIDs,]
				X <- rbind(X,par1.genos,par2.genos,prog.genos)
			}
		}
		# Join fam and X objects with cbind, write to space-delimited file with no dimnames.
		# Convert snpqtl.map object to PLINK format (chr, SNPid, cM_pos, nt_pos) by multiplying cM_pos x 1e6, save as space-delimited file.
		current.ped <- cbind(fam,X)
		write.table(current.ped,"plink_files/tmp.ped",row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
		plink.map <- data.frame(chr=snpqtl.map$chr,loci=snpqtl.map$loci,cMpos=snpqtl.map$pos,ntpos=(snpqtl.map$pos*1000000))
		write.table(plink.map,"plink_files/tmp.map",row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
		
		# Use system() to invoke plink v1.9 to convert text files to binary format
		system("/usr/local/bin/plink --file /media/hd02/Rsimulations/SimBreeder_Project/plink_files/tmp --out /media/hd02/Rsimulations/SimBreeder_Project/plink_files/tmp_bin",
				intern=FALSE)
				
		# load SNPRelate library, convert binary files to GDS format, then define GDS file as genofile and calculate KING-robust kinship 
		library(SNPRelate)
		snpgdsBED2GDS(bed.fn="plink_files/tmp_bin.bed",bim.fn="plink_files/tmp_bin.bim",
						fam.fn="plink_files/tmp_bin.fam", out.gdsfn="plink_files/tmp.gds")
		genofile <- snpgdsOpen("plink_files/tmp.gds")
		family.mtx <- matrix(unlist(strsplit((read.gdsn(index.gdsn(genofile,"sample.id"))),"-")),byrow=TRUE,ncol=2)
		king.ibd <- snpgdsIBDKING(genofile,type="KING-robust",family.id=family.mtx[,1])
		snpgdsClose(genofile)
				
		# Dimensions of king.ibd$kinship matrix are "num.parents + num.progeny" square
		# Indiv ids are in col 2 of fam object used to create .ped file, showing the 2 parents and 60 progeny IDs for each of 96 families.
		# Parent IDs are also in uniq.fam.mtx, in two columns - can make a vector of unique parent IDs and use that to remove parent rows
		# and columns from kinship matrix
		unique.parents <- unique(c(uniq.fam.mtx[[1]],uniq.fam.mtx[[2]]))
		unique.parent.rows <- which(fam[,2] %in% unique.parents)
		
		# Remove the parent rows from the king.ibd$kinship matrix and multiply by 2 to put on the same scale with g.mat
		king.mtx <- 2*(king.ibd$kinship[-unique.parent.rows,-unique.parent.rows])
		dimnames(king.mtx) <- list(rownames(fam)[-unique.parent.rows],rownames(fam)[-unique.parent.rows])

		l <- match(names(progeny.phenos$phenos), rownames(king.mtx))
		rel.mat <- king.mtx[l, l]
		rm(king.mtx);gc(full=T)

		}
  extraction.info <- list(relmat = rel.mat, delt.alleles = result, selection.snp.genos = selection.snp.genos,
                          selections = the.selections, bulmer.effect = bulmer.effect, 
                          select.EBVs = selection.EBVs, all.genetic.vals = all.genetic.vals, 
                          selection.phenos = selection.phenos, ped = ped, prog.inbred.level = progeny.inbreeding, 
                          select.inbred.level = selections.inbreeding, genos.3d = new.parent.genos, 
                          num.parents = numselections, select.genval = new.pars.genval, 
                          fullped = full.ped, select.ped.ids = select.ped.ids, 
                          all.markers = all.markers, all.phenos = all.phenos, cumulative.total = cross.design$cumul.total)
  cat("The returned object is a list containing a matrix of phenotypic data with\n")
  cat("the specified heritability, a vector of unscaled true genetic values,\n")
  cat("to the same total variance as the phenotypic values, and a vector of\n")
  cat("selected individuals, four per family, with the highest phenotype value.\n")
  cat("Both phenotypes and true genetic values include values for parents first followed by progeny.\n")
  return(extraction.info)
}
