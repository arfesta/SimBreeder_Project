###Set Working Directory, Save loctation, & Source file####
setwd("/Volumes/200GB/Google Drive/Adam Festa/") #MAC
#setwd("C:/Users/PC/Google Drive/Adam Festa")#Windows
#setwd("K:/Desktop")
Rep="data/miniTest" #SBMP or SELF
source("~/Desktop/FINAL.create.map.R")
#source("scripts/g.simulator.with.gblup.edit.R") #or 192 simulator to mimic Harry Wu paper 
Save=FALSE

###User Defined Inputs for Map####
Haldane=TRUE        # Use either Haldane Function or 
map.length=1800     # Total map length in centimorgans
num.chromos=12      # Number of chromosomes
total.loci= 1200   # Total number of loci to be distributed to equally distributed to chromosomes (must be divisble by num.chromos)
TotalQTL=1000       # Total number of QTL involved in a trait
NumSNPQTL=960     # Total number of QTL which are affected by dominance (must be less than TotalQTL)
dom.coeff=1        # Dominance level: Range (0-1)
NumParents=50       # Number of parents which are created to be inital founder population
NumMarkers=10   # Number of markers that are used in breeding strategies 
Heritability=c(.3)  # Heritablity of the trait under consideration
set.seed(23423)     # Set seed so that snp qtl MAF are reproducible
SNPQTLMAFs <- runif(NumSNPQTL,min=.01,max=.02)
#SNPeffects= rep(0, 100) # This can be set to specify the effects of SNPS (length SNP effects must match NumSNPQTL)
#useSNPeffects=T         # Operator switch for using the SNP effects

###Indicating the type of Mating Design and Selection Method#####
SBPM=F
SELFL=F
SELFP=F
DP=F
SBPW=F
GBLUP=F
ABLUP=F
SELFTEST=F
HD=F

###Create Founder Population####
map<- createMap()
parents <- createParents()
parents.TGV <- createTGV(prefix="parents.TGV")
parents.phenos <- createUnscaledParentPhenos()

#####Simulate Generations####
CrossFile1="g.study.1stgen.txt"
#CrossFile1="/Volumes/200GB/Google Drive/Adam Festa/self.txt"
progeny1 <- createProgeny(prefix="progeny1.pop")
progeny1.TGV <- createTGV(prefix="progeny1.TGV",data1=progeny1)
NumSelections=2
SBPW=T
progeny1.phenos <- createUnscaledPhenos(progenyinfo=progeny1,progenyTGV=progeny1.TGV, prefix="selection1.pop") #Selections made here based on phenotype

#CrossFile1="~/Google Drive/Adam Festa/g.study.2ndgen.txt"
CrossFile1="g.3nd.gen.ped.txt"
SBPW=F
progeny2 <- createProgeny(parentinfo=progeny1.phenos, gen=2, prefix="progeny2.pop")
progeny2.TGV <- createTGV(prefix="progeny2.TGV",data1=progeny2)
SBPW=T
progeny2.phenos <- createUnscaledPhenos(gen=2,progenyinfo=progeny2, progenyTGV=progeny2.TGV, parentphenos = progeny1.phenos, prefix="selection2.pop")

pdf("~/Desktop/practice2.pdf",width=21,height=7,pointsize=30)
old.par <- par(mfrow=c(1,3))
####### Plot Straightness#########
#map.length=1800     # Total map length in centimorgans
#num.chromos=12      # Number of chromosomes
#total.loci= 1200   # Total number of loci to be distributed to equally distributed to chromosomes (must be divisble by num.chromos)
#TotalQTL=1000        # Total number of QTL involved in a trait
#NumSNPQTL=900       # Total number of QTL which are affected by dominance (must be less than TotalQTL)
#dom.coeff=0  
#A = 1
#a=0
#user defined snpqtl MAF
r1 <- rep(NA,10)
r2<- rep(NA,10)
r3<- rep(NA,10)
r4<- rep(NA,10)
a <-1
for (i in 1:10){
  r1[i] <- mean(progeny2.phenos$phenos[a:(a+49)]) # outbred
  r2[i] <- mean(progeny2.phenos$phenos[(a+50):(a+99)]) #half-sib
  r3[i] <- mean(progeny2.phenos$phenos[(a+100):(a+149)]) # full-sib
  r4[i] <- mean(progeny2.phenos$phenos[(a+150):(a+199)]) # self
  a <- a+200}

outcross <- rep(0,10)
half<- rep(.125,10)
full<- rep(.25,10)
self<- rep(.5,10)
coastal.strt <- data.frame(par=rep(1:10,times=4),x=c(outcross,half,full,self),y=c(r1,r2,r3,r4))
# Create Line Chart
a<-.37
b<-.6
range01 <- function(x){(((b-a)*(x-min(x)))/(max(x)-min(x)))+a}
cy <- range01(coastal.strt$y)
coastal.strt$y <- cy
# convert factor to numeric for convenience 
ntrees <- max(coastal.strt$par)

# get the range for the x and y axis 
xrange <- c(0.0,0.5)
yrange <- c(0,1)

# set up the plot 
plot(xrange, yrange, axes=F,type="n", xlab="Inbreeding Coefficient",
     ylab="Probability of Straightness", bty='L',main="Straightness" ) 
color <- rainbow(ntrees)
color[3] <- colors()[641]
linetype <- c(1:ntrees) 
plotchar <- seq(1,1+ntrees,1)

# add lines 
for (i in 1:ntrees) { 
  tree <- subset(coastal.strt, par==i) 
  lines(tree$x, tree$y, type="b", lwd=2,
        lty=linetype[i], col=color[i], pch=plotchar[i], cex=.5) 
} 
xticks <- c(0,.125,.25,.5)
yticks <- seq(0.0,1.0,.2)
axis(2, at = yticks, labels = yticks, col.axis="black",las=1)
axis(1, at = xticks, labels = xticks, col.axis="black",las=1,lwd.ticks = .5)
par(xpd=TRUE)
####### Plot Height#########
#totalloci=1200
#TotalQTL=1000
#NumSNPQTL=950 #SNPs involved in trait 
#dom.coeff=1 #Range (0-1)
#user MAF specification
#A 1
#a -1.5
r1 <- rep(NA,10)
r2<- rep(NA,10)
r3<- rep(NA,10)
r4<- rep(NA,10)
a <-1
for (i in 1:10){
  r1[i] <- mean(progeny2.phenos$phenos[a:(a+49)]) # outbred
  r2[i] <- mean(progeny2.phenos$phenos[(a+50):(a+99)]) #half-sib
  r3[i] <- mean(progeny2.phenos$phenos[(a+100):(a+149)]) # full-sib
  r4[i] <- mean(progeny2.phenos$phenos[(a+150):(a+199)]) # self
  a <- a+200}

outcross <- rep(0,10)
half<- rep(.125,10)
full<- rep(.25,10)
self<- rep(.5,10)
coastal.ht <- data.frame(par=rep(1:10,times=4),x=c(outcross,half,full,self),y=c(r1,r2,r3,r4))
# Create Line Chart
a<-6.75
b<-11
range01 <- function(x){(((b-a)*(x-min(x)))/(max(x)-min(x)))+a}
cy <- range01(coastal.ht$y)
coastal.ht$y <- cy
# convert factor to numeric for convenience 
ntrees <- max(coastal.ht$par)

# get the range for the x and y axis 
xrange <- c(0.0,0.5)
yrange <- range(coastal.ht$y)

# set up the plot 
plot(xrange, yrange, axes=F,type="n", xlab="Inbreeding Coefficient",
     ylab="Height (m)", bty='L', main= "Height") 
color <- rainbow(ntrees)
color[3] <- colors()[641]
linetype <- c(1:ntrees) 
plotchar <- seq(1,1+ntrees,1)

# add lines 
for (i in 1:ntrees) { 
  tree <- subset(coastal.ht, par==i) 
  lines(tree$x, tree$y, type="b", lwd=2,
        lty=linetype[i], col=color[i], pch=plotchar[i], cex=.5) 
} 

xticks <- c(0,.125,.25,.5)
yticks <- seq(7,11,1)
axis(2, at = yticks, labels = yticks, col.axis="black",las=1)
axis(1, at = xticks, labels = xticks, col.axis="black",las=1,lwd.ticks = .5)
par(xpd=TRUE)
#legend("bottomleft", c("P1","P2","P3","P4","P5",
#                               "P6","P7","P8","P9","P10"), cex=0.8, col=color,
#       pch=plotchar, lty=linetype, title="Parent")
####### Plot Volume#########
#totalloci=1200
#TotalQTL=1000
#NumSNPQTL=960 #SNPs involved in trait  #last
#dom.coeff=1 #Range (0-1)
#user defined snpsqtl MAF
#A is 1
#a is -.1

r1 <- rep(NA,10)
r2<- rep(NA,10)
r3<- rep(NA,10)
r4<- rep(NA,10)
a <-1
for (i in 1:10){
  r1[i] <- mean(progeny2.phenos$phenos[a:(a+49)]) # outbred
  r2[i] <- mean(progeny2.phenos$phenos[(a+50):(a+99)]) #half-sib
  r3[i] <- mean(progeny2.phenos$phenos[(a+100):(a+149)]) # full-sib
  r4[i] <- mean(progeny2.phenos$phenos[(a+150):(a+199)]) # self
  a <- a+200}

outcross <- rep(0,10)
half<- rep(.125,10)
full<- rep(.25,10)
self<- rep(.5,10)
coastal.vol <- data.frame(par=rep(1:10,times=4),x=c(outcross,half,full,self),y=c(r1,r2,r3,r4))
# Create Line Chart
a<-38
b<-140
range01 <- function(x){(((b-a)*(x-min(x)))/(max(x)-min(x)))+a}
cy <- range01(coastal.vol$y)
coastal.vol$y <- cy
# convert factor to numeric for convenience 
ntrees <- max(coastal.vol$par)

# get the range for the x and y axis 
xrange <- c(0.0,0.5)
yrange <- range(coastal.vol$y)

# set up the plot 
plot(xrange, yrange, axes=F,type="n", xlab="Inbreeding Coefficient",
     ylab="Volume (cubic dm)", bty='L',main="Volume" ) 
color <- rainbow(ntrees)
color[3] <- colors()[641]
linetype <- c(1:ntrees) 
plotchar <- seq(1,1+ntrees,1)

# add lines 
for (i in 1:ntrees) { 
  tree <- subset(coastal.vol, par==i) 
  lines(tree$x, tree$y, type="b", lwd=2,
        lty=linetype[i], col=color[i], pch=plotchar[i], cex=.5) 
} 
xticks <- c(0,.125,.25,.5)
yticks <- seq(40,140,20)
axis(2, at = yticks, labels = yticks, col.axis="black",las=1)
axis(1, at = xticks, labels = xticks, col.axis="black",las=1,lwd.ticks = .5)
par(xpd=TRUE)
#legend("bottomleft", c("P1","P2","P3","P4","P5",
#                               "P6","P7","P8","P9","P10"), cex=0.8, col=color,
#       pch=plotchar, lty=linetype, title="Parent")

par(old.par)
dev.off()



#########mean###########
q <- 1
a<- rep(NA,10)
b<- rep(NA,10)
c<- rep(NA,10)
d<- rep(NA,10)
#range01 <- function(x){(x-min(x*1.25))/(max(x*.9)-min(x*1.25))}
#cy <- range01(coastal$y)
for (i in 1:10){
  a[i] <- (coastal$y[q])
  b[i] <- (coastal$y[q+10])
  c[i] <- (coastal$y[q+20])
  d[i] <- (coastal$y[q+30])
  q <- q+1
}

mean(a/a-1)
mean(b/a-1)
mean(c/a-1)
mean(d/a-1)


