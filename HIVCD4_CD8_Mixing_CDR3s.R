# HIVCD4_CD8_Mixing.R v1
# 8th dec
# tried doing it with CDR3s, wasn't as big a difference as DCRs, switch back

# NB - INDEX MEANS ARE NOT CORRECT

# built on ExploreSorted.R (see MiSeq 19)

# Takes the decombined CD4,CD8 data (in same manner as rest of paper samples)
# and produces the graphs of the diversity observed when mixing populations

# Takes .exp files (produced by Expand.py) as input
# This allows this script to sample to different sizes from whole repertoire

#8th sept
# people kept finding the 0.25:1 value confusing, so going to change it to .31:1
  # (the actual mean) to make it easier to understand

defaults <- par("mar")
savedir <- "~/TCR/WRITE_UP/THESIS/WorkingPlots/"
savepath <- paste(savedir, gsub("-", " ", Sys.Date()), sep="")
len <- length

setwd("/media/jme/SAMSUNG/ThesisAnalysis/BLR")

library(ineq)
library(entropy)
library(vegan)

defaults <- par("mar")

savepath <- paste(savedir, gsub("-", " ", Sys.Date()), sep="")

# want a way to sample from whole repertoires in R, so I can more easily sample in a loop
# wrote Expand.py to take freq-collapsed file and write out one line per original detected moleculre

# 8th dec - realised old version had been working from DCRs, not CDR3s
# needs to be CDR3s to be comparable with other figures!

allfiles <- list.files(pattern = "\\.expcdr3")

files <- c()

for(i in allfiles){
  nam <- paste("exp", strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1], sep="")
  assign(nam, as.data.frame(read.csv(i, header=FALSE, stringsAsFactors=FALSE)))
  dat <- read.csv(i, header=FALSE, stringsAsFactors=FALSE)
  files <- c(files,nam)
  
}
remove(i)

hv1a4 <- as.data.frame(get("expHV01,CD4+a"))
hv1a8 <- as.data.frame(get("expHV01,CD8+a"))

hv1b4 <- as.data.frame(get("expHV01,CD4+b"))
hv1b8 <- as.data.frame(get("expHV01,CD8+b"))

d4a4 <- as.data.frame(get("expHVD4,CD4+a"))
d4a8 <- as.data.frame(get("expHVD4,CD8+a"))

d4b4 <- as.data.frame(get("expHVD4,CD4+b"))
d4b8 <- as.data.frame(get("expHVD4,CD8+b"))

###### FIX - 
# ALSO NEED TO CHANGE THE LINES PRODUCED

# These are the mean number of total CDR3s for alpha/beta chains
#   Only counting healthy, non,CD4/8 sorted data

# maxa <- 14222   DCRs (and original tags)
# maxb <- 19482

# CDR3s, extended tags
maxa = 19074
maxb = 23780

ratios <- c("1:0","3:1","2.5:1","2:1","1.5:1","1:1",".5:1","*","0:1")
# wil refer to these below as a,b,c,d,e,f,g,h like so:
#           a     b      c      d      e      f      g     h
# later added a .25:1 enty, can go unlabelled           g2

# 8th dec - can simplify all of this using table, as CDR3s are only one item per line (not 5!)

#hv1a

# a, 1:0
# smplhv1a4 <- hv1a4[sample(nrow(hv1a4),maxa),]
# cnthv1a4 <- aggregate(rep (1, len(smplhv1a4)), by = as.list(smplhv1a4), FUN = sum)
# aa = ineq(cnthv1a4)

smplhv1a4 <- hv1a4[sample(nrow(hv1a4),maxa),]
cnthv1a4 <- table(smplhv1a4)
aa = ineq(cnthv1a4)

# b, 3:1 = into 4
sample4 <- as.integer((maxa/4)*3)
sample8 <- as.integer((maxa/4))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
ba = ineq(c(cnthv1a4,cnthv1a8))

# c, 2.5:1
sample4 <- as.integer((maxa/3.5)*2.5)
sample8 <- as.integer((maxa/3.5))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
ca = ineq(c(cnthv1a4,cnthv1a8))

# d, 2:1
sample4 <- as.integer((maxa/3)*2)
sample8 <- as.integer((maxa/3))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
da = ineq(c(cnthv1a4,cnthv1a8))

# e, 1.5:1
sample4 <- as.integer((maxa/2.5)*1.5)
sample8 <- as.integer((maxa/2.5))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
ea = ineq(c(cnthv1a4,cnthv1a8))

# f, 1:1
sample4 <- as.integer((maxa/2))
sample8 <- as.integer((maxa/2))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
fa = ineq(c(cnthv1a4,cnthv1a8))

# g, .5:1
sample4 <- as.integer((maxa/3))
sample8 <- as.integer((maxa/3)*2)
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
ga = ineq(c(cnthv1a4,cnthv1a8))

# g2 v2, .31:1
sample4 <- as.integer(maxa*(.31/1.31))
sample8 <- as.integer(maxa*(1/1.31))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
g2a <- ineq(c(cnthv1a4,cnthv1a8))



# h, 0:1
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),maxa),]
cnthv1a8 <- table(smplhv1a8)
ha = ineq(cnthv1a8)

hv1agin <-c(aa,ba,ca,da,ea,fa,ga,g2a,ha)


#hv1b

# a, 1:0
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),maxb),]
cnthv1b4 <- table(smplhv1b4)
ab <- ineq(cnthv1b4)

# b, 3:1
sample4 <- as.integer((maxb/4)*3)
sample8 <- as.integer((maxb/4))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
bb <- ineq(c(cnthv1b4,cnthv1b8))

# c, 2.5:1
sample4 <- as.integer((maxb/3.5)*2.5)
sample8 <- as.integer((maxb/3.5))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
cb <- ineq(c(cnthv1b4,cnthv1b8))

# d, 2:1
sample4 <- as.integer((maxb/3)*2)
sample8 <- as.integer((maxb/3))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
db <- ineq(c(cnthv1b4,cnthv1b8))

# e, 1.5:1
sample4 <- as.integer((maxb/2.5)*1.5)
sample8 <- as.integer((maxb/2.5))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
eb <- ineq(c(cnthv1b4,cnthv1b8))

# f, 1:1
sample4 <- as.integer((maxb/2))
sample8 <- as.integer((maxb/2))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
fb <- ineq(c(cnthv1b4,cnthv1b8))

# g, .5:1
sample4 <- as.integer((maxb/3))
sample8 <- as.integer((maxb/3)*2)
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
gb <- ineq(c(cnthv1b4,cnthv1b8))

# g2 v2, .31:1
sample4 <- as.integer(maxb*(.31/1.31))
sample8 <- as.integer(maxb*(1/1.31))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
g2b <- ineq(c(cnthv1b4,cnthv1b8))


# h, 0:1
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),maxb),]
cnthv1b8 <- table(smplhv1b8)
hb <- ineq(cnthv1b8)

hv1bgin <-c(ab,bb,cb,db,eb,fb,gb,g2b,hb)

#d4a

# a, 1:0
smpld4a4 <- d4a4[sample(nrow(d4a4),maxa),]
cntd4a4 <- table(smpld4a4)
aa = ineq(cntd4a4)

# b, 3:1
sample4 <- as.integer((maxa/4)*3)
sample8 <- as.integer((maxa/4))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
ba = ineq(c(cntd4a4,cntd4a8))

# c, 2.5:1
sample4 <- as.integer((maxa/3.5)*2.5)
sample8 <- as.integer((maxa/3.5))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
ca = ineq(c(cntd4a4,cntd4a8))

# d, 2:1
sample4 <- as.integer((maxa/3)*2)
sample8 <- as.integer((maxa/3))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
da = ineq(c(cntd4a4,cntd4a8))

# e, 1.5:1
sample4 <- as.integer((maxa/2.5)*1.5)
sample8 <- as.integer((maxa/2.5))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
ea = ineq(c(cntd4a4,cntd4a8))

# f, 1:1
sample4 <- as.integer((maxa/2))
sample8 <- as.integer((maxa/2))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
fa = ineq(c(cntd4a4,cntd4a8))

# g, .5:1
sample4 <- as.integer((maxa/3))
sample8 <- as.integer((maxa/3)*2)
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
ga = ineq(c(cntd4a4,cntd4a8))

# g2 v2, .31:1
sample4 <- as.integer(maxa*(.31/1.31))
sample8 <- as.integer(maxa*(1/1.31))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
g2a <- ineq(c(cntd4a4,cntd4a8))


# h, 0:1
smpld4a8 <- d4a8[sample(nrow(d4a8),maxa),]
cntd4a8 <- table(smpld4a8)
ha = ineq(cntd4a8)

d4agin <-c(aa,ba,ca,da,ea,fa,ga,g2a,ha)


#d4b

# a, 1:0
smpld4b4 <- d4b4[sample(nrow(d4b4),maxb),]
cntd4b4 <- table(smpld4b4)
ab <- ineq(cntd4b4)

# b, 3:1
sample4 <- as.integer((maxb/4)*3)
sample8 <- as.integer((maxb/4))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
bb <- ineq(c(cntd4b4,cntd4b8))

# c, 2.5:1
sample4 <- as.integer((maxb/3.5)*2.5)
sample8 <- as.integer((maxb/3.5))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
cb <- ineq(c(cntd4b4,cntd4b8))

# d, 2:1
sample4 <- as.integer((maxb/3)*2)
sample8 <- as.integer((maxb/3))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
db <- ineq(c(cntd4b4,cntd4b8))

# e, 1.5:1
sample4 <- as.integer((maxb/2.5)*1.5)
sample8 <- as.integer((maxb/2.5))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
eb <- ineq(c(cntd4b4,cntd4b8))

# f, 1:1
sample4 <- as.integer((maxb/2))
sample8 <- as.integer((maxb/2))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
fb <- ineq(c(cntd4b4,cntd4b8))

# g, .5:1
sample4 <- as.integer((maxb/3))
sample8 <- as.integer((maxb/3)*2)
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
gb <- ineq(c(cntd4b4,cntd4b8))

# g2 v2, .31:1
sample4 <- as.integer(maxb*(.31/1.31))
sample8 <- as.integer(maxb*(1/1.31))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
g2b <- ineq(c(cntd4b4,cntd4b8))

# h, 0:1
smpld4b8 <- d4b8[sample(nrow(d4b8),maxb),]
cntd4b8 <- table(smpld4b8)
hb <- ineq(cntd4b8)

d4bgin <-c(ab,bb,cb,db,eb,fb,gb,g2b,hb)

svg(paste(savepath, "CD4 CD8 TCR mixing Gini.svg"), width=6, height=5, pointsize=15)

thick = 2
med = 2

par(mar=c(4.1,4.1,1.1,.1))
plot(rep(0,8), ylim=c(.1,.7), xlim=c(.5,8.5), ylab="Gini Index", xlab="Ratio CD4:CD8 TCRs", 
     axes=FALSE)#, main="D4/HV1")

# Plot mean HV/HIV values as bars in background - NB WILL REQUIRE UPDATING WHEN ALL IN

# hv mean + sds
hvm = 0.2891093949
hvsd = 0.0395112203

segments(.5,hvm,8.5,hvm, col=rgb(0,0,1,0.5), lwd=thick)
segments(.5,hvm+hvsd,8.5,hvm+hvsd, col=rgb(0,0,1,0.4), lwd=1)
segments(.5,hvm-hvsd,8.5,hvm-hvsd, col=rgb(0,0,1,0.4), lwd=1)

#h hiv mean
hivm = 0.5681027253
hivsd = 0.080572715


segments(.5,hivm,8.5,hivm, col=rgb(1,0,0,0.5), lwd=thick)
segments(.5,hivm+hivsd,8.5,hivm+hivsd, col=rgb(1,0,0,0.4), lwd=1)
segments(.5,hivm-hivsd,8.5,hivm-hivsd, col=rgb(1,0,0,0.4), lwd=1)

# ratios
hivrm = 8-(2*0.310987069)
hivrsd = 0.1486751688
segments(hivrm,.1,hivrm,.7, lwd=2, col=rgb(1,0,0,0.4), lty=3)
segments(hivrm+hivrsd,.1,hivrm+hivrsd,.7, lwd=1, lty=3, col=rgb(1,0,0,0.3))
segments(hivrm-hivrsd,.1,hivrm-hivrsd,.7, lwd=1, lty=3, col=rgb(1,0,0,0.3))

pts <- c(1:7,hivrm,8)

points(pts, hv1agin, pch=3, lwd=med, col="green")
points(pts, hv1bgin, pch=3, lwd=med, col="purple")

points(pts, d4agin, pch=4, lwd=med, col="green")
points(pts, d4bgin, pch=4, lwd=med, col="purple")

lines(pts[2:9], hv1agin[2:9], lty=2, lwd=1, col=rgb(0,255,0,100, max=255))
lines(pts[2:9], hv1bgin[2:9], lty=2, lwd=1, col=rgb(160,32,240,100, max=255))

lines(pts[2:9], d4agin[2:9], lty=2, lwd=1, col=rgb(0,255,0,100, max=255))
lines(pts[2:9], d4bgin[2:9], lty=2, lwd=1, col=rgb(160,32,240,100, max=255))


legend("topleft", c("Alpha","Beta", "HV01", "HVD4", "Healthy ± sd", "HIV ± sd"), 
       pch=c(15,15,3,4,NA,NA), box.col="white",lty=c(NA,NA,NA,NA,1,1), 
       lwd=thick, col=c("green","purple", "black", "black", rgb(0,0,1,.5), rgb(1,0,0,.5)))


axis(2, at=c(.1,.2,.3,.4,.5,.6,.7))
axis(1, at=c(1,2,3,4,5,6,7,hivrm,8), labels=ratios, las=1, cex.axis=.9)
axis(1, at=hivrm, labels="x̄")#

box()

dev.off()

#########################################
######################################### ######################################### 
######################################### 

# SHANNON ENTROPY


#hv1a

# a, 1:0
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),maxa),]
cnthv1a4 <- table(smplhv1a4)
aa = entropy(cnthv1a4, unit="log2")

# b, 3:1 = into 4
sample4 <- as.integer((maxa/4)*3)
sample8 <- as.integer((maxa/4))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
ba = entropy(c(cnthv1a4,cnthv1a8), unit="log2")

# c, 2.5:1
sample4 <- as.integer((maxa/3.5)*2.5)
sample8 <- as.integer((maxa/3.5))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
ca = entropy(c(cnthv1a4,cnthv1a8), unit="log2")

# d, 2:1
sample4 <- as.integer((maxa/3)*2)
sample8 <- as.integer((maxa/3))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
da = entropy(c(cnthv1a4,cnthv1a8), unit="log2")

# e, 1.5:1
sample4 <- as.integer((maxa/2.5)*1.5)
sample8 <- as.integer((maxa/2.5))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
ea = entropy(c(cnthv1a4,cnthv1a8), unit="log2")

# f, 1:1
sample4 <- as.integer((maxa/2))
sample8 <- as.integer((maxa/2))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
fa = entropy(c(cnthv1a4,cnthv1a8), unit="log2")

# g, .5:1
sample4 <- as.integer((maxa/3))
sample8 <- as.integer((maxa/3)*2)
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
ga = entropy(c(cnthv1a4,cnthv1a8), unit="log2")

# g2 v2, .31:1
sample4 <- as.integer(maxa*(.31/1.31))
sample8 <- as.integer(maxa*(1/1.31))
smplhv1a4 <- hv1a4[sample(nrow(hv1a4),sample4),]
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),sample8),]
cnthv1a4 <- table(smplhv1a4)
cnthv1a8 <- table(smplhv1a8)
g2a <- entropy(c(cnthv1a4,cnthv1a8), unit="log2")


# h, 0:1
smplhv1a8 <- hv1a8[sample(nrow(hv1a8),maxa),]
cnthv1a8 <- table(smplhv1a8)
ha = entropy(cnthv1a8, unit="log2")

hv1ase <-c(aa,ba,ca,da,ea,fa,ga,g2a,ha)


#hv1b

# a, 1:0
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),maxb),]
cnthv1b4 <- table(smplhv1b4)
ab <- entropy(cnthv1b4, unit="log2")

# b, 3:1
sample4 <- as.integer((maxb/4)*3)
sample8 <- as.integer((maxb/4))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
bb <- entropy(c(cnthv1b4,cnthv1b8), unit="log2")

# c, 2.5:1
sample4 <- as.integer((maxb/3.5)*2.5)
sample8 <- as.integer((maxb/3.5))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
cb <- entropy(c(cnthv1b4,cnthv1b8), unit="log2")

# d, 2:1
sample4 <- as.integer((maxb/3)*2)
sample8 <- as.integer((maxb/3))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
db <- entropy(c(cnthv1b4,cnthv1b8), unit="log2")

# e, 1.5:1
sample4 <- as.integer((maxb/2.5)*1.5)
sample8 <- as.integer((maxb/2.5))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
eb <- entropy(c(cnthv1b4,cnthv1b8), unit="log2")

# f, 1:1
sample4 <- as.integer((maxb/2))
sample8 <- as.integer((maxb/2))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
fb <- entropy(c(cnthv1b4,cnthv1b8), unit="log2")

# g, .5:1
sample4 <- as.integer((maxb/3))
sample8 <- as.integer((maxb/3)*2)
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
gb <- entropy(c(cnthv1b4,cnthv1b8), unit="log2")

# g2 v2, .31:1
sample4 <- as.integer(maxb*(.31/1.31))
sample8 <- as.integer(maxb*(1/1.31))
smplhv1b4 <- hv1b4[sample(nrow(hv1b4),sample4),]
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),sample8),]
cnthv1b4 <- table(smplhv1b4)
cnthv1b8 <- table(smplhv1b8)
g2b <- entropy(c(cnthv1b4,cnthv1b8), unit="log2")

# h, 0:1
smplhv1b8 <- hv1b8[sample(nrow(hv1b8),maxb),]
cnthv1b8 <- table(smplhv1b8)
hb <- entropy(cnthv1b8, unit="log2")


hv1bse <-c(ab,bb,cb,db,eb,fb,gb,g2b,hb)

#d4a

# a, 1:0
smpld4a4 <- d4a4[sample(nrow(d4a4),maxa),]
cntd4a4 <- table(smpld4a4)
aa = entropy(cntd4a4, unit="log2")

# b, 3:1
sample4 <- as.integer((maxa/4)*3)
sample8 <- as.integer((maxa/4))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
ba = entropy(c(cntd4a4,cntd4a8), unit="log2")

# c, 2.5:1
sample4 <- as.integer((maxa/3.5)*2.5)
sample8 <- as.integer((maxa/3.5))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
ca = entropy(c(cntd4a4,cntd4a8), unit="log2")

# d, 2:1
sample4 <- as.integer((maxa/3)*2)
sample8 <- as.integer((maxa/3))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
da = entropy(c(cntd4a4,cntd4a8), unit="log2")

# e, 1.5:1
sample4 <- as.integer((maxa/2.5)*1.5)
sample8 <- as.integer((maxa/2.5))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
ea = entropy(c(cntd4a4,cntd4a8), unit="log2")

# f, 1:1
sample4 <- as.integer((maxa/2))
sample8 <- as.integer((maxa/2))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
fa = entropy(c(cntd4a4,cntd4a8), unit="log2")

# g, .5:1
sample4 <- as.integer((maxa/3))
sample8 <- as.integer((maxa/3)*2)
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
ga = entropy(c(cntd4a4,cntd4a8), unit="log2")

# g2 v2, .31:1
sample4 <- as.integer(maxa*(.31/1.31))
sample8 <- as.integer(maxa*(1/1.31))
smpld4a4 <- d4a4[sample(nrow(d4a4),sample4),]
smpld4a8 <- d4a8[sample(nrow(d4a8),sample8),]
cntd4a4 <- table(smpld4a4)
cntd4a8 <- table(smpld4a8)
g2a <- entropy(c(cntd4a4,cntd4a8), unit="log2")

# h, 0:1
smpld4a8 <- d4a8[sample(nrow(d4a8),maxa),]
cntd4a8 <- table(smpld4a8)
ha = entropy(cntd4a8, unit="log2")

d4ase <-c(aa,ba,ca,da,ea,fa,ga,g2a,ha)


#d4b

# a, 1:0
smpld4b4 <- d4b4[sample(nrow(d4b4),maxb),]
cntd4b4 <- table(smpld4b4)
ab <- entropy(cntd4b4, unit="log2")

# b, 3:1
sample4 <- as.integer((maxb/4)*3)
sample8 <- as.integer((maxb/4))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
bb <- entropy(c(cntd4b4,cntd4b8), unit="log2")

# c, 2.5:1
sample4 <- as.integer((maxb/3.5)*2.5)
sample8 <- as.integer((maxb/3.5))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
cb <- entropy(c(cntd4b4,cntd4b8), unit="log2")

# d, 2:1
sample4 <- as.integer((maxb/3)*2)
sample8 <- as.integer((maxb/3))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
db <- entropy(c(cntd4b4,cntd4b8), unit="log2")

# e, 1.5:1
sample4 <- as.integer((maxb/2.5)*1.5)
sample8 <- as.integer((maxb/2.5))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
eb <- entropy(c(cntd4b4,cntd4b8), unit="log2")

# f, 1:1
sample4 <- as.integer((maxb/2))
sample8 <- as.integer((maxb/2))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
fb <- entropy(c(cntd4b4,cntd4b8), unit="log2")

# g, .5:1
sample4 <- as.integer((maxb/3))
sample8 <- as.integer((maxb/3)*2)
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
gb <- entropy(c(cntd4b4,cntd4b8), unit="log2")

# g2 v2, .31:1
sample4 <- as.integer(maxb*(.31/1.31))
sample8 <- as.integer(maxb*(1/1.31))
smpld4b4 <- d4b4[sample(nrow(d4b4),sample4),]
smpld4b8 <- d4b8[sample(nrow(d4b8),sample8),]
cntd4b4 <- table(smpld4b4)
cntd4b8 <- table(smpld4b8)
g2b <- entropy(c(cntd4b4,cntd4b8), unit="log2")


# h, 0:1
smpld4b8 <- d4b8[sample(nrow(d4b8),maxb),]
cntd4b8 <- table(smpld4b8)
hb <- entropy(cntd4b8, unit="log2")

d4bse <-c(ab,bb,cb,db,eb,fb,gb,g2b,hb)

svg(paste(savepath, "CD4 CD8 TCR mixing Shannon.svg"), width=6, height=5, pointsize=15)

par(mar=c(4.1,4.1,1.1,.1))
plot(rep(0,8), ylim=c(10,14), xlim=c(.5,8.5), ylab="Shannon Entropy (bits)", xlab="Ratio CD4:CD8 TCRs", 
     axes=FALSE)#, main="D4/HV1")

# Plot mean HV/HIV values as bars in background - NB WILL REQUIRE UPDATING WHEN ALL IN

# hv mean + sds
hvm = 13.3026484813
hvsd = 0.5859745366

segments(.5,hvm,8.5,hvm, col=rgb(0,0,1,0.5), lwd=thick)
segments(.5,hvm+hvsd,8.5,hvm+hvsd, col=rgb(0,0,1,0.4), lwd=1)
segments(.5,hvm-hvsd,8.5,hvm-hvsd, col=rgb(0,0,1,0.4), lwd=1)


# hiv mean
hivm = 10.6447356723
hivsd = 0.651029893

segments(.5,hivm,8.5,hivm, col=rgb(1,0,0,0.5), lwd=thick)
segments(.5,hivm+hivsd,8.5,hivm+hivsd, col=rgb(1,0,0,0.4), lwd=1)
segments(.5,hivm-hivsd,8.5,hivm-hivsd, col=rgb(1,0,0,0.4), lwd=1)

# ratios
hivrm = 8-(2*0.310987069)
hivrsd = 0.1486751688
segments(hivrm,10,hivrm,14, lwd=2, col=rgb(1,0,0,0.4), lty=3)
segments(hivrm+hivrsd,10,hivrm+hivrsd,14, lwd=1, lty=3, col=rgb(1,0,0,0.3))
segments(hivrm-hivrsd,10,hivrm-hivrsd,14, lwd=1, lty=3, col=rgb(1,0,0,0.3))

pts <- c(1:7,hivrm,8)

points(pts, hv1ase, pch=3, lwd=med, col="green")
points(pts, hv1bse, pch=3, lwd=med, col="purple")

points(pts, d4ase, pch=4, lwd=med, col="green")
points(pts, d4bse, pch=4, lwd=med, col="purple")

lines(pts, hv1ase, lty=2, lwd=1, col=rgb(0,255,0,100, max=255))
lines(pts, hv1bse, lty=2, lwd=1, col=rgb(160,32,240,100, max=255))

lines(pts, d4ase, lty=2, lwd=1, col=rgb(0,255,0,100, max=255))
lines(pts, d4bse, lty=2, lwd=1, col=rgb(160,32,240,100, max=255))

legend("bottomleft", c("Alpha","Beta", "HV01", "HVD4", "Healthy ± sd", "HIV ± sd"), 
       pch=c(15,15,3,4,NA,NA), box.col="white",lty=c(NA,NA,NA,NA,1,1), 
       lwd=thick, col=c("green","purple", "black", "black", rgb(0,0,1,.5), rgb(1,0,0,.5)))

axis(2, at=c(10,11,12,13,14))
axis(1, at=c(1,2,3,4,5,6,7,hivrm,8), labels=ratios, las=1, cex.axis=.9)
axis(1, at=hivrm, labels="x̄")

box()

dev.off()
