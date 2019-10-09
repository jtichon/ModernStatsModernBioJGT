library(here)
###2.3
load(here("data", "e100.RData"))
#remove outlier
e99 = e100[-which.max(e100)]
#see picture of distribution
png((here("Class3", "e99barplot.png")))
barplot(table(e99), space = 0.8, col = "chartreuse")
dev.off()

library(vcd)
gf1 = goodfit( e99, "poisson")
png((here("Class3", "poissonrootogram.png")))
rootogram(gf1, xlab = "", rect_gp = gpar(fill = "chartreuse"))
dev.off()

##Q2.1
pois.100<-rpois(100,0.5)
png((here("Class3", "randompoissonrootogram.png")))
gf2 = goodfit(pois.100, "poisson")
rootogram(gf2, xlab="", rect_gp = gpar(fill = "chartreuse"))
dev.off()

table(e100)
table(rpois(100,3))

## Q2.3


#Function for trying different m's
pois100prob<-function(m){
  prod(dpois(c(0,1,2,7),lambda = m)^(c(58, 34, 7, 1)))
}

pois100prob(0)
pois100prob(1)
pois100prob(2)
pois100prob(0.4)

loglikelihood = function(lambda, data = e100){
  sum(log(dpois(data,lambda)))
}

png((here("Class3", "loglikpoisson.png")))
lambdas = seq(0.05, 0.95, length = 100)
loglik = vapply(lambdas, loglikelihood, numeric(1))
plot(lambdas, loglik, type = "l", col = "red", ylab = "", lwd = 2, xlab = expression(lambda))
m0 = mean(e100)
abline(v = m0, col = "blue", lwd = 2)
abline(h = loglikelihood(m0), col = "purple", lwd = 2)
m0
dev.off()

gf = goodfit(e100, "poisson")
names(gf)
gf$par

png((here("Class3", "rootogram55.png")))
pois.100<-rpois(100,0.55)
gf2 = goodfit(pois.100, "poisson")
rootogram(gf2, xlab="", rect_gp = gpar(fill = "chartreuse"))
dev.off()

## 2.4

#Make up matching dataset
cb<-c(rep(0,110), rep(1,10))
table(cb)

probs = seq(0, 0.3, by = 0.005)
likelihood = dbinom(sum(cb), prob = probs, size = length(cb))
png((here("Class3", "likelihoodbinom.png")))
plot(probs, likelihood, pch = 16, xlab = "probability of success", ylab = "likelihood", cex=0.6)
dev.off
probs[which.max(likelihood)]

#loglikelihood binomial

loglikelihood.binom = function(theta, n = 300, k = 40){
  115 + k * log(theta) + (n - k) * log(1 - theta)
}

thetas = seq(0, 1, by = 0.001)
png((here("Class3", "loglikelihoodbinom.png")))
plot(thetas, loglikelihood.binom(thetas), xlab = expression(theta), 
     ylab = expression(paste("log f(", theta, " | y")), type = "l")
dev.off()

#2.5
library("Biostrings")
here("data", "e100.RData")
staph = readDNAStringSet(here("data","staphsequence.ffn.txt"), "fasta")
staph[1]
letterFrequency(staph[[1]], letters = "ACGT", OR = 0)


#Find letter frequency
letterFrq = vapply(staph, letterFrequency, FUN.VALUE = numeric(4),
                   letters = "ACGT", OR = 0)
colnames(letterFrq) = paste0("gene", seq(along = staph))
#Compute frequencies in first 10 genes and convert to proportions
tab10 = letterFrq[, 1:10]
computeProportions = function(x) { x/sum(x) }
prop10 = apply(tab10, 2, computeProportions)
round(prop10, digits = 2)
p0 = rowMeans(prop10)
p0

cs = colSums(tab10)
cs
expectedtab10 = outer(p0, cs, FUN = "*")
round(expectedtab10)

randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) } )
all(colSums(randomtab10) == cs)

stat = function(obsvd, exptd = 20 * pvec) {
  sum((obsvd - exptd)^2 / exptd)
}
B = 1000
simulstat = replicate(B, {
  randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) })
  stat(randomtab10, expectedtab10)
})
S1 = stat(tab10, expectedtab10)
sum(simulstat >= S1)

png((here("Class3", "multinomsamplingdistrn.png")))
hist(simulstat, col = "lavender", breaks = seq(0, 75, length.out=50))
abline(v = S1, col = "red")
abline(v = quantile(simulstat, probs = c(0.95, 0.99)),
       col = c("darkgreen", "blue"), lty = 2)
dev.off()
#2.6 

#2.6.1

##Q2.10
#Compare quantiles of actual chi-squared vs our random generated test 
#statistics for multinomial p0
qs = ppoints(100)
png((here("Class3", "quantilecompare.png")))
par(mfrow=c(1,2))
hist(quantile(simulstat, qs), main = "Test stats under H0 Multinomial")
hist(quantile(qchisq(qs, df = 30), qs), main = "Randomly generated chi-squared")
dev.off()

png((here("Class3", "qqplote100vschi.png")))
qqplot(quantile(simulstat, qs),quantile(qchisq(qs, df = 30), qs))
dev.off()

png((here("Class3", "textqqplote100vschi.png")))
qqplot(qchisq(ppoints(B), df = 30), simulstat, main = "",
       xlab = expression(chi[nu==30]^2), asp = 1, cex = 0.5, pch = 16)
abline(a = 0, b = 1, col = "red")
dev.off()

#2.7

load(here("data", "ChargaffTable.RData"))
ChargaffTable

png((here("Class3", "chargaffbarplots.png")))
par(mfrow=c(2,4))
barplot(ChargaffTable[1,])
barplot(ChargaffTable[2,])
barplot(ChargaffTable[3,])
barplot(ChargaffTable[4,])
barplot(ChargaffTable[5,])
barplot(ChargaffTable[6,])
barplot(ChargaffTable[7,])
barplot(ChargaffTable[8,])

#Calculate proposed test stat
statChf = function(x){
  sum((x[, "C"] - x[, "G"])^2 + (x[, "A"] - x[, "T"])^2)
}
chfstat = statChf(ChargaffTable)

#Calculate test stats under a permutation test of relabeling C/G and A/T
permstat = replicate(100000, {
  permuted = t(apply(ChargaffTable, 1, sample))
  colnames(permuted) = colnames(ChargaffTable)
  statChf(permuted)
})

#p-value
pChf = mean(permstat <= chfstat)
pChf

hist(permstat, breaks = 100, main = "", col = "lavender")
abline(v = chfstat, lwd = 2, col = "red")

# 7.3.1

load(here("Data", "Deuteranopia.RData"))
Deuteranopia
#Chi-square Test for Independence
chisq.test(Deuteranopia)


library("HardyWeinberg")
data("Mourant")
Mourant[214:216,]

nMM = Mourant$MM[216]
nMN = Mourant$MN[216]
nNN = Mourant$NN[216]
loglik = function(p, q = 1 - p) {
  2 * nMM * log(p) + nMN * log(2*p*q) + 2 * nNN * log(q)
}
xv = seq(0.01, 0.99, by = 0.01)
yv = loglik(xv)
png(here("Class3","loglikelihoodtahiti.png"))
plot(x = xv, y = yv, type = "l", lwd = 2,
     xlab = "p", ylab = "log-likelihood")
imax = which.max(yv)
abline(v = xv[imax], h = yv[imax], lwd = 1.5, col = "blue")
abline(h = yv[imax], lwd = 1.5, col = "purple")
dev.off()

#af is a function from HardyWeinberg package to calculate predicted proportions
phat  =  af(c(nMM, nMN, nNN))
phat
pMM   =  phat^2
qhat  =  1 - phat
pHW = c(MM = phat^2, MN = 2*phat*qhat, NN = qhat^2)
sum(c(nMM, nMN, nNN)) * pHW

pops = c(1, 69, 128, 148, 192)
genotypeFrequencies = as.matrix(Mourant[, c("MM", "MN", "NN")])
png(here("Class3", "TernaryPlot.png"))
HWTernaryPlot(genotypeFrequencies[pops, ],
              markerlab = Mourant$Country[pops],
              alpha = 0.0001, curvecols = c("red", rep("purple", 4)),
              mcex = 0.75, vertex.cex = 1)
dev.off()
