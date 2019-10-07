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
<<<<<<< HEAD
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
=======
dev.off()
>>>>>>> 054bff047acba9186fece3916c2a0396b8aeea7b
