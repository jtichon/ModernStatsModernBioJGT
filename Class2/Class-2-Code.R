library(here)

#1.2
dpois(x=0:12, lambda=5)

png((here("Class2", "poissonbarplot.png")))
barplot(dpois(0:12,5), names.arg = 0:12, col = "red")
dev.off()

#1.3
genotype = c("AA", "AO", "BB", "AO", "OO", "AO", "AA", "BO", "BO", 
             "AO", "BB", "AO", "BO", "AB", "OO", "AB", "BB", "AO", "AO")
table(genotype)

#Access factor levels
genotypeF = factor(genotype)
levels(genotypeF)
table(genotypeF)

## 1.3.1
rbinom(15, prob= 0.5, size = 1)
rbinom(12, prob = 2/3, size = 1)
rbinom(1, prob = 2/3, size = 12)

set.seed(235569515)
rbinom(1, prob = 0.3, size = 15)

#Q 1.3
for(j in 1:10)
{
  print(rbinom(1, prob = 0.3, size = 15))
}

probabilities = dbinom(0:15, prob = 0.3, size = 15)
round(probabilities, 2)

png((here("Class2", "binomialbarplot.png")))
barplot(probabilities, names.arg = 0:15, col ="red")
dev.off()

#Q 1.4
dbinom(3, size = 4, prob = 2/3)
choose(4,3)*(2/3)^3*(1-2/3)^(4-3)


##1.3.2

#Q5 
dpois(3,5)
5^3*exp(-5)/factorial(3)

simulations = rbinom(n = 30000, prob = 5e-4, size = 10000)
png((here("Class2", "poissonbarplot2.png")))
barplot(table(simulations), col = "lavender")
dev.off()

##1.3.4
load(here("data", "e100.RData"))

png(here("Class2", "barplote100.png"))
barplot(e100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5, 100.5),
        names.arg = seq(along = e100), col = "darkolivegreen")
dev.off()

1-ppois(6,0.5)
ppois(6, 0.5, lower.tail=FALSE)

1-ppois(6,0.5)^100

#Code for simulating this probability
maxes = replicate(100000, {
  max(rpois(100,0.5))
})
table(maxes)

## Can't replicate text's 0.00016
#This is the average of a logical vector of TRUE/FAlSE for is >=7 or not
mean( maxes >= 7 )

#1.4

rmultinom(1, size = 4, prob=c((1/8),(3/8),(3/8),(1/8)))

#q1.7
pvec = rep(1/4, 4)
dmultinom(c(4,2,0,0), size=6, prob = pvec)
t(rmultinom(1, prob = pvec, size = 8))

#q1.8 t() stands for transpose
rmultinom(1, prob = pvec, size = 8)

#q1.9 
rmultinom(n = 8, prob = pvec, size = 1) #8 samples of size 1
rmultinom(n = 1, prob = pvec, size = 8) #1 sample of size 8

#1.4.1

#RAndomly generated under null hypothesis, all are equally likely
obsunder0 = rmultinom(1000, prob = pvec, size = 20)
dim(obsunder0)
obsunder0[, 1:11]

#Expected if all 0
expected0 = pvec * 20

#Function to create test stat - chi-squared
stat = function(obsvd, exptd = 20 * pvec) {
  sum((obsvd - exptd)^2 / exptd)
}

#Give sum contribution of first observation
stat(obsunder0[, 1])

#Apply the test stat to each column element of obsunder0
#apply second option is (1,2) for (row,column)
S0 = apply(obsunder0, 2, stat)
summary(S0)

#Quantiles for distribution of deviations (contributions to test stat)
q95 = quantile(S0, probs = 0.95)
q95

q05 = quantile(S0, probs = 0.05)
q05

#Reject if test stat is greater than 7.6 (i.e. if weighted mean is larger than
#what would be the 95th percentiles of contributions under H0)

#Find power for 3/8, 1/4, 3/12, 1/8 distribution by simulation. (i.e.
#probability of correctly reject H0)

pvecA = c(3/8, 1/4, 3/12, 1/8)
observed = rmultinom(1000, prob = pvecA, size = 20)
dim(observed)

#display first seven, calculate contribution of first column
observed[, 1:7]
apply(observed, 1, mean)

#calculate expected counts under new distribution
expectedA = pvecA * 20
expectedA

stat(observed[, 1])
S1 = apply(observed, 2, stat)
q95
sum(S1 > q95)

power = mean(S1 > q95)
power

#1.1 
#Geometric rgeom(n, p)
#Hypergeometric rhyper(N,A,B,k)

#1.2
dbinom(0,10,.3)+dbinom(1,10,.3)+dbinom(2,10,.3)
pbinom(2,10,.3)

#1.3

#Function to find probability the max of n poisson samples with mean lambda is
#equal to m
max.pois<-function(m,n,lambda)
{
  #Probability max is less than m (i.e. all n values less than m)
  prob.lessm<- ppois(m-1,lambda)^n
  
  #Probability max is less than m+1 (i.e. all n values less than m+1)
  prob.lessm1<- ppois(m,lambda)^n
  
  #Probability max is equal to m (P(max<m+1)-P(max<m))
  probm<-prob.lessm1-prob.lessm
}

asbigmax.pois<-function(m,n,lambda){
  #Find probability of all maxes less than m
  max.prob<-c(0:m-1)
  for(j in 0:m-1)
  {
    max.prob[j]<-max.pois(j,n,lambda)
  }
  
  #Find probability of max at least as large as m (i.e. 1 - prob less than m)
  prob<-1-sum(max.prob)
  prob
}

asbigmax.pois(9,20,4)

#Take 2:

maxormore<-function(m,n,lambda)
{
  1-ppois(m-1,lambda)^n
}

maxormore(9,20,4)


#1.4

#Function to find probability the max of n poisson samples with mean lambda is
#equal to m
max.pois<-function(m,n,lambda)
{
  #Probability max is less than m (i.e. all n values less than m)
  prob.lessm<- ppois(m-1,lambda)^n
  
  #Probability max is less than m+1 (i.e. all n values less than m+1)
  prob.lessm1<- ppois(m,lambda)^n
  
  #Probability max is equal to m (P(max<m+1)-P(max<m))
  probm<-prob.lessm1-prob.lessm
}

asbigmaxd.pois<-function(m=10,n=20,lambda=5){
  #Find probability of all maxes less than m
  max.prob<-c(0:m-1)
  for(j in 0:m-1)
  {
    max.prob[j]<-max.pois(j,n,lambda)
  }
  
  #Find probability of max at least as large as m (i.e. 1 - prob less than m)
  prob<-1-sum(max.prob)
  prob
}

asbigmaxd.pois()

#Take 2:

maxormored<-function(m=10,n=20,lambda=5)
{
  1-ppois(m-1,lambda)^n
}

maxormored()

#1.5

#Not sure what this means. 100 trials mean 100 samples or 100 positions
#to test? What do they mean by "number of simulations to prove"?

#Not sure what this means. 

#1.6

?Distributions

#Not discrete: chi-squared, exponential, F, gamma, normal

#Geometric
dgeom(x=0:12, p=.6)

png((here("Class2", "geometricbarplot.png")))
barplot(dgeom(0:12,p=0.6), names.arg = 0:12, col = "red")
dev.off()

#Standard Normal
dnorm(x=seq(-3,3, by=0.1)) 

png((here("Class2", "geometricbarplot.png")))
barplot(dnorm(x=seq(-3,3, by=0.5)), names.arg = seq(-3,3, by=0.5), col = "red")
dev.off()


#1.7
x.pois<-rpois(100,3)
mean(x.pois)
var(x.pois)

#1.8

## Let's leave this for together








