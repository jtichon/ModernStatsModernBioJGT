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
