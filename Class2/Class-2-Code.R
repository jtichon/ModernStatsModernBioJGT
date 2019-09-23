library(here)

#1.2
dpois(x=0:12, lambda=5)

png((here("Class2", "poisonbarplot.png")))
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
