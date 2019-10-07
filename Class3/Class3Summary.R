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