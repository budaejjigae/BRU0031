library(lattice)

url <- "https://raw.githubusercontent.com/steviep42/youtube/master/YOUTUBE.DIR/wines.csv"
my.wines <- read.csv(url, header=TRUE)

my.prc <- prcomp(my.wines[,-1], center=TRUE, scale=TRUE)

summary(my.prc)
screeplot(my.prc, main="Scree Plot", xlab="Components")
screeplot(my.prc, main="Scree Plot", type="line" )
abline(1, 0, lty = 2)

biplot(my.prc, cex=c(1, 0.7))