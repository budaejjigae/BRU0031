library(rgdal)
library(spdep)

setwd("~/Downloads/BRU0031/columbus")

columbus = readOGR(dsn = ".", layer = "columbus")
columbus_queen = read.gal(file = "columbus.gal", region.id = columbus$POLYID)
plot(columbus)
plot(columbus_queen, coordinates(columbus), add = TRUE)

W = nb2listw(neighbours = columbus_queen)
W = listw2mat(W)
View(W)