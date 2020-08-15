library(dplyr)
library(magrittr)

a = readRDS("gmap.RDS")
# cm/bp
d = mean(unlist(lapply(a, function(x){ tail(x$map,n=1)/tail(x$ppos,n=1)})))

print(d)
top = ( 1 - exp(-d/50) )
bot = 2 
print(top)
print(top/bot)
