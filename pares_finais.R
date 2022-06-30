rm(list=ls(all=T))

#descomentar para instalar
#install.packages("dplyr")

library(dplyr)    

#colocar a diretoria certa
#setwd("~/cancro")
setwd("~/Downloads/cancro")

pares = read.csv(file=paste("output/1-output of fractions of gene expression patterns.csv",sep=""))

# eliminar p values acima (ou = ) de 0.05
pares = pares %>% filter(p.vaule < 0.05)

#os 10 primeiros pares com maior valor de up up
indices_pares_up_up = order(pares$up.up, decreasing = TRUE)[1:10]
pares_up_up = pares[indices_pares_up_up,]

#os 10 primeiros pares com maior valor de down down
indices_pares_down_down = order(pares$down.down, decreasing = TRUE)[1:10]
pares_down_down = pares[indices_pares_down_down,]

View(pares_down_down)
View(pares_up_up)