library(dplyr)
x <- readxl::read_excel("sources/Data_for_Ben_Meta.xlsx")[,c(1,14:15)]
names(x) <- c('originalId', 'source', 'chip')
x$wasDerivedFrom <- 'Estonians'
y <- read.csv("intermediate/gvar.indiv_prov", as.is=T)
q <- left_join(y, x)
is.est <- q$wasDerivedFrom=='Estonians'
q$wasDerivedFrom[is.est] <- q$source[is.est]
q %>% select(-chip, -source) %>% 
    write.csv('pgs/gvar.indiv_prov', row.names=F)
