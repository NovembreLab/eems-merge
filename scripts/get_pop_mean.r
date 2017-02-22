#' simple script to get mean from pops with multiple coords, particularly in
#' HUGO data

a <- read.csv("meta/pop.csv")                       
pop2 <- aggregate(cbind(a$lat,a$long,a$admixed), list(a$POP), mean)
names(pop2) <- names(a)
write.csv(pop2, "meta/pop2.csv", quote=F, row.names=F)

