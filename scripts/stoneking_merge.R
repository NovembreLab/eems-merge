a <- read.csv("regions/location_coords.csv")           
b <- read.csv("sources/Stoneking.pops.txt", sep="\t")  
x <- merge(b, a, by.x=c("pop.ID", "sampling.location"),
           by.y=c("Population", "Group_Population"))

x <- x[,c('sample.ID', 'ID')]
names(x) <- c('SID', 'PID')
x <- data.frame(x, Source='Stoneking', Permission='Private_Demo')

write.csv(x, "intermediate/sample_pop_stoneking.csv",
          row.names=F, quote=F)


