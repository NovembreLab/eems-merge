a <- read.table("meta/individuals3.txt", fill=T, header=T, strings=F)       
b <- read.csv("meta/pop2.csv")                                              
names(b)[1] <- 'POPULATION'                                                 
q <- merge(a, b)                                                            
df <- data.frame(POP=q$POPULATION, LAT=q$lat, LONG=q$long, DATA=q$dataset)  
write.csv(file="meta/used_pops.csv", unique(df), row.names=F, quote=F)      

