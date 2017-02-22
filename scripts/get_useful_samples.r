a <- read.table("meta/individuals.txt", fill=T, strings=F, header=T)
a <- a[a$dataset != '1000g',]

q <- read.csv("meta/pop2.csv")
q$admixed[is.na(q$admixed)] = 0
q2 <- q[which(q$admixed==0 & !is.na(q$lat)),]
pops <- as.character(q2[,1])

l <- aggregate(a$ID, list(a$POPULATION), length)

present_pops <- l[l[,1] %in% pops,]$Group.1

inds2 <- a[a$POPULATION %in% present_pops,] 


pop_per_ind <- split(inds2$POPULATION, inds2$ID)

dups <- 1:2; 
for(i in pop_per_ind) 
    if(length(i)>1){
        if (length(i) != 2){
            print(length(i));
            print(i)
        }
        dups<-rbind(dups,i)
    }

dups <- unique(dups)

for(i in 1:nrow(dups)){
    d <- dups[i,]
    print(d[2])
    inds2$POPULATION[inds2$POPULATION == d[1]] <- d[2] 
}



write.table(inds2, "meta/individuals2.txt", row.names=F, quote=F)

