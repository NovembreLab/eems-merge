require(data.table)
ref <- data.frame(fread("chip/GenomeWideSNP_6.na32.annot.csv", select=2:4))
bim <- read.table("raw/MARITIME_ROUTE.bim")

a1 <- as.factor(bim$V5)
a2 <- as.factor(bim$V6)
levels(a1) <- c('A', 'C', 'G')
levels(a2) <- c('A', 'C', 'G')

bim$V5 <- a1
bim$V6 <- a2

names(bim)[2] <- 'dbSNP.RS.ID'
bim <- cbind(bim, n=1:nrow(bim))


#ref <- chip[,c(2,3,4)]
ref <- ref[ref$Physical.Position != '---',]
ref <- ref[!duplicated(ref),]
ref$Physical.Position <- as.numeric(ref$Physical.Position)

q <- merge(bim, ref, all.x=T, sort=F)  

ERROR = is.na(q$Chromosome)

write.table(q$dbSNP.RS.ID[ERROR],
	    'tmp/PASCHOU_exclude.txt', quote=F,
	    row.names=F, col.names=F)
q$Chromosome[ERROR] <- q$V1[ERROR]
q$Physical.Position[ERROR] <- q$V4[ERROR]    

#system('cp raw/MARITIME_ROUTE.bed data/MARITIME_ROUTE.bed')
#system('cp raw/MARITIME_ROUTE.fam data/MARITIME_ROUTE.fam')
q <- q[order(q$n),c(8,1,3,9,5,6)]
write.table(q, 'tmp/tmpPASCHOU.bim', quote=F,
	    row.names=F, col.names=F)

s = paste('plink --bfile raw/MARITIME_ROUTE --make-bed',
	'--out data/MARITIME_ROUTE ',
	'--exclude tmp/PASCHOU_exclude.txt',
	'--bim tmp/tmpPASCHOU.bim')
system(s)


