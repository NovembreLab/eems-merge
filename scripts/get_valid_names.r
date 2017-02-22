a <- read.table("meta/individuals2.txt", header=T, fill=T, strings=F)
b <- unique(a[,1])
write.table(cbind(b,b), "input_files/valid_names.txt", col.names=F, row.names=F,
            quote=F)

