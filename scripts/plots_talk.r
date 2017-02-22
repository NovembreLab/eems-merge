require(rworldmap)

WRAP_PT = -31

plot_worldmap <- function(...){
    m <- getMap('low')
    for(i in 1:244){
        x <- m@polygons[[i]]@Polygons
        for(j in 1:length(x)){
            k <- x[[j]]@coords 
            r <- range(k[,1])
            if(r[[1]] < WRAP_PT && r[[2]] > WRAP_PT){
                k[,1] <- pmin(k[,1]+360, -WRAP_PT+360)
                print(c("gaga", i, j))
            }
            if(r[[2]] < WRAP_PT){
                k[,1] <- k[,1] + 360
            }
            m@polygons[[i]]@Polygons[[j]]@coords  <- k
        }
    }

    plot(NA, xlim=c(-31, -31+360), ylim=c(-50,75), xaxs='i', 
         xlab='', ylab='', xaxt='n', yaxt='n', ...) 
    plot(m, add=T, ...)                                           
}


wrap.america <- function(x){
    x$LONG[x$LONG < WRAP_PT] <- x$LONG[x$LONG < WRAP_PT] + 360
    x
}


fig1.1 <- function(){
    pdf("fig1.1.pdf", width=11, height=6)
    plot_worldmap(col='lightgray', asp=1)   
    a <- wrap.america(read.csv("meta/used_pops.csv"))               
    points(a$LONG, a$LAT, col="black", pch=16)  
    dev.off()
}

fig1.1 <- function(name='human_origins', of='ho'){
    #setEPS()
    #postscript(sprintf("fig1.%s.eps", of), width=8, height=4)
    png(sprintf("fig1.%s.png", of), width=1100, height=500)
    par(mar=c(0,0,0,0), bty='n')
    plot_worldmap(col='#eeeeee', asp=1, axes=0)   
    a <- wrap.america(read.csv("meta/used_pops.csv"))               
    others = a$DATA %in% name
    cv <- rep('black', nrow(a))
    cx <- rep(1, nrow(a))
    cv[others] <- '#7f0000'
    cv[others] <- '#ffb000'
    cx[others] <- 1.3
    points(a$LONG, a$LAT, col=cv, pch=16, cex=cx)  
    dev.off()
}

fig1 <- function(){
fig1.1(of='eb', name=levels(a$DATA)[c(1,2,3,6,7,12,13)]) 
fig1.1(of='ho', name='human_origins') 
fig1.1(of='pr', name='popres') 
fig1.1(of='pa', name='panam') 
fig1.1(of='xi', name='xing2010') 
fig1.1(of='he', name='henn2012') 
fig1.1(of='re', name='reich2011') 
fig1.1(of='1', name='1')
}
