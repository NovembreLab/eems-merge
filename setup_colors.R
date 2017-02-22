require(rworldmap)
library(grid)
library(scales)
library(sp)
library(maps)
library(rgeos)
library(maptools)

pop_geo <- read.csv("pgs/gvar.pop_geo")


pop_display <- read.csv("pgs/gvar.pop_display")


coords2country = function(points)
{  
      countriesSP <- getMap(resolution='high')
  #countriesSP <- getMap(resolution='high') #you could use high res map from
#rworldxtra if you were concerned about detail

  # convert our list of points to a SpatialPoints object

  # pointsSP = SpatialPoints(points, proj4string=CRS(" +proj=longlat
# +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

  #setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  


    # use 'over' to get indices of the Polygons object containing each point 
    indices = over(pointsSP, countriesSP)

    # return the ADMIN names of each country
    #indices$ADMIN  
    indices$continent
      #indices$ISO3 # returns the ISO3 code 
      #indices$continent   # returns the continent (6 continent model)
      #indices$REGION   # returns the continent (7 continent model)
}
library(scales)

four.color.matrix <-
    function( mycols ){

        m <- matrix( NA , 100 , 100 )

        m[ 1 , 1 ] <- mycols[ 1 ] 
        m[ 1 , 100 ] <- mycols[ 2 ]
        m[ 100 , 1 ] <- mycols[ 3 ]
        m[ 100 , 100 ] <- mycols[ 4 ]

        m[ 1 , 1:100 ] <- gradient_n_pal( c( mycols[ 1 ] , 'black' , mycols[ 2 ] ) , values = c( 1 , 50 , 100 ) )(1:100)
        m[ 1:100 , 1 ] <- gradient_n_pal( c( mycols[ 1 ] , 'black' , mycols[ 3 ] ) , values = c( 1 , 50 , 100 ) )(1:100)
        m[ 1:100 , 100 ] <- gradient_n_pal( c( mycols[ 2 ] , 'black' , mycols[ 4 ] ) , values = c( 1 , 50 , 100 ) )(1:100)
        m[ 100 , 1:100 ] <- gradient_n_pal( c( mycols[ 3 ] , 'black' , mycols[ 4 ] ) , values = c( 1 , 50 , 100 ) )(1:100)

        a <- gradient_n_pal( c( mycols[ 1 ] , 'black' , mycols[ 4 ] ) , values = c( 1 , 50 , 100 ) )
        diag(m)<-a(1:100)

        b <- gradient_n_pal( c( mycols[ 3 ] , 'black' , mycols[ 2 ] ) , values = c( 1 , 50 , 100 ) )
        for(i in 1:(nrow(m) - 1)){ 
          for (j in 1:nrow(m)) if (i + j == nrow( m )+1){
              m[i,j] <- b(j)
            }
        }

        for ( i in 2:50 ){

            m[ i , i:(101-i) ] <- 
                gradient_n_pal( c( mycols[ 1 ] , 'black' , mycols[ 2 ] ) , values = c( 0 , 50 , 100 ) )(  i:(101-i) )

            m[ i:(101-i) , i ] <- 
                gradient_n_pal( c( mycols[ 3 ] , 'black' , mycols[ 1 ] ) , values = c( 0 , 50 , 100 ) )( (101-i):i )

        }



        for ( i in 51:99 ){

            m[ i , i:(101-i) ] <- 
                gradient_n_pal( c( mycols[ 3 ] , 'black' , mycols[ 4 ] ) , values = c( 0 , 50 , 100 ) )(  i:(101-i) )

            m[ i:(101-i) , i ] <- 
                gradient_n_pal( c( mycols[ 4 ] , 'black' , mycols[ 2 ] ) , values = c( 0 , 50 , 100 ) )( (101-i):i )

        }

        m
    }



plot_map <- function(){
    to_wrap <- pop_geo$longitude < -34
    pop_geo$longitude[to_wrap] <-pop_geo$longitude[to_wrap]  + 360

#    long = 100 * (pop_geo$longitude - min(pop_geo$longitude))/diff(range(pop_geo$longitude)) 
#    lat = 100 * (pop_geo$latitude - min(pop_geo$latitude))/diff(range(pop_geo$latitude)) 
#    cols <- four.color.matrix(c('red','blue','green', 'purple'))
#    cv <- cols[cbind(long, lat)]
    

grad_top <- rev(gradient_n_pal(c("gray80", "purple", "darkblue", "darkgreen", 'yellow'))(0:100/100))
grad_bot <- rev(gradient_n_pal(c("black", "pink", 'red', 'orange'))(0:100/100))
require(RColorBrewer)
#grad_top <- gradient_n_pal(brewer.pal(8, "Dark2")[c(8,7,2,6,5,1,3,4)])(0:100/100)
#grad_bot <- gradient_n_pal(brewer.pal(8, "Pastel2")[c(8,7,2,6,5,1,3,4)])(0:100/100)
cols = t(mapply(function(x,y)gradient_n_pal(c(x,y))(0:100/100), grad_bot, grad_top))
    long = floor(100 * (pop_geo$longitude - min(pop_geo$longitude))/diff(range(pop_geo$longitude)) )+1
    lat = floor(100 * (pop_geo$latitude - min(pop_geo$latitude))/diff(range(pop_geo$latitude)) ) + 1
    cv <- cols[cbind(long, lat)]
    map <- getMap('low') 
    to_wrap <- pop_geo$longitude > 180
    pop_geo$longitude[to_wrap] <-pop_geo$longitude[to_wrap]  - 360

    pdf(file="proposed_colors.pdf", width=11, height=8)
    plot(map)
    points(pop_geo$longitude, pop_geo$latitude, col=cv, lwd=2, pch=16, cex=2) 
    dev.off()
}







