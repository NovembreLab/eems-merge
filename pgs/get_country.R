require(rworldmap)
require(SDMTools)
m <- getMap("high")
pop_geo <- read.csv("pgs/gvar.pop_geo")
pop_display <- read.csv("pgs/gvar.pop_display")

get_coords <- function(m, i){
    print(length(m@polygons[[i]]@Polygons))
    lapply(
	m@polygons[[i]]@Polygons,
	function(p)p@coords)
}

get_pops_in_country <- function(m, pop_geo, i){
    crds <- get_coords(m, i)
    x <- lapply(crds, function(crd)
    which(pnt.in.poly(pop_geo[,3:2], crd)[,3]==1)
    )
    do.call(c, x)
}

get_relevant_pops <- function(pop_display, m, pop_geo, i){
    pop_no <- get_pops_in_country(m, pop_geo, i)
    pops <- pop_display[pop_no,]
    print(c(as.character(m$NAME[i]),
	    as.character(m$ISO_A2[i])))
    print(pops)
}
