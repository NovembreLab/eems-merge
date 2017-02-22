suppressPackageStartupMessages({
library(scales)
library(RColorBrewer)
library(dplyr)
})
get_cols <- function(axis1, axis2, cmat=color_mat()){
    x = floor(100 * (axis1 - min(axis1))/diff(range(axis1)) )+1
    y = floor(100 * (axis2 - min(axis2))/diff(range(axis2)) )+1
    cv <- cmat[cbind(x, y)]
}

color_mat <- function(col_top=c("black", "darkblue", "darkgreen", 'yellow'),
                      col_bot=c("gray80", "pink", 'red', 'orange')
                      ){
    grad_top <- rev(gradient_n_pal(col_top)(0:100/100))
    grad_bot <- rev(gradient_n_pal(col_bot)(0:100/100))
    cols <- t(mapply(function(x,y)gradient_n_pal(c(x,y))(0:100/100),
                    grad_bot, grad_top))
}
get_cols_wrap <-function(pop_geo, cmat=color_mat()){
    to_wrap <- pop_geo$longitude < -34
    pop_geo$longitude[to_wrap] <-pop_geo$longitude[to_wrap]  + 360
    cols <- get_cols(pop_geo$longitude, pop_geo$latitude)
    cols <- get_cols(pmin(pop_geo$longitude, 180), pop_geo$latitude)
    to_wrap <- pop_geo$longitude > 180
    pop_geo$longitude[to_wrap] <-pop_geo$longitude[to_wrap]  - 360
    print(length(cols))
    cols
}



p <- read.csv("pgs/gvar.pop_geo") %>% arrange(popId)
d <- read.csv("pgs/gvar.pop_display")%>% arrange(popId)
cols <- get_cols_wrap(p, color_mat())
d$color <- cols

d$name <- as.character(d$name)

write.csv(d, "pgs/gvar.pop_display", row.names=F)

#Abbrev + Language Family
#'PF',  'Polynesian', 'smo'
#'US-Ale', 'Eskimo-Aleut', 'ale'
#'CA-Hai', 'Na-Dene', 'or isolate', 'hai'
#'CA-Tsi', 'Penutian', 'tsi'
#'CA-Ni', 'Penutian', 'ncg'
#'CA-Ath', 'Na-Dene', 'ath'
#'CA-Sts', 'Salishan', 'shs'
#'CA-Spl', 'Salishan', 'shs'
#'MX-Pim', 'Uto-Aztecan', 'ood'
#'CA-Chi', 'Na-Dene','chp'
#'CA-Cre', 'Algic', 'cre'
#'MX-Mxt', 'Oto-Manguean', 'miz'
#'MX-Mxe', 'Mixe-Zoque', 'mxq'
#'MX-Zap', 'Oto-Manguean', 'zap'



