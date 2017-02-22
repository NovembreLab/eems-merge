a <- read.csv("intermediate/locations_all.csv")
u <- a$Uncertainty
u[is.na(u)] <- 0
a$Uncertainty <- u

pgswriter <- function(...){
    write.csv(...,row.names=F)
}

numeric_id <- function(df){
    df %>% mutate(popId=as.numeric(as.factor(popId))+1000)
}

if(F){
palette(c("black", topo.colors(6)))
plot(a$Longitude, a$Latitude,
     col=a$Certainty.group, cex=pmax(1, a$Uncertainty/100),
       pch=16, 
       xlim=c(40,160),
       ylim=c(10,80))                                      
plot(b$Longitude, b$Latitude,
     col=b$Certainty.group, cex=pmax(1, b$Uncertainty/100),
       pch=16, 
       xlim=c(80,130),
       ylim=c(50,80))                                      
require(maps)
map(add=T)
}


require(dplyr)
tib_geo <- read.csv("tib/HGDP_Tibetan_Merged_160509.pop_geo", as.is=T)
pops.to.keep <- c(
    "Tsum",
    "UM",
    "TBB_Nachu",
    "TBB_Nyingchi",
    "TBB_Shannan",
    "TBX_Chamdo",
    "TBX_Lhasa",
    "TBX_Shannon",
    "TBX_Shigatse")
tib_ind <- read.csv("tib/HGDP_Tibetan_Merged_160509.indiv_label", as.is=T)
tib_prov <- read.csv("tib/HGDP_Tibetan_Merged_160509.indiv_prov", as.is=T)

tib <- tib_ind %>% left_join(tib_geo) %>% 
    left_join(tib_prov) %>% 
    filter(popId %in% pops.to.keep)

tib %>% write.csv("tib/tibetan.csv", quote=F) 

tib[,c(2,1)] %>% write.table("tib/tib.plink", quote=F, row.names=F)

tib <- numeric_id(tib)

tib %>% select(sampleId, wasDerivedFrom, used, 
               originalId, permissions) %>% 
    pgswriter("tib/tibetan.indiv_prov")

tib %>% select(popId, sampleId) %>%
    pgswriter("tib/tibetan.indiv_label")

tib_geo %>% filter(popId %in% pops.to.keep) %>%
    numeric_id %>%
    pgswriter(file="tib/tibetan.pop_geo")
disp <- read.csv("tib/HGDP_Tibetan_Merged_160509.pop_display", as.is=T) %>% 
    filter(popId %in% pops.to.keep) %>%
    numeric_id %>% 
    mutate(name=sapply(name, function(s)strsplit(s, ".", fixed=T)[[1]][1]))

disp$abbrev = c("CN-NB",
                "CN-YB",
                "CN-SB",
                "CN-CX",
                "CN-LX",
                "CN-SH",
                "CN-TS",
                "CN-UM"
                )

pgswriter(disp,  file="tib/tibetan.pop_display")

#disp$name <- sapply(disp$name, function(s)strsplit(s, ".", fixed=T)[[1]][1])
#pgswriter(disp, file="tibetan.pop_display_raw")

