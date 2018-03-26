require(dplyr)
require(yaml)
P <- function(y, x)print(sprintf("%s : %s", x,y))

x <- read.csv("pgs/gvar3.indiv_label")    
y <- read.csv("pgs/gvar3.pop_geo")         
pd <- read.csv("pgs/gvar3.pop_display")
data <- left_join(x,y)                       

C <- yaml.load_file("~/n/eems_tib/config/data.yaml")$filter

PP <- function(s){data %>% filter(popId %in% C[[s]]) %>% nrow %>% P(s)}

sapply(names(C), PP)
P("Americans", sum(data$longitude <(-35)))

a <- data %>% filter(longitude > -35, !popId %in% C$clean)

indiv_meta <- readRDS("~/n/eems_tib/paper/panels.rds")
in_paper_panels <- do.call(bind_rows, indiv_meta) 

#get samples that are not yet accounted for
 a %>% filter(!sampleId %in% in_paper_panels$sampleId) %>% left_join(pd)-> data2   

#counts of unaccounted
cnts <- data2 %>% group_by(name, popId) %>% summarize(n=n()) %>% arrange(n) %>% as.data.frame

