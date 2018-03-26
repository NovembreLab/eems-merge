library(dplyr)

#merge guys
il <- read.csv("pgs/gvar3.indiv_label")
to_merge <- read.csv("pgs/merge.csv")
to_remove <- to_merge %>% select(popId)

il %>% left_join(to_merge) %>%
    mutate(popId=ifelse(is.na(target), popId, target)) %>% 
    select(-target) %>%
    write.csv("pgs/gvar3.indiv_label", row.names=F, quote=T)

print("done 1")

#update pos
pg <- read.csv("pgs/gvar3.pop_geo")
fix_pos <- read.csv("pgs/update_pos.csv")
pg %>% left_join(fix_pos) %>%
    mutate(latitude=ifelse(is.na(new.lat), latitude, new.lat)) %>% 
    mutate(longitude=ifelse(is.na(new.long), longitude, new.long)) %>%
    mutate(latitude=round(latitude,2))%>%
    mutate(longitude=round(longitude,2))%>%
    select(-new.lat, -new.long) %>% 
    filter(! popId %in% to_remove$popId) %>%
    write.csv("pgs/gvar3.pop_geo", row.names=F, quote=T)
print("done 2")

#update names
pg <- read.csv("pgs/gvar3.pop_display", as.is=T)
better_names <- read.csv("pgs/gvar3.names", as.is=T)
pg %>% left_join(better_names) %>%
    mutate(name=ifelse(is.na(newName), name, newName)) %>% 
    select(-newName) %>% 
    filter(! popId %in% to_remove$popId) %>%
    write.csv("pgs/gvar3.pop_display", row.names=F, quote=T)
