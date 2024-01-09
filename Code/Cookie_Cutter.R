set.seed(888)
library(rworldmap)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(sp)
require(terra)

setwd("~/Dropbox/Lindell/Island_Languages/LANGUAGE DATA")
lang_map <- readRDS("global_intersection_resolvedNEW.rds")
head(lang_map@data)
require(dplyr)

#mainlands <- readRDS('mainlands_25km2.rds')
mainlands <- readRDS('continents.rds')
mainlands <- gUnion(mainlands, mainlands)
mainlands
#plot(mainlands)
crs(mainlands)

get_radius <- function(x){
  (sqrt(x/pi))
}

get_radius(120)



eaproj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"

#so we will need centroids for our islands
require(maptools)
islands <- subset(lang_map, island11kkm2 == 'Island')#[1:10,]


get_cookies_w_mainland_coverage <- function(mainlands, island, tol = 0.1, it_limit = 200, y_pos = NULL, y_buffer = 0){
  control = 0
  while (control < it_limit) {
    if(is.null(y_pos) == F){
      bb <- bbox(mainlands)
      bb[2,1:2] <- y_pos
      bb[2,1] <- bb[2,1]-y_buffer #here we add our latitudinal buffer
      bb[2,1] <- bb[2,2]+y_buffer
      #and we same a point at random from our crop
      point <- try(spsample(crop(mainlands, extent(bb)), n = 1, 'random'))
    } else {
      point <- try(spsample(mainlands, n = 1, 'random'))
    }
    if(class(point) == 'try-error') {
      cookie = point
      return(cookie)
      next
    }
    point
    translation_vector <- as.vector(coordinates(point) - coordinates(gCentroid(island)))
    cookie <- shift(island, dx = translation_vector[1], dy = translation_vector[2])
    
    
    x <- st_intersection(st_as_sf(cookie) , st_as_sf(mainlands)) %>%
      mutate(intersect_area = st_area(.)) %>%
      st_drop_geometry()  # drop geometry as we don't need it
    if((as.numeric(x$intersect_area)/area(cookie) >= (1-tol)) == T) {
      control = it_limit
    } else {
      control = control+1
      cookie = 'timed-out'
      print("too much ocean, cutting new cookie")
    }
  }
  return(cookie)
  
}



get_cookies_handler <- function(island){
  y_pos <- coordinates(gCentroid(island))[[2]]
  out <- get_cookies_w_mainland_coverage(mainlands = mainlands, island = island, y_pos = y_pos, y_buffer = 500)
  return(out)
}



#### You are up to here ####



require(rgdal)
require(pbmcapply)
#we want a function that just takes the run ID and radii vector as input
run_experiment <- function(islands, mainlands = mainlands, map = lang_map, it_limit = 200, y_buffer = 500, tol = 0.1, nick = 0.04, cave = NULL) {
  get_cookies_handler <- function(island){
    y_pos <- coordinates(gCentroid(island))[[2]]
    isle <- island
    out <- get_cookies_w_mainland_coverage(mainlands = mainlands, island = isle, y_pos = y_pos, y_buffer = y_buffer, it_limit = it_limit)
    out$LM_ID <- island$LM_ID
    return(out)
  }
  synthetic_islands <- lapply(split(islands, 1:length(islands)), FUN = try(get_cookies_handler))
  names(synthetic_islands) <- NULL
  
  #drop all the ones which did not converge
  islands <- islands[which(unlist(lapply(synthetic_islands, class)) == 'SpatialPolygonsDataFrame'),]
  synthetic_islands <- synthetic_islands[which(unlist(lapply(synthetic_islands, class)) == 'SpatialPolygonsDataFrame')]
  synthetic_islands <- do.call(bind, synthetic_islands)
  areas <- area(synthetic_islands)
  
  co <- coordinates(synthetic_islands)
  co <- SpatialPoints(co, CRS(eaproj))
  co <- spTransform(co, CRS("+proj=longlat"))
  
   
  synthetic_islands$Cookie_ID = synthetic_islands$LM_ID
  synthetic_islands$cookie_area = areas
  synthetic_islands$island_area = area(islands)
  synthetic_islands$long = co@coords[,1]
  synthetic_islands$lat = co@coords[,2]
  
  inter <- raster::intersect(synthetic_islands, lang_map)
  inter <- inter[inter@data$island11kkm2 == 'Mainland',]  #on the off chance our cookie samples a mainland and an island we drop the island component
  inter@data$intersect_area <- area(inter)
  
  
  if(is.null(nick) == F) {
  #we probably want to run some control for 'nicks' languages which just intersect with a island
  #inter$inter_cookie_area_km2 <- area(inter) #we get the are of every intersection
  inter$prop_cookie <- inter@data$intersect_area/inter@data$cookie_area #and calculate their proportion of the cookies
  #and now we can use this to drop any languages below a certain threshold
  
  if(is.null(cave) == F) {
  inter$prop_cave <- inter@data$intersect_area/inter@data$LANG_area_km2
  inter <- inter[-which(inter$prop_cookie < nick & inter$prop_cave < cave),]
  
  } else {
    inter <- inter[-which(inter$prop_cookie < nick),]
    }
  
  }
  
  
  return(inter)
}



#okay for each run we want to do the following
#1 - Sample N islands
#2 - Cut N Cookies and intersect with lang_map with the run_experiment function
#3 - Intersect the sample of islands with langmap
#4 - Save the cookies' intercept
#5 - Save the islands' intercept
#remake our full islands spatial object
require(bayesbio)
lang_map$special_IDs <- factor(lang_map$LM_ID)
levels(lang_map$special_IDs) <- createStrings(number = length(levels(lang_map$special_IDs)), length = 4)
levels(lang_map$special_IDs)
lang_map$LM_ID <- as.integer(lang_map$LM_ID)

islands <- subset(lang_map, island11kkm2 == 'Island')#[1:10,]
LM_ID <- unique(islands@data$LM_ID)
special_IDs <- unique(islands@data$special_IDs)
islands <- unionSpatialPolygons(islands, IDs = islands@data$special_IDs)
islands <- islands[match(special_IDs, names(islands)),]
data <- data.frame(LM_ID, special_IDs)
rownames(data) <- data$special_IDs
islands <- SpatialPolygonsDataFrame(islands, data = data)


run_experiment_handler <- function(run, N,islands = islands, mainlands = mainlands, tol = 0.1, y_buffer = 500, map = lang_map, nick = 0.04, cave = NULL, it_limit = 200){
#1 - Sample N islands
  islands <- islands[(sample(1:length(islands), N)),]
  islands$area <- area(islands)
  islands$expID <- 1:nrow(islands)
#2 - Cut N Cookies and intersect with lang_map with run experiment function
  exp <- run_experiment(islands = islands, mainlands = mainlands, tol = tol, y_buffer = y_buffer, map = map, nick = nick, it_limit = it_limit)
#  exp$Run <- run
#3 - Intersect the islands sample with lang map  
  islands <- raster::intersect(islands, lang_map)
#4 - we remove languages which account for less than our nick size of an island to be consistent with the cookies
  islands$intersect_area <- area(islands)
  islands$intersection_prop <- islands@data$intersect_area/islands@data$area #and calculate their proportion of the cookies    
  islands <- islands[which(islands$intersection_prop >= nick),]
#5 - Save the cookies intersect
  saveRDS(exp, file = paste0('cookies_run',run,'.rds'))
#6 - Save the islands intersect
  saveRDS(islands, file = paste0('islands_run',run,'.rds'))
#7 - return a list with both outputs
  output <- list()
  output$islands <- islands
  output$cookies <- exp
  return(output)
}

set.seed(888)

test <- run_experiment_handler(run = 1, N = 13160, islands = islands, mainlands = mainlands, 
                               tol = 0.1, y_buffer = 500, map = lang_map, nick = .04, it_limit =100)

#so if we ask which IDs in the islands are not present in the cookies we will know which islands couldn't be turned into cookies
fails <- unique(test$islands$LM_ID)[which(unique(test$islands$LM_ID) %in% unique(test$cookies$LM_ID) == F)]
test$islands$failed <- F
test$islands$failed[which(test$islands$LM_ID %in% fails)] <- T
fails <- test$islands[which(test$islands$LM_ID %in% fails),]
fails_cent <- gCentroid(fails, byid = T)
fails_cent <- coordinates(fails_cent)
 

#lets drop the failure
test$islands <- test$islands[-which(test$islands$LM_ID %in% fails$LM_ID),]

#great we have got it down to a single island that failed

test
# save(test, file = 'Cookie_Cutter_Completish_25kmNEW.Rdata')
# load('Cookie_Cutter_Completish_25kmNEW.Rdata')
test
save(test, file = 'Cookie_Cutter_Completish_Continent.Rdata')
load('Cookie_Cutter_Completish_Continent.Rdata')


which(unique(test$cookies@data$Cookie_ID) %in% unique(test$islands@data$LM_ID) == F)

missing_islands <- test$islands[which(test$islands$LM_ID %in% test$cookies$Cookie_ID == F),]
# plot(mainlands)
# plot(missing_islands, add = T, bor = 'red')

#okay so even after the while loop we have 219 islands for which we have no cookie
#honestly that is so few I think we can just go ahead and remove them from the dataset
test$islands <- test$islands[which(test$islands$LM_ID %in% test$cookies$Cookie_ID == T),]


# # Mainlands <- lang_map[which(lang_map$LM_area_km2 > 25000),]
test


#next we will summarise the results of our analysis 

  cookie <- test$cookies
  cookie$WE <- 1/cookie$LANG_area_km2
  cookie$PoR <- area(cookie)/cookie$LANG_area_km2
  
  islands <- test$islands
  islands$WE <- 1/islands$LANG_area_km2
  islands$PoR <- area(islands)/islands$LANG_area_km2
  
  
  require(tidyr)
  require(dplyr)
  sum_cookie <- cookie@data %>%
    group_by(LM_ID) %>%
    summarise(lang_rich = length(unique(LANG_ISO)), 
              WE = sum(WE), area = mean(cookie_area),
              sPoR = sum(PoR),
              long = mean(long), lat = mean(lat))
  sum_cookie$treatment <- 'Cookie'
  sum_cookie
  
  coordinates <- gCentroid(islands, byid = T)
  coordinates <- spTransform(x = coordinates, CRS("+proj=longlat"))
  coordinates <- coordinates(coordinates)
  islands$long <- coordinates[,1]
  islands$lat <- coordinates[,2]
  
  sum_island <- islands@data %>%
    group_by(LM_ID) %>%
    summarise(lang_rich = length(unique(LANG_ISO)), 
              WE = sum(WE), area = mean(area),
              sPoR = sum(PoR),
              long = mean(long), lat = mean(lat))
  sum_island$treatment <- 'Island'

  summary <- bind_rows(sum_island, sum_cookie)
  summary$joint_IDs <- summary$LM_ID

  
require(ggplot2)
islands_cookies <-   ggplot(summary, aes(long,lat)) +
                      geom_point(aes(color=treatment), size=1, alpha = .3) +
                      geom_line(aes(group = joint_IDs), linewidth = .01)  +
                      theme_classic() +
                      xlab('Longitude') +
                      ylab('Latitude') +
                      theme(legend.title = element_blank(),
                            text = element_text(size=20)) +
  coord_equal()
islands_cookies            
ggsave('~/Dropbox/Lindell/Island_Languages/LANGUAGE DATA/cookie_cutter_ContinentsNEW.pdf', device = 'pdf')
length(unique(sum_island$LM_ID) )
  

#write.csv(summary, 'cookie_cutter_summary_25km.csv')
write.csv(summary, 'cookie_cutter_summary_Continent.csv')


#build tables summarising the results
summary25km <- read.csv('cookie_cutter_summary_25km.csv')
summaryContinent <- read.csv('cookie_cutter_summary_Continent.csv')

summary25km$model <- '25km2'
summaryContinent$model <- 'Continent'

mainlands_lat_long <- st_transform(st_as_sf(mainlands),  CRS('+proj=longlat'))  
mainlands_lat_long <- fortify(mainlands_lat_long)


summary11k$X <- NULL
summaryContinent$X <- NULL

summary25km$X <- NULL
cookie_cutter_summary <- rbind(summary25km, summaryContinent)
cookie_cutter_summary
write.csv(cookie_cutter_summary, file = 'cookie_cutter_summary.csv')
colnames(summary25km)
colnames(summaryContinent)





### ANALYSIS STARTS HERE #######
set.seed(888)
library(rworldmap)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(sp)
require(terra)

setwd("~/Dropbox/Lindell/Island_Languages/LANGUAGE DATA")


library("rnaturalearth")
library("rnaturalearthdata")
library(ggthemes)
require(ggplot2)
world <- map_data("world")
class(world)


df <- read.csv('cookie_cutter_summary.csv')
summary25k <- read.csv('cookie_cutter_summary_25km.csv')
summaryContinent <- read.csv('cookie_cutter_summary_Continent.csv')

summary25k$model <- '25km2'
summaryContinent$model <- 'Continent'
head(df)


#outlier detection, so i noticed one of the cookies
hist(df$lang_rich)
plot(df$area, df$lang_rich)

df <- df[order(df$lang_rich, decreasing = T),]
head(df)

#its absolutely an outlier so we may drop it, but check with Lindell 

#df <- df[-1,]


#the coordinates of that cookie are   -6.414658, 145.4960 and sure enough 
#that falls right within Papua New Guinea

mainlands <- readRDS('global_mainlands.rds')
#plot(mainlands)


# we need to trim to 10,000 comparisons, but we need to make sure its the same 10,000 between the islands and cookies
sample <- sample(1:(nrow(summary25k)/2), 10000)
sample <- c(sample, sample+(nrow(summary25k)/2))

summary25k <- summary25k[sample,]

sample <- sample(1:(nrow(summaryContinent)/2), 10000)
sample <- c(sample, sample+(nrow(summaryContinent)/2))

summaryContinent <- summaryContinent[sample,]


df <- rbind(summary25k, summaryContinent)

df
# samp <- sample(1:nrow(summary11k), 5000)
# samp



islands_cookiesWE <-   ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region), 
    col = '#616365', fill = 'lightgrey', lwd = .1
  ) +   #coord_map("moll") +
  theme_map() +
  geom_point(data = summary25k, aes(x = long, y = lat, color=scale(log10(WE)), shape = treatment), size=.5, alpha = .7) +
#  geom_line(aes(group = joint_IDs), linewidth = .01)  +
  theme_classic() +
  xlab('longitude') +
  ylab('latitude') +
  scale_color_gradient(low = '#2D95F6', high = 'red', name = 'log10(WE)') +
  xlim(-180, 180) +
  ylim(-90, 90) +
  scale_shape_discrete(name = '') +
  coord_equal()
  #theme(legend.title = element_blank())
islands_cookiesWE


islands_cookiesLR <-  ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region), 
    col = '#616365', fill = 'lightgrey', lwd = .1
  ) +   #coord_map("moll") +
  theme_map() +
  geom_point(data = summary25k, aes(x = long, y = lat, color=scale((lang_rich)), shape = treatment), size=.5, alpha = .7) +
  #  geom_line(aes(group = joint_IDs), linewidth = .01)  +
  theme_classic() +
  xlab('longitude') +
  ylab('latitude') +
  scale_color_gradient(low = '#2D95F6', high = 'red', name = LR) +
  xlim(-180, 180) +
  ylim(-90, 90) +
  scale_shape_discrete(name = '')+
  coord_equal()
islands_cookiesLR

islands_cookiesPoR <-  ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region), 
    col = '#616365', fill = 'lightgrey', lwd = .1
  ) +   #coord_map("moll") +
  theme_map() +
  geom_point(data = summary25k, aes(x = long, y = lat, color=scale(log10(sPoR)), shape = treatment), size=.5, alpha = .7) +
  #  geom_line(aes(group = joint_IDs), linewidth = .01)  +
  theme_classic() +
  xlab('longitude') +
  ylab('latitude') +
  scale_color_gradient2(low = '#2D95F6', mid = '#2D95F6', high = 'red', midpoint = -1, name = 'log10(PoR)') +
  xlim(-180, 180) +
  ylim(-90, 90) +
  scale_shape_discrete(name = '')+
  coord_equal()
islands_cookiesPoR





WE_stripplot <- ggplot(data = summary25k, aes(x = treatment, y = WE, shape = treatment, col = scale(log10(WE))))+
  geom_jitter(position=position_jitter(0.2), alpha = .7)+
  geom_violin(width = .3, col = '#616365', alpha = .7) +
  scale_y_continuous(trans = 'log10', name = 'log10(WE)') +
  scale_color_gradient(low = '#2D95F6', high = 'red') +
  theme_classic() +
  xlab('') +
  theme(legend.position = 'none')
WE_stripplot


LR_stripplot <- ggplot(data = summary25k, aes(x = treatment, y = lang_rich, shape = treatment, col = scale(log10(lang_rich))))+
  geom_jitter(position=position_jitter(0.2), alpha = .7)+
  geom_violin(width = 1, outlier.shape = NA, col = '#616365', alpha = .7) +
  scale_y_continuous(breaks = c(1:10)) +
  scale_color_gradient(low = '#2D95F6', high = 'red') +
  theme_classic() +
  xlab('') +
  ylab('LR') +
  theme(legend.position = 'none')
LR_stripplot



PoR_stripplot <- ggplot(data = summary25k, aes(x = treatment, y = sPoR, shape = treatment, col = scale(log10(sPoR))))+
  geom_jitter(position=position_jitter(0.2), alpha = .7)+
  geom_violin(width = 1, outlier.shape = NA, col = '#616365', alpha = .7) +
  #scale_y_continuous(trans = 'log10', name = 'log10(PoR)') +
  scale_color_gradient2(low = '#2D95F6', mid = '#2D95F6', high = 'red', midpoint = -1, name = 'log10(PoR)') +
  theme_classic() +
  xlab('') +
  ylab('Sum PoR') +
  theme(legend.position = 'none')
PoR_stripplot


ggplot(summary25k, aes(x = log10(sPoR), y = log10(WE))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  facet_wrap(~ treatment)

cor.test(log10(summary25k$sPoR), log10(summary25k$WE))

lm <- lm(log10(summary25k$sPoR) ~ log10(summary25k$WE))
summary(lm)

#so we probably want to do some statistical tests to show that these are indeed different

require(dplyr)
summary <- df %>% 
            group_by(model, treatment) %>%
            summarise(meanLR = mean(lang_rich), sdLR = sd(lang_rich),
                      meanWE = mean(WE), sdWE = sd(WE),
                      meanPoR = mean(sPoR), sdPoR = sd(sPoR))

summary

#lets get ecdf for our indices
#we need to create a new variable that encodes model and treatment

df$model_treatment <- paste0(df$treatment, '_', df$model)

LR_ecdf <- lapply(split(df$lang_rich, df$model_treatment), ecdf)


LR_cdf <- ggplot(data = subset(df, model == '25km2'), aes(x = lang_rich, col = treatment)) + 
  stat_ecdf() +
  scale_color_manual(values = c('black', 'darkgrey'), name = "") +
  xlab('Language Richness')+
  ylab('F[Language Richness]')+
 scale_x_continuous(trans = 'log10') +
  #xlim(0, 5e+06) +
  theme_classic() 
  #theme(legend.position = 'NULL')+
#require(ggpubr)
#ggarrange(p1, p1_ln, nrow = 1, common.legend = T)
LR_cdf

WE_cdf <- ggplot(data = subset(df, model == '25km2'), aes(x = WE,  col = treatment)) + 
  stat_ecdf() +
  scale_color_manual(values = c('black', 'darkgrey'), name = "") +
  xlab('Weighted Endemism')+
  ylab('F[Weighted Endemism]')+
  scale_x_continuous(trans = 'log10') +
  #xlim(0, 5e+06) +
  theme_classic() 
#theme(legend.position = 'NULL')+
#require(ggpubr)
#ggarrange(p1, p1_ln, nrow = 1, common.legend = T)
WE_cdf

PoR_cdf <- ggplot(data = subset(df, model == '25km2'), aes(x = sPoR,  col = treatment)) + 
  stat_ecdf() +
  scale_color_manual(values = c('black', 'darkgrey'), name = "") +
  xlab('PoR')+
  ylab('F[PoR]')+
  scale_x_continuous(trans = 'log10') +
  #xlim(0, 5e+06) +
  theme_classic() 
#theme(legend.position = 'NULL')+
#require(ggpubr)
#ggarrange(p1, p1_ln, nrow = 1, common.legend = T)
PoR_cdf



require(ggpubr)

figure1 <- ggarrange(islands_cookiesLR, LR_stripplot,LR_cdf,
                     islands_cookiesWE,  WE_stripplot, WE_cdf,
          ncol = 3, nrow = 2, widths = c(1, .6, .6),
          labels = c('a', 'c', 'e', 'b', 'd', 'f'))
figure1
ggsave('figure1_25km.pdf', plot = figure1, units = 'mm')

figure1.1 <- ggarrange(islands_cookiesLR, LR_stripplot,LR_cdf,
                     islands_cookiesWE,  WE_stripplot, WE_cdf,
                     islands_cookiesPoR, PoR_stripplot,PoR_cdf,
                      
                     ncol = 3, nrow = 3, widths = c(1, .6, .6),
                     labels = c('a', 'd', 'g', 'b', 'e', 'h', 'c', 'f', 'i'))
figure1.1
ggsave('figure1_25km.pdf', plot = figure1.1, units = 'mm')


stripplots <- ggarrange(LR_stripplot, WE_stripplot, PoR_stripplot, nrow = 1, ncol = 3, 
                        labels = c('d', 'e', 'f'))
stripplots

figure1.2 <- ggarrange(islands_cookiesLR, 
                       islands_cookiesWE,   
                       islands_cookiesPoR, 
                       stripplots,
                       ncol = 1, nrow = 4, heights = c(1,1,1,.8),
                       labels = c('a', 'b', 'c', '')) 
figure1.2
ggsave('~/Dropbox/Lindell/Island_Languages/Figures/figure1_25kmNEW.pdf', plot = figure1.2, units = 'mm', width = 210, height = 270, 
      scale = 1)
save(figure1.2, file = '~/Dropbox/Lindell/Island_Languages/Figures/figure1.2.Rdata')

#next we will make the same fiugure but for the continental treatment, this will probably just go in the supps
islands_cookiesWE <-   ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region), 
    col = '#616365', fill = '#616365', lwd = .1
  ) +   #coord_map("moll") +
  theme_map() +
  geom_point(data = subset(df, model == 'Continent'), aes(x = long, y = lat, color=WE, shape = treatment), size=.5, alpha = .7) +
  #  geom_line(aes(group = joint_IDs), linewidth = .01)  +
  theme_classic() +
  xlab('longitude') +
  ylab('latitude') +
  scale_color_gradient(low = '#2D95F6', high = '#FFAFCC', trans = 'log10', name = 'log10(WE)') +
  xlim(-180, 180) +
  ylim(-90, 90) +
  scale_shape_discrete(name = '')
#theme(legend.title = element_blank())
islands_cookiesWE


islands_cookiesLR <-  ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region), 
    col = '#616365', fill = '#616365', lwd = .1
  ) +   #coord_map("moll") +
  theme_map() +
  geom_point(data = subset(df, model == 'Continent'), aes(x = long, y = lat, color=lang_rich, shape = treatment), size=.5, alpha = .7) +
  #  geom_line(aes(group = joint_IDs), linewidth = .01)  +
  theme_classic() +
  xlab('longitude') +
  ylab('latitude') +
  scale_color_gradient(low = '#2D95F6', high = '#FFAFCC', trans = 'log10', name = 'log10(LR)') +
  xlim(-180, 180) +
  ylim(-90, 90) +
  scale_shape_discrete(name = '')
islands_cookiesLR


WE_stripplot <- ggplot(data = subset(df, model == 'Continent'), aes(x = treatment, y = WE, shape = treatment, col = WE))+
  geom_jitter(position=position_jitter(0.2), alpha = .7)+
  geom_violin(width = .3, col = '#616365', alpha = .7) +
  scale_y_continuous(trans = 'log10', name = 'log10(Weighted Endemism)') +
  scale_color_gradient(low = '#2D95F6', high = '#FFAFCC', trans = 'log10') +
  theme_classic() +
  xlab('') +
  theme(legend.position = 'none')
WE_stripplot


LR_stripplot <- ggplot(data = subset(df, model == 'Continent'), aes(x = treatment, y = lang_rich, shape = treatment, col = lang_rich))+
  geom_jitter(position=position_jitter(0.2), alpha = .7)+
  geom_violin(width = 1, outlier.shape = NA, col = '#616365', alpha = .7) +
  scale_y_continuous(breaks = c(1:10)) +
  scale_color_gradient(low = '#2D95F6', high = '#FFAFCC',  trans = 'log10') +
  theme_classic() +
  xlab('') +
  ylab('Language Richness') +
  theme(legend.position = 'none')
LR_stripplot

LR_cdf <- ggplot(data = subset(df, model == 'Continent'), aes(x = lang_rich, col = treatment)) + 
  stat_ecdf() +
  scale_color_manual(values = c('black', 'darkgrey'), name = "") +
  xlab('Language Richness')+
  ylab('F[Language Richness]')+
  scale_x_continuous(trans = 'log10') +
  #xlim(0, 5e+06) +
  theme_classic() 
#theme(legend.position = 'NULL')+
#require(ggpubr)
#ggarrange(p1, p1_ln, nrow = 1, common.legend = T)
LR_cdf

WE_cdf <- ggplot(data = subset(df, model == 'Continent'), aes(x = WE,  col = treatment)) + 
  stat_ecdf() +
  scale_color_manual(values = c('black', 'darkgrey'), name = "") +
  xlab('Weighted Endemism')+
  ylab('F[Weighted Endemism]')+
  scale_x_continuous(trans = 'log10') +
  #xlim(0, 5e+06) +
  theme_classic() 
#theme(legend.position = 'NULL')+
#require(ggpubr)
#ggarrange(p1, p1_ln, nrow = 1, common.legend = T)
WE_cdf


require(ggpubr)

figure1 <- ggarrange(islands_cookiesLR, LR_stripplot,LR_cdf,
                     islands_cookiesWE,  WE_stripplot, WE_cdf,
                     ncol = 3, nrow = 2, widths = c(1, .6, .6),
                     labels = c('a', 'c', 'e', 'b', 'd', 'f'))
figure1
ggsave('figure1_Continent.pdf', plot = figure1,  units = 'mm', scale = 1.5)


#next we want to actually run some statistical test


df_split <- split(df, df$model_treatment)

write.csv(summary, 'cookie_cutter_summary.csv')

ks.test(df_split$Cookie_25km2$lang_rich, df_split$Island_25km2$lang_rich)
ks.test(df_split$Cookie_25km2$WE, df_split$Island_25km2$WE)
ks.test(df_split$Cookie_25km2$sPoR, df_split$Island_25km2$sPoR)

mean(df_split$Cookie_25km2$lang_rich-df_split$Island_25km2$lang_rich)
sd(df_split$Cookie_25km2$lang_rich-df_split$Island_25km2$lang_rich)

mean(df_split$Cookie_25km2$WE-df_split$Island_25km2$WE)
sd(df_split$Cookie_25km2$WE-df_split$Island_25km2$WE)

mean(df_split$Cookie_25km2$sPoR-df_split$Island_25km2$sPoR)
sd(df_split$Cookie_25km2$sPoR-df_split$Island_25km2$sPoR)

ks.test(df_split$Cookie_Continent$lang_rich, df_split$Island_Continent$lang_rich)
ks.test(df_split$Cookie_Continent$WE, df_split$Island_Continent$WE)
ks.test(df_split$Cookie_Continent$sPoR, df_split$Island_Continent$sPoR)

mean(df_split$Cookie_Continent$lang_rich-df_split$Island_Continent$lang_rich)
sd(df_split$Cookie_Continent$lang_rich-df_split$Island_Continent$lang_rich)

mean(df_split$Cookie_Continent$WE-df_split$Island_Continent$WE)
sd(df_split$Cookie_Continent$WE-df_split$Island_Continent$WE)

mean(df_split$Cookie_Continent$sPoR-df_split$Island_Continent$sPoR)
sd(df_split$Cookie_Continent$sPoR-df_split$Island_Continent$sPoR)



any(df_split$Cookie_25km2$Cookie_ID %in% df_split$Island_25km2$LM_ID == F)
any(df_split$Island_25km2$LM_ID %in% df_split$Cookie_25km2$Cookie_ID == F)
which(df_split$Island_25km2$LM_ID %in% df_split$Cookie_25km2$Cookie_ID == F)


#so for now lest remove it so we can do pairwise comparisons
#df_split$Island_25km2 <- df_split$Island_25km2[-613,]
#and just check we did it correctly 
#which(df_split$Island_25km2$LM_ID %in% df_split$Cookie_25km2$Cookie_ID == F)

#we didn't cut anything from the Continent treatment so we should be Gucci 

df_split$Cookie_25km2 <- df_split$Cookie_25km2[match(df_split$Island_25km2$LM_ID, df_split$Cookie_25km2$Cookie_ID),]
head(df_split$Island_25km2$LM_ID)
head(df_split$Cookie_25km2$Cookie_ID)
df_split$Cookie_Continent <- df_split$Cookie_Continent[match(df_split$Island_Continent$LM_ID, df_split$Cookie_Continent$Cookie_ID),]


#paired LR tests
wilcox.test(df_split$Island_25km2$lang_rich, df_split$Cookie_25km2$lang_rich, paired = T)
wilcox.test(df_split$Island_Continent$lang_rich, df_split$Cookie_Continent$lang_rich, paired = T)

wilcox.test(df_split$Island_25km2$WE, df_split$Cookie_25km2$WE, paired = T)
wilcox.test(df_split$Island_Continent$WE, df_split$Cookie_Continent$WE, paired = T)

wilcox.test(df_split$Island_25km2$sPoR, df_split$Cookie_25km2$sPoR, paired = T)
wilcox.test(df_split$Island_Continent$sPoR, df_split$Cookie_Continent$sPoR, paired = T)


#yep they are also all significant

dif_25km2LR <- df_split$Island_25km2$lang_rich - df_split$Cookie_25km2$lang_rich
dif_ContLR <- df_split$Island_Continent$lang_rich - df_split$Cookie_Continent$lang_rich
dif_25km2WE<- df_split$Island_25km2$WE - df_split$Cookie_25km2$WE
dif_ContWE <- df_split$Island_Continent$WE - df_split$Cookie_Continent$WE
prop_25km2WE<- (df_split$Cookie_25km2$WE/df_split$Island_25km2$WE)
prop_ContWE <- (df_split$Cookie_Continent$WE/df_split$Island_Continent$WE)



c(mean(dif_25km2LR), sd(dif_25km2LR))
c(mean(dif_ContLR), sd(dif_ContLR))

c(mean(dif_25km2WE), sd(dif_25km2WE))
c(mean(dif_ContWE), sd(dif_ContWE))

c(mean(prop_25km2WE), sd(prop_25km2WE))
c(mean(prop_ContWE), sd(prop_ContWE))
