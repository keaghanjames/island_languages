set.seed(888)
library(rworldmap)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(sp)
require(maptools)
require(terra)

setwd("~/Dropbox/Lindell/Island_Languages/LANGUAGE DATA")
df <- as_tibble(read.csv('island_summary.csv'))
#we need to construct a distance matrix and phylo distance matrix
df

### for trimming by occupancy sizes ###
#if we want to delete islands with occupancy sizes less than 0.1km2
islands <- readRDS('island_polygons_united.Rds')

island_areas <- terra::area(islands)
small_islands <- which(island_areas < 0.1)

df <- df[-small_islands,]


require(ggplot2)
require(ggpubr)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggthemes)
require(ggplot2)
world <- map_data("world")
class(world)
require(RColorBrewer)
g <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region), 
    col = '#616365', fill = 'lightgrey', lwd = .1
  ) +   #coord_map("moll") +
  theme_map() +
  geom_point(data = df, aes(x = long, y = lat, col = max_elevation+1)) +
  scale_color_viridis(trans = 'log10', name = 'Max Elv.', option = 'magma')#+
g

g2 <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region), 
    col = '#616365', fill = 'lightgrey', lwd = .1
  ) +   #coord_map("moll") +
  theme_map() +
  geom_point(data = subset(df, island_endemics > 0), aes(x = long, y = lat, col = max_elevation+1)) +
  scale_color_viridis(trans = 'log10', name = 'Max Elv.', option = 'magma')#+
g2
ggsave('~/Dropbox/Lindell/Island_Languages/Figures/max_elevation_island_endemics.pdf')

phyloSor_scores <- readRDS('phyloSor_scores_trim.rds')
distance_matrix <- readRDS('geographic_island_distances_trim.rds')
#if it worked out diaganol should be all 1s
diag(phyloSor_scores)
diag(distance_matrix)

df$treatment <- NULL

#we want to make a scaled version of our dataset for analysis
df_scale <- apply(X = df[,c(4, 7:14, 17)], MARGIN = 2, scale)
#I think LR should be log transformed
df_scale <- as.data.frame(df_scale)
df_scale$WE <- scale(log10(df$WE))
df_scale$LR <- df$lang_rich
df_scale$island_endemics <- df$island_endemics


#i think for now we will remove log area 
df_scale$log_area <- NULL

predictors <- colnames(df_scale)[1:(ncol(df_scale)-3)]
predictors <- predictors[-9]
# fml <- paste0('LR ~ ', paste(predictors, collapse = '+'), '+ weighted_sums + res')
fml <- paste0('island_endemics ~ ', paste(predictors, collapse = '+'), '+ weighted_sums + res')
# fml <- paste0('WE ~ ', paste(predictors, collapse = '+'), '+ weighted_sums + res')

fml <- formula(fml)
fml


# fml_intercept <- paste0('LR ~ 1 + weighted_sums + res')
fml_intercept <- paste0('island_endemics ~ 1 + weighted_sums + res')
# fml_intercept <- paste0('WE ~ 1 + weighted_sums + res')

fml_intercept <- formula(fml_intercept)


df_scale$lat <- df$lat
df_scale$long <- df$long

#we will use fitspphylo and log transform the LR data

#distance_matrix <- readRDS('geographic_island_distances.rds')
#phyloSor_scores <- readRDS('phyloSor_scores.rds')

#we need to get rid of the islands with NA values for MGS from the data but also our two covariance matrices
phylomatrix_cull <- phyloSor_scores[complete.cases(df_scale),complete.cases(df_scale)]
spmatrix_cull <- distance_matrix[complete.cases(df_scale),complete.cases(df_scale)]
df_scale_cull <- df_scale[complete.cases(df_scale),]

dim(df_scale_cull)

#okay so once we take all out complete cases into account we only have 604 data points
plot(df_scale_cull$long, df_scale_cull$lat)

#lets see who is driving this 

df_scale[df_scale == 'NaN'] <- NA

na_counts <- colSums(is.na(df_scale))
na_counts

#excellent, no missing data 


Wsp <- spmatrix_cull
Wphy <- phylomatrix_cull
#we need to standardize the matrices such that they are equivalent in direction (more distant pairs have higher values)
#and scale (all rows sum to one)
Wsp <- spmatrix_cull
diag(Wsp) <- 0
Wrows <- rowSums(Wsp)
for(i in 1:length(Wrows)) {
  Wsp[i,] <- Wsp[i,]/Wrows[i]
}


Wphy <- phylomatrix_cull
Wphy <- phylomatrix_cull/max(phylomatrix_cull)
diag(Wphy)
Wphy <- 1-Wphy #to flip the direction
Wrows <- rowSums(Wphy)
for(i in 1:length(Wrows)) {
  Wphy[i,] <- Wphy[i,]/Wrows[i]
}
rowSums(Wphy)
diag(Wphy)

X <- df_scale_cull[,1:9]
# y <- df_scale_cull$LR
y <- df_scale_cull$island_endemics
# y <- df_scale_cull$WE

X <- as.matrix(X)
X

colnames(X)
colnames(df_scale_cull)


#probably here to loop

# res <- readRDS('res_LR.rds')
res <- readRDS('res_IE.rds')
# res <- readRDS('res_WE.rds')

a <- res$par
a

## now that we know the best fitting model for each outcome variabel lets get the coefficients ###
# load('model_combosLR.Rdata')
load('model_combosIE.Rdata')
# load('model_combosWE.Rdata')

# load('BIC_LR.rdata')
load('BIC_IE.rdata')
# load('BIC_WE.rdata')

tail(model_combos)
fml_full <- model_combos[[length(model_combos)]]

model_combos[[2]]
model_combos <- model_combos[order(unlist(BICs), decreasing = F)]
fml <- model_combos[[1]]


autoglm_diversity <- function (a, y, X, Wsp, Wphy, fml) {
  tmp <- which(!is.na(rowSums(X)))  # identifies rows in X that do not have missing values
  X.tmp <- X[tmp,]  # subsets X to exclude rows with missing values
  y.tmp <- y[tmp]  # subsets y to exclude rows with missing values
  W <- a*Wsp+(1-a)*Wphy
  W <- W[tmp, tmp]  # subsets W to exclude rows and columns with missing values
  Wy <- W %*% as.numeric(y.tmp)  # calculates a weighted sum of y.tmp using W
  X2 <- W %*% X.tmp  # calculates a weighted sum of X.tmp using W
  X2 <- cbind(X.tmp, X2)  # concatenates X.tmp and X2 as columns in a new matrix
  res <- lm(Wy~X2)$residuals  # fits a linear model of Wy on X2 and extracts the residuals
  X1 <- cbind(y.tmp,Wy, X.tmp, res)  # concatenates Wy, X.tmp, and the residuals as columns in a new matrix
  colnames(X1)[1:2] <- c(as.character(fml)[2], 'weighted_sums')
   out2 <- glm(formula = fml, data = as.data.frame(X1), family=poisson(link="log"))
  # out2 <- glm(formula = fml, data = as.data.frame(X1), family=gaussian) #for WE we will use a gaussian distribution
}


# model_summary <- read.csv('model_summariesLR.csv')
model_summary <- read.csv('model_summariesIE.csv')
# model_summary <- read.csv('model_summariesWE.csv')
head(model_summary)

model_summary <- model_summary[which(model_summary$delta <= 6),]
model_summary$X <- NULL
model_summary

model_combos <- model_combos[which(model_summary$delta <= 6)]
model_combos

fitted_models <- lapply(model_combos, function(x) autoglm_diversity(a = res$par, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = x))
fitted_models

# save(fitted_models, file = 'LR_best_models.Rdata')
save(fitted_models, file = 'IE_best_models.Rdata')
# save(fitted_models, file = 'WE_best_models.Rdata')

model_coefficients <- lapply(fitted_models, coef)
model_coefficients

model_confint <- lapply(fitted_models, confint)
model_confint

summaries <- lapply(fitted_models, summary)
summaries[[1]]$coefficients[which(rownames(summaries[[1]]$coefficients) %in%  colnames(model_summary)),4]
summaries

for(i in 1:length(fitted_models)) {
  model_summary[i,which(colnames(model_summary) %in% names(model_coefficients[[i]]))] <- 
    round(model_coefficients[[i]][which(names(model_coefficients[[i]]) %in% colnames(model_summary))], 3)
  
  model_summary[i,which(colnames(model_summary) %in% names(model_coefficients[[i]]))] <- 
    paste0(model_summary[i,which(colnames(model_summary) %in% names(model_coefficients[[i]]))], ' (',
           round(model_confint[[i]][,1][which(rownames(model_confint[[i]]) %in% colnames(model_summary))], 3),
           ' - ',
           round(model_confint[[i]][,2][which(rownames(model_confint[[i]]) %in% colnames(model_summary))], 3),
           ')')
}

model_summary$model_id <- NULL
model_summary
# write.csv(model_summary, 'decent_modelsLR.csv')
write.csv(model_summary, 'decent_modelsIE.csv')
# write.csv(model_summary, 'decent_modelsWE.csv')

#now we will combine all three

LR_models <- read.csv('decent_modelsLR.csv')
IE_models <- read.csv('decent_modelsIE.csv')
WE_models <- read.csv('decent_modelsWE.csv')

Outcome <- c(rep('LR', nrow(LR_models)),
             rep('IE', nrow(IE_models)),
             rep('WE', nrow(WE_models)))

models <- bind_rows(LR_models, IE_models, WE_models)
models$X <- Outcome
colnames(models)[1] <- 'Outcome'


write.csv(models, 'decent_models_diversity.csv', row.names = F)


require(modelsummary)

load('LR_best_models.Rdata')
# load('IE_best_models.Rdata')
# load('WE_best_models.Rdata')

sum <- summary(fitted_models[[1]])
with(sum, 1 - deviance/null.deviance) #get McFadden's psuedo R2

fitted_models[[1]]$formula

#to quantify the individual contribution of phy-spatial auto
fml <- formula('island_endemics ~ dist25k + area_km2 + max_elevation + w_mean_temperature_seasonality ')
# 
glm <- glm(fml, data = fitted_models[[1]]$data, family = 'poisson')
glmSum <- summary(glm)
glmSum
with(glmSum, 1 - deviance/null.deviance) #get McFadden's psuedo R2


is_sig <- function(string){
  string <- strsplit(string, '\\(')[[1]][2]
  string <- unlist(strsplit(string, '-'))
  string[2] <- as.numeric(gsub(x = string[2], pattern = ')', replacement = ''))
  string[1] <- as.numeric(string[1])
  string <- as.numeric(string)
  sig <- between(0, string[1], string[2])
  return(sig == F) 
}

siggy <- lapply(models$area_km2, is_sig)
siggy


