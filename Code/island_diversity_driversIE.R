
       ##############################################
########### Analysis #### Analysis #### Analysis ###########
       ##############################################
set.seed(888)
library(rworldmap)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(sp)
require(maptools)
require(terra)
require(dplyr)

setwd("~/Dropbox/Lindell/Island_Languages/LANGUAGE DATA")
df <- as_tibble(read.csv('island_summary.csv'))
#we need to construct a distance matrix and phylo distance matrix
df

### for trimming by occupancy sizes ###
#if we want to delete islands with occupancy sizes less than 0.1km2
islands <- readRDS('island_polygons_united.Rds')

island_areas <- area(islands)
small_islands <- which(island_areas < 0.1)

df <- df[-small_islands,]

require(rgeos)
require(fields)



phyloSor_scores <- readRDS('phyloSor_scores_trim.rds')
distance_matrix <- readRDS('geographic_island_distances_trim.rds')
#if it worked out diaganol should be all 1s
diag(phyloSor_scores)
diag(distance_matrix)

dim(df)
dim(distance_matrix)
dim(phyloSor_scores)
any(colnames(phyloSor_scores) != df$special_IDs)
any(colnames(distance_matrix) != df$special_IDs)

#brilliant, the matrices are in exactly the same order as the df


colnames(df)
df$treatment <- NULL
df$geo_islolation <- NULL
df$phy_islolation <- NULL
colnames(df)

#we want to make a scaled version of our dataset for analysis
df_scale <- apply(X = df[,c(4, 7:14, 17)], MARGIN = 2, scale)
#I think LR should be log transformed
df_scale <- as.data.frame(df_scale)
df_scale$WE <- df$WE
df_scale$LR <- df$lang_rich
df_scale$island_endemics <- df$island_endemics


#i think for now we will remove log area 
df_scale$log_area <- NULL

predictors <- colnames(df_scale)[1:(ncol(df_scale)-3)]
predictors
fml <- paste0('island_endemics ~ ', paste(predictors, collapse = '+'), '+ weighted_sums + res')
fml <- formula(fml)
fml


fml_intercept <- paste0('island_endemics ~ 1 + weighted_sums + res')
fml_intercept <- formula(fml_intercept)


df_scale$lat <- df$lat
df_scale$long <- df$long


#we need to get rid of the islands with NA values for MGS from the data but also our two covariance matrices
phylomatrix_cull <- phyloSor_scores[complete.cases(df_scale),complete.cases(df_scale)]
spmatrix_cull <- distance_matrix[complete.cases(df_scale),complete.cases(df_scale)]
df_scale_cull <- df_scale[complete.cases(df_scale),]

dim(df_scale_cull)


df_scale[df_scale == 'NaN'] <- NA

na_counts <- colSums(is.na(df_scale))
na_counts

#excellent, no missing data 


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
#y <- round(log(df_scale_cull$WE), digits = 5)
y <- df_scale_cull$island_endemics

X <- as.matrix(X)
X

colnames(X)
colnames(df_scale_cull)
a <- 0.5

fml <- formula('island_endemics ~ weighted_sums + res')

autoglm_diversity <- function (a, y, X, Wsp, Wphy, fml) {
  tmp <- which(!is.na(rowSums(X)))  # identifies rows in X that do not have missing values
  X.tmp <- X[tmp,]  # subsets X to exclude rows with missing values
  y.tmp <- y[tmp]  # subsets y to exclude rows with missing values
  #W <- a*Wsp+(1-a)*Wphy  # calculates weights using values in a, Wsp, Wnb, and Wphy
  #W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
  W <- a*Wsp+(1-a)*Wphy
  W <- W[tmp, tmp]  # subsets W to exclude rows and columns with missing values
  Wy <- W %*% as.numeric(y.tmp)  # calculates a weighted sum of y.tmp using W
  X2 <- W %*% X.tmp  # calculates a weighted sum of X.tmp using W
  X2 <- cbind(X.tmp, X2)  # concatenates X.tmp and X2 as columns in a new matrix
  res <- lm(Wy~X2)$residuals  # fits a linear model of Wy on X2 and extracts the residuals
  X1 <- cbind(y.tmp,Wy, X.tmp, res)  # concatenates Wy, X.tmp, and the residuals as columns in a new matrix
  colnames(X1)[1:2] <- c(as.character(fml)[2], 'weighted_sums')
  #z <- rep(0, nrow(X1))
  out2 <- glm(formula = fml, data = as.data.frame(X1), family=poisson(link="log"))
  #out2 <- glm(formula = fml, data = as.data.frame(X1), family=gaussian) #for WE we will use a gaussian distribution
  out2 <- -logLik(out2)[[1]]# returns the negative log-likelihood of the model fit by ordinalNet()
}

res <- try(optim(0.5, autoglm_diversity, method="Brent", lower=0, upper=1, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = fml))
res
saveRDS(res, file = 'res_IE.rds')

#we need a version of the model that outputs the BIC

autoglm_diversity <- function (a, y, X, Wsp, Wphy, fml) {
  tmp <- which(!is.na(rowSums(X)))  # identifies rows in X that do not have missing values
  X.tmp <- X[tmp,]  # subsets X to exclude rows with missing values
  y.tmp <- y[tmp]  # subsets y to exclude rows with missing values
  #W <- a*Wsp+(1-a)*Wphy  # calculates weights using values in a, Wsp, Wnb, and Wphy
  #W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
  W <- a*Wsp+(1-a)*Wphy
  W <- W[tmp, tmp]  # subsets W to exclude rows and columns with missing values
  Wy <- W %*% as.numeric(y.tmp)  # calculates a weighted sum of y.tmp using W
  X2 <- W %*% X.tmp  # calculates a weighted sum of X.tmp using W
  X2 <- cbind(X.tmp, X2)  # concatenates X.tmp and X2 as columns in a new matrix
  res <- lm(Wy~X2)$residuals  # fits a linear model of Wy on X2 and extracts the residuals
  X1 <- cbind(y.tmp,Wy, X.tmp, res)  # concatenates Wy, X.tmp, and the residuals as columns in a new matrix
  colnames(X1)[1:2] <- c(as.character(fml)[2], 'weighted_sums')
  #z <- rep(0, nrow(X1))
  out2 <- glm(formula = fml, data = as.data.frame(X1), family=poisson(link="log"))
  #out2 <- glm(formula = fml, data = as.data.frame(X1), family=gaussian) #for WE we will use a gaussian distribution

  out2 <- BIC(out2)# returns the negative log-likelihood of the model fit by ordinalNet()
}

model <- autoglm_diversity(a = res$par, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = fml)
model

autoglm_handler <- function(fml) { autoglm_diversity(a = res$par, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = fml)}
autoglm_handler(fml)


#now lets code up the most complex version of our model 

fml <- formula('island_endemics ~ dist25k + dist_continents + area_km2 + w_mean_growing_season + 
   max_elevation  + w_mean_temperature_seasonality + w_mean_Seasonality_Precipitation + 
   weighted_sums + res')

#we will use this to create all possible model combinations

#we need to generate some models
all_vars <- unlist(strsplit(as.character(fml)[[3]], '+ '))
all_vars <- all_vars[which(all_vars != "+")]
#all_vars <- all_vars[-c(10:11)]
all_vars <- all_vars[-c(8:10)]
all_vars


seq <- 1:length(all_vars)
combs <- lapply(seq, function(x) combn(all_vars, x))
combs


model_combos <- lapply(combs, function(x) as.list(as.data.frame(x)))
length(model_combos)
model_combos <- Reduce(c, model_combos)
names(model_combos) <- NULL
model_combos

for(i in 1:length(model_combos)) {
  model_combos[[i]] <- paste0('island_endemics ~ ', paste0(c(model_combos[[i]]), collapse = ' + ' ),'+ weighted_sums + res')
  model_combos[[i]] <- formula(model_combos[[i]])
  
}
model_combos[[1]]
length(model_combos)

#and we just need to add the intercept only model
model_combos <- c(model_combos, formula(island_endemics ~ weighted_sums + res))

require(pbmcapply)
BICs <- pbmclapply(model_combos, autoglm_handler, mc.cores = 3)
save(BICs, file = 'BIC_IE.rdata')
save(model_combos, file = 'model_combosIE.Rdata')




### FROM HERE KEAGHAN #######

# 
# load('model_combosLR.Rdata')
# load('BIC_LR.rdata')
# 
load('model_combosIE.Rdata')
load('BIC_IE.rdata')
# 
# #load('model_combos_lnWE.Rdata')
# #load('BIC_lnWE.rdata')
# 
model_combos[[2]]

model_combos <- model_combos[order(unlist(BICs), decreasing = F)]
model_combos[[1]]

BICs <- unlist(BICs)[order(unlist(BICs), decreasing = F)]

compute_bic_weights <- function(bic_values) {
  delta_bic <- bic_values - min(bic_values)
  weights <- exp(-0.5 * delta_bic) / sum(exp(-0.5 * delta_bic))
  return(list(delta_bic, weights))
}

bic_weights <- (compute_bic_weights(BICs))

require(dplyr)
model_summaries <- tibble(BIC = BICs, delta = bic_weights[[1]], weights = bic_weights[[2]])
model_summaries
model_matrix <- matrix(ncol = length(all_vars), nrow = 0)
colnames(model_matrix) <- all_vars
model_matrix


for(i in 1:length(model_combos)){
  model_matrix <- rbind(model_matrix, colnames(model_matrix) %in%
                          Reduce(c, strsplit(as.character(model_combos[[i]])[[3]], ' '))
  )
  print(i)
}
head(model_matrix)

row_sums <- rowSums(model_matrix)

model_matrix <- data.frame(apply(model_matrix, 2,function(x) ifelse(x, "X", " ")))
head(model_matrix)
head(model_summaries)
model_matrix <- as_tibble(model_matrix)
model_summaries <- bind_cols(model_summaries, model_matrix,)
model_summaries$model_id <- 1:nrow(model_summaries)
model_summaries$Nparam <- row_sums

model_summaries <- arrange(model_summaries, BIC)
model_summaries$model_id
model_summaries
#write.csv(model_summaries, 'model_summariesLR.csv')
write.csv(model_summaries, 'model_summariesIE.csv')
#write.csv(model_summaries, 'model_summariesWE.csv')

# modelsLR <- read.csv('model_summariesLR.csv')
# head(modelsLR)
# modelsWE <- read.csv('model_summariesWE.csv')
# head(modelsWE)
# modelsEnd <- read.csv('model_summariesEnd.csv')
# head(modelsEnd)







