setwd('~/Dropbox/Lindell/Island_Languages/LANGUAGE DATA/')
library(readr)
library(stringr)
library(dplyr)
library(knitr)
library(ggplot2)
library(rworldmap)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(sp)
require(maptools)
require(terra)


#we know that phenome counts is count data, so we should use a Poisson regression
setwd('~/Dropbox/Lindell/Island_Languages/LANGUAGE DATA')
require(nlme)
df <- read.csv('phoneme_data_by_lang.csv')
df$log_area <- df$area

any(df$min_dist25k != df$min_dist11k)
#okay, so it looks like our min distances for 25k and 11k are identifical so we will drop one of these from 
#the analysis

colnames(df)

#we need to scale our our continous predictor variables

#scatter plot by region

regional_scater <- ggplot(df, aes(x = PC_Phoible, y = PC_Creanza, col = macroarea)) +
  geom_point() +
  scale_x_continuous(trans = 'log10', name = 'log10(PC Creanza)') +
  scale_y_continuous(trans = 'log10', name = 'log10(PC Phoible)') +
  scale_color_discrete(name = '') +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey') +
  theme_classic() +
  coord_equal() +
  theme(legend.position = 'bottom')
regional_scater
df
ggsave('Creanza_v_Phoible.pdf', scale = .8)

df[,c(2,3,4,8,9,10,21)] <- scale(df[,c(2,3,4,8,9, 10,21)])
df

#we need to add latitude

lang_map

full_model <- formula('PC_Creanza ~ L1_pop + min_dist11k + min_distCont +
                 area + IE')

glm_fullC <- glm(full_model, data = df, family = 'poisson')
summary(glm_fullC)

#right so it looks like L1_pop, the size of the island and beloning to particular regions is 
#driving a lot of the phoneme count

#Lets look for our other dataset
full_model <- formula('PC_Phoible ~ L1_pop + min_dist11k + min_distCont +
                 area + IE + Bordering_LR')
glm_fullP <- glm(full_model, data = df, family = 'poisson')
summary(glm_fullP)
#very different results, with most of the predictors being significant including 
#whether a language is an island endemic or not

df_clean <- df[,c(1:2, 3:10)]

df_clean
load('~/Dropbox/Lindell/Polysynthesis/For_Xia/Wphy.Rdata')
load('~/Dropbox/Lindell/Polysynthesis/For_Xia/Wsp.Rdata')

#this just makes sure we only include languages we have phylogenetic info for

dim(df_clean)


Wsp <- Wsp[colnames(Wsp) %in% df$LANG_ISO,rownames(Wsp) %in% df$LANG_ISO]
Wphy <- Wphy[colnames(Wphy) %in% df$LANG_ISO,rownames(Wphy) %in% df$LANG_ISO]


df_clean <- df_clean[which(df_clean$LANG_ISO %in% colnames(Wphy) == T),]

fml <- formula('PC_Phoible ~ L1_pop + min_dist11k + min_distCont + area + IE  + res + weighted_sums')

autoglm <- function (a, y, X, Wsp, Wphy, fml) {
  tmp <- which(!is.na(rowSums(X)))  # identifies rows in X that do not have missing values
  X.tmp <- X[tmp,]  # subsets X to exclude rows with missing values
  y.tmp <- y[tmp]  # subsets y to exclude rows with missing values
  #W <- a*Wsp+(1-a)*Wphy  # calculates weights using values in a, Wsp, Wnb, and Wphy
 # W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
  W <- a*Wsp+(1-a)*Wphy
  W <- W[tmp, tmp]  # subsets W to exclude rows and columns with missing values
  Wy <- W %*% as.numeric(y.tmp)  # calculates a weighted sum of y.tmp using W
  X2 <- W %*% X.tmp  # calculates a weighted sum of X.tmp using W
  X2 <- cbind(X.tmp, X2)  # concatenates X.tmp and X2 as columns in a new matrix
  res <- lm(Wy~X2)$residuals  # fits a linear model of Wy on X2 and extracts the residuals
  X1 <- cbind(y.tmp,Wy, X.tmp, res)  # concatenates Wy, X.tmp, and the residuals as columns in a new matrix
  colnames(X1)[1:2] <- c(as.character(fml)[2], 'weighted_sums')
  out2 <- glm(formula = fml, data = as.data.frame(X1), family='poisson')
  summary(out2)
  
  #-out$loglik  # returns the negative log-likelihood of the model fit by ordinalNet()
  out2 <- -logLik(out2)[[1]]
}

X <- df_clean[,-1]
#X <- X[,-c(7:12)]
X$IE <- as.numeric(X$IE)
X$PC_Creanza <- NULL
X <- as.matrix(X)
y <- df_clean$PC_Phoible

res <- try(optim((1), autoglm, method="Brent", lower=c(0), upper=c(1), y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = fml))
res

aPhoible <- res$par

autoglm <- function (a, y, X, Wsp, Wphy, fml) {
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
  out2 <- BIC(out2)# returns the negative log-likelihood of the model fit by ordinalNet()
}

autoglm2 <- function (a, y, X, Wsp, Wphy, fml) {
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
  #out2 <- BIC(out2)# returns the negative log-likelihood of the model fit by ordinalNet()
}

model <- autoglm(a = res$par, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = fml)
model

autoglm_handler <- function(fml) { autoglm(a = res$par, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = fml)}
autoglm_handler(fml)

#now we can make our list of models
all_vars <- unlist(strsplit(as.character(fml)[[3]], '+ '))
all_vars <- all_vars[which(all_vars != "+")]
#all_vars <- all_vars[-c(15:16)]
all_vars <- all_vars[-c(which(all_vars == 'res' | all_vars ==  'weighted_sums'))]
seq <- 1:length(all_vars)
combs <- lapply(seq, function(x) combn(all_vars, x))
combs


model_combos <- lapply(combs, function(x) as.list(as.data.frame(x)))
length(model_combos)
model_combos <- Reduce(c, model_combos)
names(model_combos) <- NULL
model_combos

for(i in 1:length(model_combos)) {
  model_combos[[i]] <- paste0('PC_Phoible ~ ', paste0(c(model_combos[[i]]), collapse = ' + ' ),'+ weighted_sums + res')
  model_combos[[i]] <- formula(model_combos[[i]])
  
}
model_combos[[length(model_combos) + 1]] <- formula('PC_Phoible ~ 1')
length(model_combos)

require(pbmcapply)
BICs <- pbmclapply(model_combos, autoglm_handler, mc.cores = 3)
save(BICs, file = 'BIC_PC_Phoible.Rdata')
save(model_combos, file = 'model_combos_PC_Phoible.Rdata')

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
write.csv(model_summaries, 'model_summaries_PC_Phoible.csv')


#lets make our best fitting model tables 
model_summaries <- model_summaries[which(model_summaries$delta <= 6),]
model_summaries

model_combos <- model_combos[which(model_summaries$delta <= 6)]
model_combos

fitted_models <- lapply(model_combos, function(x) autoglm2(a = aPhoible, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = x))
fitted_models
save(file = 'phoible_best_fitting_models.Rdata', fitted_models)

model_coefficients <- lapply(fitted_models, coef)
model_coefficients

model_confint <- lapply(fitted_models, confint)
model_confint

summaries <- lapply(fitted_models, summary)
summaries[[1]]$coefficients[which(rownames(summaries[[1]]$coefficients) %in%  colnames(model_summaries)),4]
summaries

model_summaries <- as.data.frame(model_summaries)

for(i in 1:length(fitted_models)) {
  model_summaries[i,which(colnames(model_summaries) %in% names(model_coefficients[[i]]))] <- 
    round(model_coefficients[[i]][which(names(model_coefficients[[i]]) %in% colnames(model_summaries))], 3)
  
  model_summaries[i,which(colnames(model_summaries) %in% names(model_coefficients[[i]]))] <- 
    paste0(model_summaries[i,which(colnames(model_summaries) %in% names(model_coefficients[[i]]))], ' (',
           round(model_confint[[i]][,1][which(rownames(model_confint[[i]]) %in% colnames(model_summaries))], 3),
           ' - ',
           round(model_confint[[i]][,2][which(rownames(model_confint[[i]]) %in% colnames(model_summaries))], 3),
           ')')
}

model_summaries$model_id <- NULL
model_summaries

sum <- summary(fitted_models[[1]])
model_summaries$pR2 <- with(sum, 1 - deviance/null.deviance) #get McFadden's psuedo R2

sum

write.csv(model_summaries, 'decent_modelsPhoible.csv')

summaries

### and now we will do it for Craenza

fml <- formula('PC_Creanza ~ L1_pop + min_dist11k + min_distCont + area + IE + res + weighted_sums')

X <- df_clean[,-1]
#X <- X[,-c(7:12)]
X$IE <- as.numeric(X$IE)
X$PC_Phoible <- NULL
X <- as.matrix(X)
y <- df_clean$PC_Creanza

res <- try(optim((1), autoglm, method="Brent", lower=c(0), upper=c(1), y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = fml))
res

aCreanza <- res$par

autoglm <- function (a, y, X, Wsp, Wphy, fml) {
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
  out2 <- BIC(out2)# returns the negative log-likelihood of the model fit by ordinalNet()
}

model <- autoglm(a = res$par, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = fml)
model

autoglm_handler <- function(fml) { autoglm(a = res$par, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = fml)}
autoglm_handler(fml)


model_combos <- lapply(combs, function(x) as.list(as.data.frame(x)))
length(model_combos)
model_combos <- Reduce(c, model_combos)
names(model_combos) <- NULL
model_combos

for(i in 1:length(model_combos)) {
  model_combos[[i]] <- paste0('PC_Creanza ~ ', paste0(c(model_combos[[i]]), collapse = ' + ' ),'+ weighted_sums + res')
  model_combos[[i]] <- formula(model_combos[[i]])
  
}

model_combos[[length(model_combos) +1]] <- formula('PC_Creanza ~ 1')
require(pbmcapply)
BICs <- pbmclapply(model_combos, autoglm_handler, mc.cores = 3)
save(BICs, file = 'BIC_PC_Creanza.Rdata')
save(model_combos, file = 'model_combos_PC_Creanza.Rdata')

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
model_summaries$IE
model_summaries


write.csv(model_summaries, 'model_summaries_PC_Creanza.csv', row.names = F)

#lets make our best fitting model tables 
model_summaries <- model_summaries[which(model_summaries$delta <= 6),]
model_summaries

model_combos <- model_combos[which(model_summaries$delta <= 6)]
model_combos

fitted_models <- lapply(model_combos, function(x) autoglm2(a = aCreanza, y=y, X=X, Wsp=Wsp, Wphy=Wphy, fml = x))
fitted_models
save(file = 'creanza_fitted_models.Rdata', fitted_models)

model_coefficients <- lapply(fitted_models, coef)
model_coefficients

model_confint <- lapply(fitted_models, confint)
model_confint

summaries <- lapply(fitted_models, summary)
summaries[[1]]$coefficients[which(rownames(summaries[[1]]$coefficients) %in%  colnames(model_summaries)),4]
summaries

model_summaries <- as.data.frame(model_summaries)

for(i in 1:length(fitted_models)) {
  model_summaries[i,which(colnames(model_summaries) %in% names(model_coefficients[[i]]))] <- 
    round(model_coefficients[[i]][which(names(model_coefficients[[i]]) %in% colnames(model_summaries))], 3)
  
  model_summaries[i,which(colnames(model_summaries) %in% names(model_coefficients[[i]]))] <- 
    paste0(model_summaries[i,which(colnames(model_summaries) %in% names(model_coefficients[[i]]))], ' (',
           round(model_confint[[i]][,1][which(rownames(model_confint[[i]]) %in% colnames(model_summaries))], 3),
           ' - ',
           round(model_confint[[i]][,2][which(rownames(model_confint[[i]]) %in% colnames(model_summaries))], 3),
           ')')
}

model_summaries$model_id <- NULL
model_summaries

pR2 <- vector()
for(i in 1:length(fitted_models)) {
  sum <- summary(fitted_models[[i]])
  pR2 <- c(pR2, with(sum, 1 - deviance/null.deviance)) #get McFadden's psuedo R2
}
model_summaries$pR2 <- pR2
write.csv(model_summaries, 'decent_modelsCreanza.csv')


