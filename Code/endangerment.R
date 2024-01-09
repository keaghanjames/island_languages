library(ordinalNet)
library(ivreg)
library(nloptr)

setwd("~/Dropbox/Lindell/Polysynthesis/plotting language trees")
load("bestmodel")

updated <- read.csv('~/Dropbox/Lindell/Island_Languages/languoid_data_for_analysis.csv')
which((data$id_ISO_lang == updated$ISO) == F)

data$island <- as.numeric(updated$islands11km)

fit <- lm(lang_L1.POP_lang_tr~island,data=data,na.action=na.exclude)
summary(fit)


#accounting for spatial and phylogenetic autocorrelation
X <- data$island
y <- data$lang_L1.POP_lang_tr
Wspy <- Wsp%*%as.numeric(y)
Wnby <- Wnb%*%as.numeric(y)
Wphyy <- Wphy%*%as.numeric(y)
WspX <- Wsp%*%X
WnbX <- Wnb%*%X
WphyX <- Wphy%*%X



#We need to add inteaction terms for each region
data$island_Western_Africa <- data$island*data$Western_Africa
data$island_Oceania<- data$island*data$Oceania
data$island_Europe<- data$island*data$Europe
data$island_Southern_Asia<- data$island*data$Southern_Asia
data$island_South_America<- data$island*data$South_America
data$island_Northern_America<- data$island*data$Northern_America
data$island_Africa<- data$island*data$Africa
data$island_South_Eastern_Asia<- data$island*data$`South-Eastern_Asia`
data$island_Arab<- data$island*data$Arab
data$island_Asia<- data$island*data$Asia
data$island_Central_America <- data$island*data$Central_America
data$island_Australia_and_New_Zealand <- data$island*data$Australia_and_New_Zealand


ncol(data)


autoord.lasso <- function (y,X,W) {
        Wy <- W%*%as.numeric(y)
        X2 <- W%*%X
        X2 <- cbind(X,X2)
        res <- lm(Wy~X2)$residuals
        X1 <- cbind(Wy,X,res)
        ordinalNetCV(X1,y,nFolds=10,alpha=0.5,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,tuneMethod="cvLoglik",maxiterOut=1000)
}

#X <- data[,c(model,colnames(data)[2808:2820])] #model is the best model we had for language endangerment + island endemism status + island region interactions
X <- data[,c(model,colnames(data)[2808])] #model is the best model we had for language endangerment + island endemism status

X <- as.matrix(X)

W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
tmp <- which(!is.na(rowSums(X)))
y.tmp <- y2[tmp]
X.tmp <- X[tmp,]
W <- W[tmp,tmp]
best.poly.model <- autoord.lasso(y=y.tmp,X=X.tmp,W=W)
best.poly.model$fit$coefs[,"island"]
best.poly.model$fit$coefs[20,]
data$island


best.poly.model$fit$coefs[20,]


saveRDS(best.poly.model, '~/Dropbox/Lindell/Island_Languages/best.island.model.rds')
best.poly.model$fit$coefs[20,]
write.csv(best.poly.model$fit$coefs, '~/Dropbox/Lindell/Island_Languages/endangerment_w_island_coefficients.csv')

summary(best.poly.model)
best.poly.model$fit$loglik[20]

best.poly.model$fit$bic
best.model$fit$bic

test_stat <- 2 * (best.model$fit$loglik[20] - best.poly.model$fit$loglik[20])
p_value <- pchisq(test_stat, 1, lower.tail = FALSE)
test_stat
p_value



#endangerment figure
data_old <- data
data <- data_old


data$EGIDS_tr <- as.factor(data$EGIDS_tr)
data$island <- ifelse(data$island == 1, 'Island Endemic', 'Mainland Language')


colnames(data)


data <- data[,c('id_ISO_lang', 'island', 'EGIDS_tr')]
require(dplyr)
prop_end <- data %>% 
  group_by(island, EGIDS_tr, .drop = F) %>%
  summarise(langs = length(id_ISO_lang))

prop_end

sum_isl <- sum(prop_end$langs[which(prop_end$island == 'Island Endemic')])
sum_isl

sum_main <- sum(prop_end$langs[which(prop_end$island == 'Mainland Language')])
sum_main

prop_end$proportion <- ifelse(prop_end$island == 'Island Endemic', prop_end$langs/sum_isl, prop_end$langs/sum_main)

prop_end

require(ggplot2)

# Stacked + percent
endangered_prop <- ggplot(prop_end, aes(fill=island, y=proportion, x=EGIDS_tr)) + 
  geom_bar(position="dodge", stat="identity") +
  annotate("rect", xmin = 0.3, xmax = 6.55, ymin = 0, ymax = 0.35,
           alpha = .7,fill = "white") +
  ylab('Proportion of languages') +
  xlab('EGIDS') +
  scale_fill_manual(values = c('lightblue', 'grey'), name = "") +
 # geom_text(label = 'Not threatened', x = 2, y = .34, fontface = 'italic',
  #          size = 5) +
  theme_classic() +
  theme(legend.position = 'bottom')
endangered_prop
ggsave('~/Dropbox/Lindell/Island_Languages/Figures/proportion_endangered.pdf', scale = .6)

save(endangered_prop, file = '~/Dropbox/Lindell/Island_Languages/Figures/figure3.Rdata')




