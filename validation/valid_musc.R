
rm(list=ls())

library(lme4); library(lmerTest)
se = function(x) sd(x,na.rm=T)/sqrt(length(x))

# Load data ---------------------------------------------------------------

dat  = read.csv('raw_ratings_musc.csv')
demo = read.csv('demo_musc.csv')

# Descriptives ------------------------------------------------------------

N = length(unique(dat$subject))

# Clean -------------------------------------------------------------------

dat = dat[dat$attractiveness!='\"',]
dat$size = ifelse(substring(dat$img, 1, 1) == 'T', 
                  300 - as.numeric(substring(dat$img, 2, 4)), 
                  300 + as.numeric(substring(dat$img, 2, 4)))
dat$size0          = dat$size/max(dat$size, na.rm=T)
dat$weight         = as.numeric(dat$weight)
dat$attractiveness = as.numeric(dat$attractiveness)
dat$beauty         = as.numeric(dat$beauty)
dat$harmony        = as.numeric(dat$harmony)

# Weight ----------------------------------------------------------

## Visualize ----
m      = tapply(-dat$weight, dat$size, mean, na.rm=T)
pts    = tapply(-dat$weight, list(dat$subject,dat$size), mean, na.rm=T)
ylimit = c(min(pts), max(pts))
yax    = seq(ylimit[1], ylimit[2], length.out=5)

pdf('figs/musc_weight_ratings.pdf', 6,4)
plot(names(m), m, xlab='Objective Muscularity', ylab='Avg. Muscularity Rating', main='Muscularity Ratings',
     ylim=ylimit, col='red', pch=5, cex=1.5, lwd=2, yaxt='n')
axis(2,at=yax, labels = rev(abs(yax)))
for(i in 1:ncol(pts)) {
  points(rep(names(m)[i],nrow(pts)), pts[,i], pch=16, col=scales::alpha('grey20', .05))
}

dev.off()


## Model ----
dat$neg_weight = -dat$weight

weight_mod = lmer(neg_weight ~ size0 + (size0|subject), data=dat)
summary(weight_mod)
sjPlot::tab_model(weight_mod)

# Attractiveness ----------------------------------------------------------

## Visualize ----
m    = tapply(dat$attractiveness, dat$size, mean, na.rm=T)
pts  = tapply(dat$attractiveness, list(dat$subject,dat$size), mean, na.rm=T)
labs = names(rev(tapply(dat$size, dat$img,mean)))

pdf('figs/musc_attractiveness_ratings.pdf', 6,4)
plot(names(m), m, xlab='Objective Muscularity', ylab='Avg. Attractiveness Rating', main='Attractiveness'
     , ylim=c(0,100), col='red', pch=5, cex=1.5, lwd=2)
for(i in 1:ncol(pts)) {
  points(rep(names(m)[i],nrow(pts)), pts[,i], pch=16, col=scales::alpha('grey20', .05))
}

dev.off()

## Model ----

attractiveness_mod = lmer(attractiveness ~ size0 + (size0|subject), data=dat)
summary(attractiveness_mod)


# Harmony ----------------------------------------------------------

## Visualize ----
m    = tapply(dat$harmony, dat$size, mean, na.rm=T)
pts  = tapply(dat$harmony, list(dat$subject,dat$size), mean, na.rm=T)

pdf('figs/musc_harmony_ratings.pdf', 6,4)
plot(names(m), m, xlab='Objective Muscularity', ylab='Avg. Harmony Rating', main='Harmony',
     ylim=c(0,100), col='red', pch=5, cex=1.5, lwd=2)
for(i in 1:ncol(pts)) {
  points(rep(names(m)[i],nrow(pts)), pts[,i], pch=16, col=scales::alpha('grey20', .05))
}

dev.off()

## Model ----

harmony_mod = lmer(harmony ~ size0 + (size0|subject), data=dat)
summary(harmony_mod)

# Beauty ------------------------------------------------------------------

## Visualize ----
m    = tapply(dat$beauty, dat$size, mean, na.rm=T)
pts  = tapply(dat$beauty, list(dat$subject,dat$size), mean, na.rm=T)

pdf('figs/musc_beauty_ratings.pdf', 6,4)
plot(names(m), m, xlab='Objective Muscularity', ylab='Avg. Beauty Rating', main='Beauty',
     ylim=c(0,100), col='red', pch=5, cex=1.5, lwd=2)
for(i in 1:ncol(pts)) {
  points(rep(names(m)[i],nrow(pts)), pts[,i], pch=16, col=scales::alpha('grey20', .05))
}

dev.off()

## Model ----

beauty_mod = lmer(beauty ~ size0 + (size0|subject), data=dat)
summary(beauty_mod)






