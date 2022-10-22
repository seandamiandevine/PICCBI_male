
rm(list=ls())
library(lme4); library(lmerTest)
se = function(x) sd(x,na.rm=T)/sqrt(length(x))

# Load data ---------------------------------------------------------------

files1 = paste0('data/thin_fat/',list.files('data/thin_fat/'))
files2 = paste0('data/thin_musc/',list.files('data/thin_musc/'))

dat1   = do.call(plyr::rbind.fill, lapply(files1, read.csv ))
dat2   = do.call(plyr::rbind.fill, lapply(files2, read.csv ))

dat1$stim_type = 'TF'
dat2$stim_type = 'TM'

dat = plyr::rbind.fill(dat1, dat2)
rm(dat1,dat2)


## Clean ----
dat$stimulus = dplyr::lag(dat$stimulus)
dat          = dat[dat$stimulus!='',]
dat          = dat[!is.na(dat$trial),]
dat$stimulus = sub(".*/", "", dat$stimulus)
dat$stimulus = sub(".png", "", dat$stimulus)
dat$size     = ifelse(substring(dat$stimulus, 1, 1) == 'T', 
                      300 - as.numeric(substring(dat$stimulus, 2, 4)), 
                      300 +  as.numeric(substring(dat$stimulus, 2, 4)))

## Subset and clean more ----
choice              = dat[,c('subject','stim_type', 'prolID', 'condition', 'age', 'gender', 's1self', 'chosenbody', 's1selfjudge', 's2selfjudge', 's2self', 'stimulus', 'size', 'key_press', 'rt', 'trial', 'block', 'freq')]
choice              = choice[choice$gender=='Male',] # exclude 1 enby
choice              = choice[choice$block!='Practice', ]
choice$hit          = ifelse(choice$key_press==76, 1, 0)
choice$trialinblock = choice$trial
choice$trial        = as.vector(unlist(tapply(choice$subject, choice$subject, function(i) 1:length(i))))
choice$block        = as.numeric(choice$block)
choice$timebin      = cut(choice$trial, 4, labels = F)
choice$stimbin      = cut(choice$size, 10, labels=F)
choice$rt           = as.numeric(choice$rt)

quest  = dat[,c('subject', colnames(dat)[!colnames(dat) %in% colnames(choice)])]

# Demographics ------------------------------------------------------------

N = length(unique(dat$prolID))
table(choice[choice$trial==1,'stim_type'],choice[choice$trial==1,'condition'])

Mage = mean(choice$age)
sage = sd(choice$age)
rage = range(choice$age)

cat('Mean age =', Mage, '( SD =',sage,', range =',rage[1],'-',rage[2],')\n')

z = choice[choice$trial==1,]
tapply(z$age, z$stim_type, mean)
tapply(z$age, z$stim_type, sd)
tapply(z$age, z$stim_type, range)



# Visualize Stim Distribution ---------------------------------------------

cols = c('pink', 'grey80')

pdf('figs/stim_distributions.pdf', 6, 4)

msize = tapply(choice$size, list(choice$block,choice$stim_type), mean)
matplot(1:max(choice$block),msize,type='b', lty=1, col=cols, pch=16, lwd=2, 
        xlab = 'Block Number', ylab='Mean Size of Stimulus')
legend('right',bty='n', lty=1,pch=16,col=cols,title='Stimulus Type',legend=c('Thin/Overweight','Thin/Muscular'))

dev.off()


# Check for backward mapping ----------------------------------------------

pjudge = tapply(choice$hit, list(choice$stimbin, choice$timebin, choice$prolID), mean)
cols   = c('orange', 'darkred')

pdf('figs/per_subject/per_subject.pdf', 6,4)
for(i in unique(choice$prolID)) {

  tmp = pjudge[,c(1,4),i]
  matplot(tmp, type='b', pch=16, xaxt='n', lty=1, col=cols, ylim=c(0,1), ylab=c('P(Judge)'), 
          main=paste0(i,'\n', choice$stim_type[choice$prolID==i][1], '\n',choice$condition[choice$prolID==i][1]))
  axis(1,at=c(1,max(choice$stimbin)), labels=c('\nVery\nThin', '\nVery\nOverweight'), tick = F)
  legend('topleft',bty='n',lty=1,lwd=2,col=cols,legend=c('First 200 Trials', 'Last 200 Trials'), cex=.75)
  
}

dev.off()

# visually inspect file to determine backward mappings
bw = c('s8zncc1gxl2fvaf', 'uyw46bj9hbdkgcf', 'tnel6ukuahn5lgh')

choice$hit = ifelse(choice$subject %in% bw, ifelse(choice$hit==1,0,1), choice$hit)


# Visualize PICC curve ----------------------------------------------------

pjudge  = tapply(choice$hit, list(choice$stimbin, choice$timebin, choice$condition, choice$stim_type), mean)
sejudge = tapply(choice$hit, list(choice$stimbin, choice$timebin,choice$condition, choice$stim_type), se)
cols = c('orange', 'darkred')

pdf('figs/stable_prevalence_curve_TF.pdf', 6, 4)
matplot(pjudge[,c(1,4),1, 'TF'], type='b', lty=1, lwd=2, pch=16, col=cols, ylim=c(0,1),
        xaxt='n',xlab='', ylab='P(Judge Overweight)', main='Equal P(Thin)\nStable')
axis(1,at=c(1,max(choice$stimbin)), labels=c('\nVery\nThin', '\nVery\nOverweight'), tick = F)
legend('topleft',bty='n',lty=1,lwd=2,col=cols,legend=c('First 200 Trials', 'Last 200 Trials'))

dev.off()

pdf('figs/decreasing_prevalence_curve_TF.pdf', 6, 4)
matplot(pjudge[,c(1,4),2,'TF'], type='b', lty=1, lwd=2, pch=16, col=cols, ylim=c(0,1),
        xaxt='n',xlab='', ylab='P(Judge Overweight)', main='Increasing P(Thin)\nChanging')
axis(1,at=c(1,max(choice$stimbin)), labels=c('\nVery\nThin', '\nVery\nOverweight'), tick = F)
legend('topleft',bty='n',lty=1,lwd=2,col=cols,legend=c('First 200 Trials', 'Last 200 Trials'))

dev.off()

pdf('figs/stable_prevalence_curve_TM.pdf', 6, 4)
matplot(pjudge[,c(1,4),1, 'TM'], type='b', lty=1, lwd=2, pch=16, col=cols, ylim=c(0,1),
        xaxt='n',xlab='', ylab='P(Judge Muscular)', main='Equal P(Muscular)\nStable')
axis(1,at=c(1,max(choice$stimbin)), labels=c('\nVery\nThin', '\nVery\nMuscular'), tick = F)
legend('topleft',bty='n',lty=1,lwd=2,col=cols,legend=c('First 200 Trials', 'Last 200 Trials'))

dev.off()

pdf('figs/decreasing_prevalence_curve_TM.pdf', 6, 4)
matplot(pjudge[,c(1,4),2,'TM'], type='b', lty=1, lwd=2, pch=16, col=cols, ylim=c(0,1),
        xaxt='n',xlab='', ylab='P(Judge Muscular)', main='Increasing P(Muscular)\nChanging')
axis(1,at=c(1,max(choice$stimbin)), labels=c('\nVery\nThin', '\nVery\nMuscular'), tick = F)
legend('topleft',bty='n',lty=1,lwd=2,col=cols,legend=c('First 200 Trials', 'Last 200 Trials'))

dev.off()


# Ambiguous stim diff plot ----------------------------------------------------------

idx = floor(median(unique(choice$stimbin)))
dm  = matrix(
  c((pjudge[,4,,'TF'] - pjudge[,1,,'TF'])[idx,],
    (pjudge[,4,,'TM'] - pjudge[,1,,'TM'])[idx,]),
  2,2, byrow=T
  )
dse = matrix(
  c((sejudge[,4,,'TF'] - sejudge[,1,,'TF'])[idx,],
    (sejudge[,4,,'TM'] - sejudge[,1,,'TM'])[idx,]),
  2,2, byrow=T
)
ylimit = range(pretty(c(dm-dse,dm+dse)))

dm = t(dm)
dse = t(dse)
rownames(dm) = c('Stable', 'Changing')
  
pdf('figs/ambiguous_body.pdf', 6, 6)

par(mar=c(5.1, 5.1, 4.1, 2.1))
b = barplot(dm, beside=T, ylim=ylimit,
            names.arg = c('Judge Overweight', 'Judge Muscular'),
            ylab='P(Judge Same Body)\nEnd vs. Start', 
            main='Most Ambiguous Body', 
            legend.text = T, 
            args.legend = list(x='bottomleft', bty='n', title='Prevalence Condition'))

abline(h=0)
arrows(b, dm-dse, b, dm+dse, length=0)

dev.off()

# Model Choice -------------------------------------------------------------------

choice$size_c      = choice$size - median(unique(choice$size))
choice$size_c0     = choice$size_c/max(choice$size_c)
choice$trial_c     = choice$trial - median(choice$trial)
choice$trial_c0    = choice$trial_c/max(choice$trial_c)
choice$condition_c = choice$condition - .5
choice$stim_type_c = as.numeric(choice$stim_type=='TM')-.5

m0 = glmer(hit ~ 1 + (trial_c0|prolID), data=choice, family='binomial', nAGQ = 0)
m1 = glmer(hit ~ size_c0 + (trial_c0|prolID), data=choice, family='binomial', nAGQ = 0)
m2 = glmer(hit ~ trial_c0*size_c0 + (trial_c0|prolID), data=choice, family='binomial', nAGQ = 0)
m3 = glmer(hit ~ condition_c*trial_c0*size_c0 + (trial_c0|prolID), data=choice, family='binomial', nAGQ = 0)
m4 = glmer(hit ~ stim_type_c*condition_c*trial_c0*size_c0 + (trial_c0|prolID), data=choice, family='binomial', nAGQ = 0)

capture.output(anova(m0, m1, m2, m3, m4), file='out/choice_lrt.txt')
saveRDS(m4, file='out/choice_model.rds')
sjPlot::tab_model(m4, transform = NULL, file='out/choice_model.html')

submod1 = glmer(hit ~ condition_c*trial_c0*size_c0 + (trial_c0|prolID), 
                data=choice[choice$stim_type=='TF',], family='binomial', nAGQ = 0)
saveRDS(submod1, file='out/choice_submodel1.rds')
sjPlot::tab_model(submod1, transform = NULL, file='out/choice_submodel1.html')

submod2 = glmer(hit ~ condition_c*trial_c0*size_c0 + (trial_c0|prolID), 
                data=choice[choice$stim_type=='TM',], family='binomial', nAGQ = 0)
saveRDS(submod2, file='out/choice_submodel2.rds')
sjPlot::tab_model(submod2, transform = NULL, file='out/choice_submodel2.html')


# Self-judgements ---------------------------------------------------------

## Self judgement ----
self = reshape2::melt(choice[choice$trial==1,c('prolID','stim_type','condition','s1self','s2self')], 
                      id.vars=c('prolID','stim_type','condition'))

mself  = tapply(self$value, list(self$variable, self$condition, self$stim_type), mean)
seself = tapply(self$value, list(self$variable, self$condition, self$stim_type), se)

# Thin/Overweight
pdf('figs/self_judgement_TF.pdf', 6, 6)

b = barplot(mself[,,'TF'], beside=T, ylim=c(0,1), 
            xlab='Prevalence Condition',
            ylab='Size of Chosen Body',
            names.arg = c('Stable','Changing'), 
            main='Thin/Overweight Judgements', 
            legend.text = T, 
            args.legend = list(bty='n', legend=c('Start of Task', 'End of Task'), cex=1.5))
arrows(b, mself[,,'TF']-seself[,,'TF'],b, mself[,,'TF']+seself[,,'TF'], length=0)

dev.off()

# Thin/Muscular
pdf('figs/self_judgement_TM.pdf', 6, 6)

b = barplot(mself[,,'TM'], beside=T, ylim=c(0,1), 
            xlab='Prevalence Condition',
            ylab='Size of Chosen Body',
            names.arg = c('Stable','Changing'), 
            main='Thin/Muscular Judgements', 
            legend.text = T, 
            args.legend = list(bty='n', legend=c('Start of Task', 'End of Task'), cex=1.5))
arrows(b, mself[,,'TM']-seself[,,'TM'],b, mself[,,'TM']+seself[,,'TM'], length=0)

dev.off()

### Model ----
self$condition_c = self$condition-.5
self$stim_type_c = ifelse(self$stim_type=='TF', -.5,.5)
self$time_c      = ifelse(self$variable=='s1self',-.5,.5)

self_judge_mod = lmer(value ~ condition_c*stim_type_c*time_c + (1|prolID), data=self)
summary(self_judge_mod)
sjPlot::tab_model(self_judge_mod, transform = NULL, file='out/self_judge_mod.html')

### Individual differences ----
conds = tapply(dat$condition, dat$prolID, function(x) x[1])

#### TF ----
tmp   = ranef(submod1)$prolID

r   = tmp$trial_c0[rownames(tmp) %in% names(conds[conds==1])]

tmp   = self[self$stim_type=='TF' & self$condition==1,]
mself = tapply(tmp$value, list(tmp$prolID, tmp$variable), mean, na.rm=T)
diff  = mself[,2]  - mself[,1]

pdf('figs/self_judgement_ind_diff_TF.pdf', 6, 6)

par(mar=c(5.1, 6.1, 4.1, 2.1))

plot(r, diff, pch=3, col=scales::alpha('black', 0.25),
     xlab='', xaxt='n',
     ylab='', yaxt='n',
     main='Thin/Overweight\nChanging Prevalence')
axis(1, at=c(min(r), 0, max(r)), 
     labels = c('Reverse Effect', 'No Effect', 'Strong Effect'), 
     cex.axis=0.7, las=1)
axis(2, at=c(min(diff), 0, max(diff)), 
     labels = c('Thinner', 'Same', 'Bigger'), 
     cex.axis=0.7, las=2)
title(ylab="Judgement at End", line=4, cex.lab=1)


abline(lm(diff~r), col='red', lwd=2)
r_test = cor.test(r, diff)
legend('topright', bty='n', paste0('r = ',round(r_test$estimate,2), '\np = ',round(r_test$p.value,2)))

dev.off()

#### TM ----
tmp   = ranef(submod2)$prolID

r   = tmp$trial_c0[rownames(tmp) %in% names(conds[conds==1])]

tmp   = self[self$stim_type=='TM' & self$condition==1,]
mself = tapply(tmp$value, list(tmp$prolID, tmp$variable), mean, na.rm=T)
diff  = mself[,2]  - mself[,1]

pdf('figs/self_judgement_ind_diff_TM.pdf', 6, 6)

par(mar=c(5.1, 6.1, 4.1, 2.1))

plot(r, diff, pch=3, col=scales::alpha('black', 0.25),
     xlab='', xaxt='n',
     ylab='', yaxt='n',
     main='Thin/Muscular\nChanging Prevalence')
axis(1, at=c(min(r), 0, max(r)), 
     labels = c('Strong Effect', 'No Effect', 'Reverse Effect'), 
     cex.axis=0.7, las=1)
axis(2, at=c(min(diff), 0, max(diff)), 
     labels = c('Thinner', 'Same', 'Bigger'), 
     cex.axis=0.7, las=2)
title(ylab="Judgement at End", line=4, cex.lab=1)


abline(lm(diff~r), col='red', lwd=2)
r_test = cor.test(r, diff)
legend('topright', bty='n', paste0('r = ',round(r_test$estimate,2), '\np = ',round(r_test$p.value,2)))

dev.off()

## Self-concept ----

self = reshape2::melt(choice[choice$trial==1,c('prolID','stim_type','condition','s1selfjudge','s2selfjudge')], 
                      id.vars=c('prolID','stim_type','condition'))

self$resp = ifelse(self$value==65, 0, 1)

mself  = tapply(self$resp, list(self$variable, self$condition, self$stim_type), mean)
seself = tapply(self$resp, list(self$variable, self$condition, self$stim_type), se)

# Thin/Overweight
pdf('figs/self_concept_TF.pdf', 6, 6)

b = barplot(mself[,,'TF'], beside=T, ylim=c(0,1), 
            xlab='Prevalence Condition',
            ylab='P(Judge Chosen Body Overweight)',
            names.arg = c('Stable','Changing'), 
            main='Thin/Overweight Judgements', 
            legend.text = T, 
            args.legend = list(bty='n', legend=c('Start of Task', 'End of Task'), cex=1.5))
arrows(b, mself[,,'TF']-seself[,,'TF'],b, mself[,,'TF']+seself[,,'TF'], length=0)

dev.off()

# Thin/Muscular
pdf('figs/self_concept_TM.pdf', 6, 6)

b = barplot(mself[,,'TM'], beside=T, ylim=c(0,1), 
            xlab='Prevalence Condition',
            ylab='P(Judge Chosen Body Muscular)',
            names.arg = c('Stable','Changing'), 
            main='Thin/Muscular Judgements', 
            legend.text = T, 
            args.legend = list(bty='n', legend=c('Start of Task', 'End of Task'), cex=1.5))
arrows(b, mself[,,'TM']-seself[,,'TM'],b, mself[,,'TM']+seself[,,'TM'], length=0)

dev.off()

### Model ----
self$condition_c = self$condition-.5
self$stim_type_c = ifelse(self$stim_type=='TF', -.5,.5)
self$time_c      = ifelse(self$variable=='s1selfjudge',-.5,.5)

self_concept_mod = glmer(resp ~ condition_c*stim_type_c*time_c + (1|prolID), data=self, family='binomial')
summary(self_concept_mod)

sjPlot::tab_model(self_concept_mod, transform = NULL, file='out/self_concept_mod.html')

#### logistic ordinal regression ------
s1 = as.numeric(choice$s1selfjudge!=65)
s2 = as.numeric(choice$s2selfjudge!=65)

choice$selfjudge_diff = s1-s2

self2 = reshape2::melt(choice[choice$trial==1,c('prolID','stim_type','condition','selfjudge_diff')], 
                      id.vars=c('prolID','stim_type','condition'))
self2$value = as.factor(self2$value)
self2$condition_c = self2$condition-.5
self2$stim_type_c = ifelse(self2$stim_type=='TF', -.5,.5)

cat.mord0 = MASS::polr(value ~ 1, data = self2, Hess = T)
cat.mord1 = MASS::polr(value ~ condition_c, data = self2, Hess = T)
cat.mord2 = MASS::polr(value ~ condition_c + stim_type_c, data = self2, Hess = T)
cat.mord3 = MASS::polr(value ~ condition_c * stim_type_c, data = self2, Hess = T)

cat.mord3.pvals = pnorm(abs(coef(summary(cat.mord3))[, "t value"]),lower.tail = FALSE)*2
cbind(coef(summary(cat.mord3)), "p" = cat.mord3.pvals)
confint(cat.mord3)

mod_comp = anova(cat.mord0, cat.mord1, cat.mord2, cat.mord3)

sjPlot::tab_model(cat.mord3, transform = NULL)

### Individual differences ----

conds = tapply(dat$condition, dat$prolID, function(x) x[1])

#### TF ----
tmp1 = ranef(submod1)$prolID
tmp1 = tmp1[rownames(tmp1) %in% names(conds[conds==1]),]

tmp2   = self[self$stim_type=='TF' & self$condition==1,]
mself = tapply(tmp2$resp, list(tmp2$prolID, tmp2$variable), mean, na.rm=T)
diff  = mself[,1]  - mself[,2]

tmp1$diff = diff

mr     = tapply(tmp1$trial_c0, tmp1$diff, mean)
ser    = tapply(tmp1$trial_c0, tmp1$diff, se)
ylimit = range(pretty(c(mr-ser, mr+ser)))

pdf('figs/self_concept_ind_diff_TF.pdf', 6, 6)

par(mar=c(5.1, 6.1, 4.1, 2.1))
b = barplot(mr, ylim=ylimit,ylab='', yaxt='n', 
            main='Thin/Overweight\nChanging Prevalence Condition',
            names.arg = c('Judge Thin at Start\nJudge Overweight at End', 
                          'No Change', 
                          'Judge Overweight at Start\nJudge Thin at End'), 
            cex.names = .75)
axis(2, at=c(min(mr), 0, max(mr)), 
     labels = c('Reverse Effect', 'No Effect', 'Strong Effect'), 
     cex.axis=0.7, las=1)
arrows(b, mr-ser, b, mr+ser, length=0)
abline(h=0)

##### model with ordinal logit regression ----
# https://www.analyticsvidhya.com/blog/2016/02/multinomial-ordinal-logistic-regression/#:~:text=Ordinal%20regression%20is%20used%20to,one%20or%20more%20independent%20variables.

tmp1$diff_f = factor(tmp1$diff, levels = c(-1,0,1))

cat.mord0 = MASS::polr(diff_f ~ 1, data = tmp1, Hess = T)
cat.mord1 = MASS::polr(diff_f ~ trial_c0, data = tmp1, Hess = T)
cat.mord1.pvals = pnorm(abs(coef(summary(cat.mord1))[, "t value"]),lower.tail = FALSE)*2
cbind(coef(summary(cat.mord1)), "p" = cat.mord1.pvals)

confint(cat.mord1)

mod_comp = anova(cat.mord0, cat.mord1)

legend('bottomleft', bty='n', legend=paste0('p = ',round(mod_comp$`Pr(Chi)`[2],2)))

dev.off()


#### TM ----
tmp1 = ranef(submod2)$prolID
tmp1 = tmp1[rownames(tmp1) %in% names(conds[conds==1]),]

tmp2   = self[self$stim_type=='TM' & self$condition==1,]
mself = tapply(tmp2$resp, list(tmp2$prolID, tmp2$variable), mean, na.rm=T)
diff  = mself[,1]  - mself[,2]

tmp1$diff = diff

mr     = tapply(tmp1$trial_c0, tmp1$diff, mean)
ser    = tapply(tmp1$trial_c0, tmp1$diff, se)
ylimit = range(pretty(c(mr-ser, mr+ser)))

pdf('figs/self_concept_ind_diff_TM.pdf', 6, 6)

par(mar=c(5.1, 6.1, 4.1, 2.1))
b = barplot(mr, ylim=ylimit,ylab='', yaxt='n', 
            main='Thin/Muscular\nChanging Prevalence Condition',
            names.arg = c('Judge Thin at Start\nJudge Muscular at End', 
                          'No Change', 
                          'Judge Muscular at Start\nJudge Thin at End'), 
            cex.names = .75)
axis(2, at=c(min(mr), 0, max(mr)), 
     labels = c('Strong Effect', 'No Effect', 'Reverse Effect'), 
     cex.axis=0.7, las=1)
arrows(b, mr-ser, b, mr+ser, length=0)
abline(h=0)

##### model with ordinal logit regression ----
# https://www.analyticsvidhya.com/blog/2016/02/multinomial-ordinal-logistic-regression/#:~:text=Ordinal%20regression%20is%20used%20to,one%20or%20more%20independent%20variables.

tmp1$diff_f = factor(tmp1$diff, levels = c(-1,0,1))

cat.mord0 = MASS::polr(diff_f ~ 1, data = tmp1, Hess = T)
cat.mord1 = MASS::polr(diff_f ~ trial_c0, data = tmp1, Hess = T)
cat.mord1.pvals = pnorm(abs(coef(summary(cat.mord1))[, "t value"]),lower.tail = FALSE)*2
cbind(coef(summary(cat.mord1)), "p" = cat.mord1.pvals)

confint(cat.mord1)

mod_comp = anova(cat.mord0, cat.mord1)


legend('bottomleft', bty='n', legend=paste0('p = ',round(mod_comp$`Pr(Chi)`[2],2)))

dev.off()


# RT effects --------------------------------------------------------------
## Clean  ----
cat('P(RT < 50 ms) =', mean(choice$rt < 10, na.rm=T), '\nP(RT > 3000 ms)=',mean(choice$rt > 3000, na.rm=T),'\n')
choice$rt[choice$rt < 50]   = NA
choice$rt[choice$rt > 3000] = NA

## Group plot ----
mrt    = tapply(choice$rt, list(choice$condition, choice$stim_type), mean, na.rm=T)
sert   = tapply(choice$rt, list(choice$condition, choice$stim_type), se)
ylimit = range(pretty(c(mrt-sert-50, mrt+sert+50)))

pdf('figs/rt_conditions.pdf', 6, 6)

b = barplot(mrt, beside=T, names.arg = c('Overweight Judgements', 'Mucularity Judgements'), ylab='Mean RT', 
            ylim=ylimit, xpd=F, 
            legend.text = T, args.legend = list(x='topleft',bty='n', title='Prevalence Condition',legend=c('Stable','Changing')))
arrows(b, mrt-sert, b, mrt+sert, length=0)

dev.off()

## Conditional "accuracy" plot ----
mrt    = tapply(choice$rt, list(choice$stimbin, choice$condition, choice$stim_type), mean, na.rm=T)
sert   = tapply(choice$rt, list(choice$stimbin, choice$condition, choice$stim_type), se)
ylimit = range(pretty(c(mrt-sert-50, mrt+sert+50)))

# T/F
pdf('figs/conditional_rt_TF.pdf', 6, 6)

cols = c('grey80', 'grey40')
matplot(mrt[,,1], xlab='Model Size', ylab='Mean RT', main='Thin/Overweight\nJudgements',
        ylim=ylimit,pch=16, col=cols, type='b', lty=1, lwd=2, xaxt='n')
axis(1, at=c(1,max(choice$stimbin)), labels = c('Very\nThin', 'Very\nOverweight'))
arrows(1:nrow(mrt[,,1]), mrt[,1,1]-sert[,1,1], 1:nrow(mrt[,,1]), mrt[,1,1]+sert[,2,1], col=cols[1], length=0) 
arrows(1:nrow(mrt[,,1]), mrt[,2,1]-sert[,2,1], 1:nrow(mrt[,,1]), mrt[,2,1]+sert[,2,1], col=cols[2], length=0) 
legend('bottomleft',bty='n',col=cols,lty=1,pch=16,legend=c('Stable','Changing'), title='Prevalence Condition')

dev.off()

# T/M
pdf('figs/conditional_rt_TM.pdf', 6, 6)

cols = c('grey80', 'grey40')
matplot(mrt[,,2], xlab='Model Size', ylab='Mean RT', main='Thin/Muscular\nJudgements',
        ylim=ylimit,pch=16, col=cols, type='b', lty=1, lwd=2, xaxt='n')
axis(1, at=c(1,max(choice$stimbin)), labels = c('Very\nThin', 'Very\nMuscular'))
arrows(1:nrow(mrt[,,2]), mrt[,1,2]-sert[,1,2], 1:nrow(mrt[,,2]), mrt[,1,2]+sert[,2,2], col=cols[1], length=0) 
arrows(1:nrow(mrt[,,2]), mrt[,2,2]-sert[,2,2], 1:nrow(mrt[,,2]), mrt[,2,2]+sert[,2,2], col=cols[2], length=0) 
# legend('bottomleft',bty='n',col=cols,lty=1,pch=16,legend=c('Stable','Changing'), title='Prevalence Condition')

dev.off()

## Model ----
choice$log_rt   = log(choice$rt)
choice$size_c02 = choice$size_c0^2

rt_mod  = lmer(log_rt ~ stim_type_c * condition_c + (1|prolID), data=choice)

summary(rt_mod)
sjPlot::tab_model(rt_mod, transform = NULL, file='out/rt_mod.html')


# Individual Differences --------------------------------------------------

## Difference score ----
idx = floor(median(unique(choice$stimbin)))
tmp = choice[choice$stimbin==idx & choice$condition==1 ,]

pjudge  = tapply(tmp$hit, list(tmp$prolID, tmp$timebin, tmp$stim_type), mean)
sejudge = tapply(tmp$hit, list(tmp$prolID, tmp$timebin, tmp$stim_type), se)

mdiffTF = pjudge[,4,'TF'] - pjudge[,1,'TF']
mdiffTM = pjudge[,4,'TM'] - pjudge[,1,'TM']

pdf('figs/ind_diff_diffscore.pdf', 4, 4)

par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(mdiffTF, mdiffTM, 
     xlab='P(Judge Same Body Muscular)\nEnd vs. Start',
     ylab='P(Judge Same Body Overweight)\nEnd vs. Start', 
     main='Changing Prevalence\nCondition Only')

dev.off()

## Random effects ----

conds = tapply(dat$condition, dat$prolID, function(x) x[1])

r_tf = ranef(submod1)$prolID['trial_c0']
r_tm = ranef(submod2)$prolID['trial_c0']

r_tf$prolID = rownames(r_tf)
r_tm$prolID = rownames(r_tm)

r = merge(r_tf, r_tm, by = 'prolID')
r$condition = conds[r$prolID]

cols = c('orange', 'darkred')

pdf('figs/ind_diff_raneff.pdf', 4, 4)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(r[r$condition==0,2], r[r$condition==0,3],
     xlab='Random Effect of Trial\nThin/Overweight', 
     ylab='Random Effect of Trial\nThin/Muscular', 
     pch=16, col=cols[1], xlim=range(r[,2]), ylim=range(r[,3]))

points(r[r$condition==1,2], r[r$condition==1,3],pch=16,col=cols[2])
legend('topright',bty='n',pch=16,col=cols,title='Prevalence Condition',legend=c('Stable','Changing'))

dev.off()


# Save image --------------------------------------------------------------

save.image('main.Rdata')
write.csv(choice, 'models/hddm/choice_data_clean.csv')
