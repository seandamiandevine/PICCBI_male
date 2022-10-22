rm(list=ls())
library(lme4)
se = function(x) sd(x,na.rm=T)/sqrt(length(x))

# Set constants ----------------------------------------------------------------

freq_list  = read.csv('block_freqs.csv')
conditions = c(0,1)
nblocks    = 16
ntrials    = 800
ntrialperb = ntrials/nblocks
nsubjects  = seq(50,200,by=25)
sizes      = seq(0, 600, by=10)
nsamples   = 100

# Coefficients 
# see https://journals.sagepub.com/doi/full/10.1177/09567976221082941

# og_gammas = c('(Intercept)'  = -1.9,
#               condition_c = 0.08,
#               trial0_c    = -0.62, 
#               size0_c     = 21.21, 
#               'condition_c:trial0_c' = -0.65, 
#               'condition_c:size0_c'  = -0.48, 
#               'trial0_c:size0_c'     = 2.05, 
#               'condition_c:trial0_c:size0_c' = 3.85)

# og_tau  = c('id.(Intercept)' = 1.21, 'id.trial0_c' = 1.53, 'id.trial0_c.(Intercept)'=0)  # assume independent raneffs

# try to weaken effects to test mesi
og_gammas = c('(Intercept)'  = 0,                     # indifference at baseline
              condition_c = 0.01,                     # smaller
              trial0_c    = -0.25,                    # smaller  
              size0_c     = 21.21,                    # unchanged
              'condition_c:trial0_c' = -0.1,          # smaller
              'condition_c:size0_c'  = -0.48,         # unchaned
              'trial0_c:size0_c'     = 2.05,          # unchanged
              'condition_c:trial0_c:size0_c' = 2)     # smaller

# unchanged var structure
og_tau  = c('id.(Intercept)' = 1.21, 'id.trial0_c' = 1.53, 'id.trial0_c.(Intercept)'=0)  # assume independent raneffs


# Simulate ----------------------------------------------------------------

power = data.frame()

for(n in nsubjects) {
  cat('\n\n================================Simulating with N =',n,'================================\n\n')
  out = matrix(NA, 
               nsamples, 
               2+length(og_gammas)+length(og_tau)+length(og_gammas), # n, s, coefs, tau, p
  ) 
  colnames(out) = c('n','s',
                    paste0('b_',names(og_gammas)), 
                    paste0('t_',names(og_tau)), 
                    paste0('p_',names(og_gammas)))
  
  # Design matrix
  df           = expand.grid(trialinblock=1:ntrialperb, block=1:nblocks, id=1:n)
  df$trial     = rep(1:ntrials, n)
  df$condition = ifelse(df$id %% 2==0, 1, 0)
  df$freq      = ifelse(df$condition==0, freq_list[df$block,2], freq_list[df$block,3])
  stimcat      = sapply(df$freq, function(i) sample(c(0,1), size=1, prob=c(1-i,i)))
  df$size      = sapply(stimcat, function(i) ifelse(i==0, 
                                                    sample(sizes[sizes < median(sizes)], size=1), 
                                                    sample(sizes[sizes > median(sizes)], size=1))
                        )
  
  df$trial_c   = df$trial-median(df$trial)
  df$size_c    = df$size-median(df$size)
  df$trial0_c  = df$trial_c/max(df$trial_c)
  df$size0_c   = df$size_c/max(df$size_c)
  df$condition_c = df$condition-.5
  
  # x = tapply(df$size, list(df$block, df$condition), mean)
  # matplot(x, pch=16, lty=1, type='b')
  
  # Simulate nsamples responses
  formula = response ~ condition_c*trial0_c*size0_c + (trial0_c|id)
  p = list(beta=og_gammas, theta = og_tau)
  y = simulate(formula[-2], nsim=nsamples, family='binomial', newparams=p, newdata=df)
 
  # df$response = y[,2]
  # stimbin = cut(df$size, 10 , labels = F)
  # timebin = cut(df$trial, 4, labels = F)
  # mresp = tapply(df$response, list(stimbin, timebin, df$condition), mean)
  # diff  = mresp[,4,] - mresp[,1,]
  # matplot(diff, lty=1, pch=16, type='b', col=c('orange', 'darkred'))
  
  pb = txtProgressBar(1,nsamples)
  for(s in 1:nsamples) {
    
    df$response = y[,s]
    suppressWarnings({
        mod = glmer(formula, data = df, family='binomial', 
                    nAGQ = 0, glmerControl(calc.derivs = F)) # make it go fast
      })
    
    b = as.vector(fixef(mod))
    t = as.vector(c(diag(VarCorr(mod)[[1]]), VarCorr(mod)[[1]][1,2]))
    p = as.vector(summary(mod)$coefficients[,'Pr(>|z|)'])
    out[s,] = c(n,s, b, t, p)
    
    setTxtProgressBar(pb,s)
    
  }
  
  cat('P(Sig)|N =', mean(out[,'p_condition_c:trial0_c:size0_c'] < 0.05, na.rm=T), '\n')
  power = rbind(power, as.data.frame(out))
}

write.csv(power, 'power_nsim=100_mesi.csv')


# Visualize ---------------------------------------------------------------

power = read.csv('power_nsim=100_mesi.csv')
alpha = 0.05
sig   = ifelse(power$p_condition_c.trial0_c.size0_c < alpha, 1, 0)
psig  = tapply(sig, power$n, mean)
sesig = tapply(sig, power$n, se)

pdf('power_curve.pdf', 6, 6)
plot(names(psig), psig, type='b', pch=16, xlab='Sample Size', ylab='P(Significant)',
     main='Power Simulation', ylim=c(0,1))
arrows(as.numeric(names(psig)), psig-sesig,as.numeric(names(psig)), psig+sesig, length=0)
abline(h=0.8, lty=2)
dev.off()

