import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy.stats import gaussian_kde
import hddm
import os
from tqdm import tqdm

os.chdir('/Users/sean/documents/PICCBI_male/analysis/models/hddm')

plt.style.use('grayscale')
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams['figure.figsize'] = (10.4, 6.8)
matplotlib.rcParams['lines.linewidth'] = 2.5
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=['0', '0.5']) 

pd.options.mode.chained_assignment = None  # ignore chaining warning; default='warn'               

# ****************************************************************************
# *                              Load and Clean                              *
# ****************************************************************************

dat = pd.read_csv('choice_data_clean.csv')
dat = dat.rename(columns={"prolID":'subj_idx', "hit":'response'}) 
dat = dat[~dat.rt.isna()]
dat['rt'] = dat['rt']/1000 # to seconds

# ****************************************************************************
# *                                 Fit HDDM                                 *
# ****************************************************************************

dat['stim_cond'] = dat.stim_type.astype(str) + dat.condition.astype(str)

# basic_ddm = hddm.HDDM(dat, depends_on={'a':'stim_cond', 'v':'stim_cond', 't':'stim_cond'})
# basic_ddm.find_starting_values()
# basic_ddm.sample(3000, burn=1000, db='pickle')
# basic_ddm.save(f'fits/fit_basic')
# ddm_summary = basic_ddm.gen_stats()
# ddm_summary.to_csv('stats/basic_ddm_summary.csv')

regs = ['a ~ trial_c0 * size_c0', 'v ~ trial_c0 * size_c0', 't ~ trial_c0 * size_c0']
ddm  = hddm.models.HDDMRegressor(dat, regs, bias=True, depends_on={'a':'stim_cond', 'v':'stim_cond', 't':'stim_cond'})
ddm.find_starting_values()
ddm.sample(3000,burn=1000,db='pickle')
ddm.save(f'fits/fit_depends_on')

# ddm = hddm.load('fits/fit_depends_on')

# ****************************************************************************
# *                         Add p-values and Geweke's                        *
# ****************************************************************************

ddm_summary = ddm.gen_stats()

fig    = plt.figure(1, figsize=[30,30], tight_layout=True)
pars   = ddm_summary[~ddm_summary.index.str.contains('Intercept_')].index.to_list()
traces = ddm.nodes_db.node[[i for i in pars]]

ddm_summary['p_err']  = ddm_summary['mc err']/ddm_summary['std']
ddm_summary['geweke'] = np.nan
ddm_summary['P']      = np.nan

for i,p in enumerate(traces):
    trace = p.trace()

    # ax = fig.add_subplot(round(len(pars)/9), 9, i+1)
    # ax.plot(trace)
    # ax.set(xlabel='', ylabel='', title=pars[i])

    # geweke
    m1,m2 = np.mean(trace[:int(trace.shape[0]*.5)]), np.mean(trace[int(trace.shape[0]*.9):])
    v1,v2 = np.var(trace[:int(trace.shape[0]*.5)]), np.var(trace[int(trace.shape[0]*.9):])
    z = (m2-m1)/np.sqrt(v2+v1)

    # p-value
    p_lt_0 = np.mean(trace < 0)
    p      = p_lt_0 if ddm_summary.loc[pars[i]]['mean'] < 0 else 1-p_lt_0

    # save
    ddm_summary.loc[pars[i], 'geweke'] = z
    ddm_summary.loc[pars[i], 'P'] = 1-p # to match frequentist interpretation

# fig.savefig('figs/traceplots_fullfit.png')
ddm_summary.to_csv('stats/summary_depends_on_all.csv')
tmp = ddm_summary[~ddm_summary.index.str.contains('Intercept_')]
tmp[~tmp.index.str.contains('_subj')].to_csv('stats/summary_depends_on_main.csv')


# ****************************************************************************
# *                              Error bar plot                              *
# ****************************************************************************

titles = {'a':'Decision Threshold', 'v':'Drift Rate', 't':'Non-Decision Time'}

fig,ax = plt.subplots(1,3, figsize=[16,6], constrained_layout=True)
color  = ['orange','darkred']

for i,p in enumerate(['a','v','t']):
    tmp  = ddm_summary[ddm_summary.index.str[0]==p]
    tmp  = tmp[tmp.index.str.contains('Intercept\(')]
    idx  = tmp.index.to_list()
    cats = [s[s.find("(")+1:s.find(")")] for s in idx]
    tmp['stim_cond'] = [i[:-1] for i in cats]
    tmp['prev_cond'] = [int(i[-1]) for i in cats]

    tmp['stim_cond'] = np.where(tmp.stim_cond=='TF', 'Thin/\nOverweight','Thin/\nMuscular')
    tmp['prev_cond'] = np.where(tmp.prev_cond==1,'Changing','Stable')
       
    tmp['lb'] = tmp['mean'] - tmp['2.5q']
    tmp['ub'] = tmp['97.5q'] - tmp['mean']
    ci = tmp[['lb','ub']].T.to_numpy()

   	sns.pointplot(data=tmp, x='stim_cond', y='mean', hue='prev_cond', linestyles='', ci=None, 
   		ax=ax[i],palette=color, ylabel='', xlabel='')
    ax[i].errorbar(x=[0,0,1,1], y=tmp['mean'], yerr=ci, fmt='none', c=color)
    ax[i].set(title=titles[p], xlabel='', ylabel='')

    if i ==0:
    	ax[i].set(ylabel='Parameter Value')
    	ax[i].legend(title='Prevalence Condition', loc='lower left', prop={'size':14})
    else:
    	ax[i].legend([],[], frameon=False)


plt.show()
fig.savefig('figs/summary_point_all.png')
plt.close()


# ****************************************************************************
# *                        P-values and Density Plots                        *
# ****************************************************************************

titles = {'a':'Decision Threshold', 'v':'Drift Rate', 't':'Non-Decision Time'}
pars   = ddm_summary[ddm_summary.index.str.contains('Intercept\(')].index.to_list()
traces = ddm.nodes_db.node[[i for i in pars]]

diffs  = {'comp':[],'b':[],'ci_lb':[],'ci_ub':[],'p':[]}
stim_cols = ['darkgreen','darkblue']
cond_cols = ['orange','darkred']

for i,p in enumerate(['a','v','t']):

	tmp = traces[[j for j in pars if f'{p}_' in j]]

	# compare stim conds
    fig, ax = plt.subplots(1)

	tf = np.concatenate([j.trace() for j in tmp[tmp.index.str.contains('TF')]])
	tm = np.concatenate([j.trace() for j in tmp[tmp.index.str.contains('TM')]])
    diff = tf-tm

    diffs['comp'].append(f'{p}_TF_TM')
    diffs['b'].append(diff.mean())
    lb,ub  = np.quantile(diff, [.025,.975])
    diffs['ci_lb'].append(lb)
    diffs['ci_ub'].append(ub)
    diffs['p'].append((diff>0).mean())


    xaxis    = np.linspace(np.concatenate([tf,tm]).min(),np.concatenate([tf,tm]).max(), 500)
    dens_tf  = gaussian_kde(tf)(xaxis)
    dens_tm  = gaussian_kde(tm)(xaxis)
    dens_tf0 = dens_tf/dens_tf.max()
    dens_tm0 = dens_tm/dens_tm.max()

    ax.plot(xaxis,dens_tf0, label='Thin/Overweight', c=stim_cols[0])
    ax.plot(xaxis,dens_tm0, label='Thin/Muscular', c=stim_cols[1])
    ax.legend(title='Stimulus Type', loc='upper left')
    ax.set(ylim=[0, 1.5], ylabel='Normalized Density', title=titles[p], yticks=[0,0.5,1])
    ax.fill_between(xaxis, dens_tf0, color=stim_cols[0], alpha=0.35)
    ax.fill_between(xaxis, dens_tm0, color=stim_cols[1], alpha=0.35)

    fig.savefig(f'figs/TF_vs_TM_density_{p}.png')
    plt.close()

	# compare prev conds
    fig, ax = plt.subplots(1)

	t0 = np.concatenate([j.trace() for j in tmp[tmp.index.str.contains('0')]])
	t1 = np.concatenate([j.trace() for j in tmp[tmp.index.str.contains('1')]])
    diff = t0-t1

    diffs['comp'].append(f'{p}_S_C')
    diffs['b'].append(diff.mean())
    lb,ub  = np.quantile(diff, [.025,.975])
    diffs['ci_lb'].append(lb)
    diffs['ci_ub'].append(ub)
    diffs['p'].append((diff>0).mean())

    xaxis    = np.linspace(np.concatenate([t0,t1]).min(),np.concatenate([t0,t1]).max(), 500)
    dens_t0  = gaussian_kde(t0)(xaxis)
    dens_t1  = gaussian_kde(t1)(xaxis)
    dens_t00 = dens_t0/dens_t0.max()
    dens_t10 = dens_t1/dens_t1.max()

    ax.plot(xaxis,dens_t00, label='Stable', c=cond_cols[0])
    ax.plot(xaxis,dens_t10, label='Changing', c=cond_cols[1])
    ax.set(ylim=[0, 1.5], ylabel='Normalized Density', title=titles[p])
    ax.legend(title='Prevalence Condition')
    ax.fill_between(xaxis, dens_t00, color=cond_cols[0], alpha=0.35)
    ax.fill_between(xaxis, dens_t10, color=cond_cols[1], alpha=0.35)

    fig.savefig(f'figs/0_vs_1_density_{p}.png')
    plt.close()


pd.DataFrame(p_values, index=[0]).T.to_csv('stats/stim_comp_pvals.csv')







