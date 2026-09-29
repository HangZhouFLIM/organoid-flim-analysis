"""Additional figure-specific displays and derived measurements."""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

def extra_outputs(obs,cfg,out,plots=True):
 panel=cfg['panel']
 if cfg.get('photon_group_map'):
  grouped=obs.copy();grouped['photon_group']=grouped.condition.map(cfg['photon_group_map'])
  if grouped.photon_group.isna().any():raise ValueError('An acquisition setting has no photon-group assignment')
  rows=[]
  for (group,metric),g in grouped.groupby(['photon_group','metric'],sort=False):
   v=g.value.to_numpy(float);rows.append(dict(photon_group=group,metric=metric,n_acquisition_measurements=len(v),n_settings=g.condition.nunique(),mean=float(np.mean(v)),sd=float(np.std(v,ddof=1)),sem=float(stats.sem(v)),interpretation='Descriptive pooling of acquisition-setting measurements; not an independent biological replicate count'))
  pd.DataFrame(rows).to_csv(out/f'{panel}_photon_group_summary.csv',index=False)
 if cfg['mode']=='nile':
  d=obs[obs.metric.eq('Int')&obs.channel.isin(['nile_red_green','nile_red_red'])]
  d=d.groupby(['sample','channel']).value.mean().unstack();delta=d.nile_red_red-d.nile_red_green
  r=stats.wilcoxon(delta,alternative='two-sided',method='auto')
  pd.DataFrame([dict(n_pairs=len(d),test='Wilcoxon signed-rank',p_value=r.pvalue,p_adjusted=r.pvalue,adjustment='none; separate singleton comparison',green_mean=d.nile_red_green.mean(),red_mean=d.nile_red_red.mean())]).to_csv(out/f'{panel}_emission_statistics.csv',index=False)
  d.to_csv(out/f'{panel}_equal_region_emission_means.csv')
  if plots:
   fig,ax=plt.subplots(figsize=(2.4,3));ax.plot([0,.5],d[['nile_red_green','nile_red_red']].to_numpy().T,lw=.5,color='.6');ax.scatter(np.zeros(len(d)),d.nile_red_green,s=12);ax.scatter(np.full(len(d),.5),d.nile_red_red,s=12);ax.set_xticks([0,.5],['Green','Red']);ax.set_ylabel('Nile Red intensity (exported units)');fig.tight_layout()
   for ext in ['png','svg']:fig.savefig(out/f'{panel}_emission_comparison.{ext}',dpi=300)
   plt.close(fig)
 if cfg['mode']=='calibration' and plots:
  d=obs.pivot(index=['condition','sample'],columns='metric',values='value').reset_index();d['concentration_mM']=d.condition.astype(float)
  fig,ax=plt.subplots(figsize=(4.2,3.5));theta=np.linspace(0,np.pi,300);ax.plot(.5+.5*np.cos(theta),.5*np.sin(theta),color='.5',lw=.7)
  centers=[]
  for concentration,g in d.groupby('concentration_mM'):
   ax.scatter(g.G,g.S,s=6,alpha=.4);ax.errorbar(g.G.mean(),g.S.mean(),xerr=g.G.sem(),yerr=g.S.sem(),fmt='o',ms=3,capsize=2);ax.annotate(f'{concentration:g} mM',(g.G.mean(),g.S.mean()),fontsize=6,xytext=(3,4),textcoords='offset points');centers.append((g.G.mean(),g.S.mean()))
  centers=np.asarray(centers);ax.plot(centers[:,0],centers[:,1],lw=.7,color='.4');ax.set(xlim=(0,1),ylim=(0,.5),xlabel='G',ylabel='S');ax.set_aspect('equal');fig.tight_layout()
  for ext in ['png','svg']:fig.savefig(out/f'{panel}_phasor_calibration.{ext}',dpi=300)
  plt.close(fig)

def verify_additional(cfg,data,out,plots=False):
 """Verify derived displays against their separately deposited reference tables."""
 panel=cfg['panel'];checks=[]
 if panel=='FigS07b':
  frame=pd.read_csv(data/cfg['input']);means=frame.groupby(['cell_type_label','pathway']).sample_mean_score.mean().unstack();z=(means-means.mean())/means.std(ddof=0)
  expected=pd.read_csv(data/cfg['summary_reference'],sep='\t').set_index('cell_type')
  for cell in expected.index:
   for metric in means.columns:
    for suffix,actual in [('',means),('_column_z',z)]:
     ref=float(expected.loc[cell,metric+suffix]);value=float(actual.loc[cell,metric]);checks.append(dict(item=f'{cell}/{metric}{suffix}',reference=ref,recomputed=value,passed=bool(np.isclose(ref,value,rtol=1e-8,atol=1e-10))))
 if panel=='FigS07c':
  frame=pd.read_csv(data/cfg['input']);wide=frame.pivot(index='normal_sample_label',columns='cell_type_label',values='complex_v_structural_nuclear_mean');expected=pd.read_csv(data/cfg['summary_reference'],sep='\t');rows=[]
  for _,r in expected.iterrows():
   other=r.comparison.removeprefix('Enterocyte minus ');value=float(wide.loc[r.normal_sample_label,'Enterocyte']-wide.loc[r.normal_sample_label,other]);rows.append(dict(normal_sample_label=r.normal_sample_label,comparison=r.comparison,complex_v_difference=value));checks.append(dict(item=r.normal_sample_label+'/'+r.comparison,reference=r.complex_v_difference,recomputed=value,passed=bool(np.isclose(value,r.complex_v_difference,rtol=1e-8,atol=1e-10))))
  pd.DataFrame(rows).to_csv(out/f'{panel}_enterocyte_differences.csv',index=False)
  if plots:
   fig,ax=plt.subplots(figsize=(4.2,3))
   for i,(comparison,g) in enumerate(pd.DataFrame(rows).groupby('comparison',sort=False)):
    vals=g.complex_v_difference;ax.scatter(np.full(len(vals),i*.5),vals,s=12);ax.errorbar(i*.5,vals.mean(),yerr=stats.sem(vals)*stats.t.ppf(.975,len(vals)-1),fmt='o',color='black',ms=3,capsize=2)
   ax.axhline(0,lw=.7,color='.6');labels=pd.DataFrame(rows).comparison.drop_duplicates();ax.set_xticks(np.arange(len(labels))*.5,[v.replace('Enterocyte minus ','') for v in labels],rotation=25,ha='right');ax.set_ylabel('Complex V score: enterocyte minus comparator');fig.tight_layout()
   for ext in ['png','svg']:fig.savefig(out/f'{panel}_enterocyte_differences.{ext}',dpi=300)
   plt.close(fig)
 if panel=='Fig06g':
  actual=pd.read_csv(out/f'{panel}_emission_statistics.csv').iloc[0];reference=pd.read_csv(data/cfg['emission_reference']).iloc[0]
  for field in ['p_adjusted','green_mean','red_mean','n_pairs']:
   checks.append(dict(item='emission/'+field,reference=reference[field],recomputed=actual[field],passed=bool(np.isclose(reference[field],actual[field],rtol=1e-8,atol=1e-10))))
 if checks:
  checks=pd.DataFrame(checks);checks.to_csv(out/f'{panel}_additional_verification.csv',index=False)
  if not checks.passed.all():raise ValueError('Derived display verification failed: '+str(checks[~checks.passed].iloc[0].to_dict()))
 return len(checks)
