"""Portable article reproduction from the deposited, selected measurement exports."""
from __future__ import annotations
import argparse,hashlib,json,math,re,sys,traceback
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from workflows.pixel_measurements import from_pixels, from_measurement_masks
from workflows.object_measurements import from_object_exports, verify_object_masks
from workflows.specialized_outputs import extra_outputs, verify_additional
from workflows.plotting import plot_endpoints
from workflows.raw_inputs import validate_raw_inputs
from workflows.spatial_cells import run_spatial_cells

METRIC_ALIASES={'photon_pixel':'Int','tauphase':'TauPhase','taumodulation':'TauModulation','taunormalized':'TauNormalized'}

def holm(p):
 p=np.asarray(p,float);order=np.argsort(p);result=np.empty_like(p)
 result[order]=np.minimum(1,np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
 return result

def summary(v):
 v=np.asarray(v,float);v=v[np.isfinite(v)];n=len(v);m=float(np.mean(v)) if n else np.nan
 sd=float(np.std(v,ddof=1)) if n>1 else np.nan;sem=sd/np.sqrt(n) if n else np.nan
 ci=float(stats.t.ppf(.975,n-1)*sem) if n>1 else np.nan
 return dict(n=n,mean=m,sd=sd,sem=sem,ci95_low=m-ci,ci95_high=m+ci)

def normal(v):return len(v)>=3 and bool(stats.shapiro(v).pvalue>=.05)

def welch_anova(groups):
 n=np.array([len(x) for x in groups]);means=np.array([np.mean(x) for x in groups]);var=np.array([np.var(x,ddof=1) for x in groups]);k=len(groups)
 w=n/var;total=w.sum();center=np.sum(w*means)/total
 c=np.sum((1-w/total)**2/(n-1));F=np.sum(w*(means-center)**2)/(k-1)/(1+2*(k-2)*c/(k*k-1));df2=(k*k-1)/(3*c)
 return float(F),float(stats.f.sf(F,k-1,df2)),float(k-1),float(df2)

def endpoint(groups,ref):
 keys=list(groups);arrays=list(groups.values());allnormal=all(normal(v) for v in arrays);out=[]
 for k,v in groups.items():
  r=stats.shapiro(v) if len(v)>=3 else None
  out.append(dict(kind='normality',condition=k,test='Shapiro-Wilk',statistic=r.statistic if r else np.nan,p_value=r.pvalue if r else np.nan,n=len(v)))
 if allnormal:
  r=stats.levene(*arrays,center='median');out.append(dict(kind='variance',condition='all',test='Brown-Forsythe',statistic=r.statistic,p_value=r.pvalue))
  F,p,df1,df2=welch_anova(arrays);out.append(dict(kind='omnibus',condition='all',test='Welch ANOVA',statistic=F,p_value=p,df1=df1,df2=df2))
  a=groups[ref]
  for key,b in groups.items():
   if key==ref:continue
   v1,v2=np.var(a,ddof=1)/len(a),np.var(b,ddof=1)/len(b);se=v1+v2;delta=float(np.mean(b)-np.mean(a));t=abs(delta)/np.sqrt(se);df=se**2/(v1**2/(len(a)-1)+v2**2/(len(b)-1));p=float(stats.studentized_range.sf(t*np.sqrt(2),len(groups),df))
   out.append(dict(kind='posthoc',condition=key,reference=ref,test='Games-Howell-type',statistic=t,p_value=p,p_adjusted=p,df2=df,effect=delta,n_reference=len(a),n_test=len(b),adjustment='studentized range across all groups'))
 else:
  r=stats.kruskal(*arrays);out.append(dict(kind='omnibus',condition='all',test='Kruskal-Wallis',statistic=r.statistic,p_value=r.pvalue))
  pooled=np.concatenate(arrays);ranks=stats.rankdata(pooled);N=len(pooled);_,counts=np.unique(pooled,return_counts=True);tie=1-np.sum(counts**3-counts)/(N**3-N)
  means={};offset=0
  for k,a in groups.items():means[k]=np.mean(ranks[offset:offset+len(a)]);offset+=len(a)
  post=[]
  for k,b in groups.items():
   if k==ref:continue
   a=groups[ref];z=(means[k]-means[ref])/np.sqrt(N*(N+1)/12*tie*(1/len(a)+1/len(b)));p=float(2*stats.norm.sf(abs(z)))
   post.append(dict(kind='posthoc',condition=k,reference=ref,test='Dunn',statistic=z,p_value=p,effect=float(np.mean(b)-np.mean(a)),n_reference=len(a),n_test=len(b),adjustment='Holm across treatment-versus-control comparisons'))
  for row,p in zip(post,holm([r['p_value'] for r in post])):row['p_adjusted']=p
  out.extend(post)
 return out

def paired_table(obs,cfg,nile=False):
 rows=[];regions=cfg.get('pair_regions',['crypt','villus']);a,b=regions
 for (condition,channel,metric),g in obs.groupby(['condition','channel','metric'],sort=False):
  if 'metrics' in cfg and metric not in cfg['metrics']:continue
  pairs=g.pivot(index='sample',columns='region',values='value')
  if a not in pairs or b not in pairs:continue
  pairs=pairs[[a,b]].dropna();delta=(pairs[a]-pairs[b]).to_numpy() if nile else (pairs[b]-pairs[a]).to_numpy()
  is_normal=cfg.get('paired_test')=='ttest' or normal(delta)
  if is_normal:r=stats.ttest_1samp(delta,0);test='Paired t-test'
  elif np.allclose(delta,0):r=type('R',(),{'statistic':0,'pvalue':1});test='Wilcoxon signed-rank'
  else:r=stats.wilcoxon(delta,alternative='two-sided',method='auto');test='Wilcoxon signed-rank'
  family=(str(condition)+'|primary_2' if metric=='Int' and channel in ['nile_red_red','nile_red_green'] else str(condition)+'|secondary_11') if nile else str(condition)+'|'+str(channel)
  rows.append(dict(kind='paired',condition=condition,channel=channel,region=a+'-'+b,metric=metric,test=test,statistic=float(r.statistic),p_value=float(r.pvalue),family=family,effect=float(np.mean(delta)),n_pairs=len(delta),difference=(a+' - '+b) if nile else (b+' - '+a),**{'difference_'+k:v for k,v in summary(delta).items()}))
 if cfg.get('paired_family')=='across_conditions':
  for r in rows:r['family']=r['channel']+'|'+r['metric']
 for family in set(r['family'] for r in rows):
  rr=[r for r in rows if r['family']==family]
  for r,p in zip(rr,holm([x['p_value'] for x in rr])):r['p_adjusted']=p
 return rows

def analyse(obs,cfg):
 mode=cfg['mode'];rows=[]
 if mode in ['paired','nile']:return pd.DataFrame(paired_table(obs,cfg,nile=mode=='nile'))
 if mode=='calibration':
  for m,g in obs.groupby('metric',sort=False):
   for level,d in [('image',g),('concentration',g.groupby('condition',as_index=False).value.mean())]:
    r=stats.spearmanr(d.condition.astype(float),d.value)
    rows.append(dict(kind='descriptive_correlation',metric=m,analysis_unit=level,test='Spearman',statistic=r.statistic,p_value=r.pvalue,n=len(d)))
  return pd.DataFrame(rows)
 for (region,channel,metric),g in obs.groupby(['region','channel','metric'],sort=False):
  if cfg.get('metrics') and metric not in cfg['metrics']:continue
  groups={str(k):v.value.to_numpy(float) for k,v in g.groupby('condition',sort=False)}
  if len(groups)<2:continue
  ref=cfg.get('reference','control');ref=next(iter(groups)) if ref=='first' else ref
  if ref not in groups:raise ValueError(f'{cfg["panel"]}: missing reference {ref}; groups={list(groups)}')
  if mode=='endpoint':out=endpoint(groups,ref)
  else:
   out=[]
   comparisons=[(ref,key) for key in groups if key!=ref]
   if cfg.get('all_pairs'):
    keys=list(groups);comparisons=[(keys[i],keys[j]) for i in range(len(keys)) for j in range(i+1,len(keys))]
   comparisons=cfg.get('comparisons',comparisons)
   for ref,key in comparisons:
    b=groups[key]
    if key==ref:continue
    a=groups[ref];method=cfg.get('test','welch');use_mw=method=='mannwhitney' or (method=='normality_two_group' and not (normal(a) and normal(b)))
    r=stats.mannwhitneyu(a,b,alternative='two-sided') if use_mw else stats.ttest_ind(a,b,equal_var=False,alternative='two-sided')
    out.append(dict(kind='posthoc',condition=key,reference=ref,test='Mann-Whitney U' if use_mw else 'Welch t-test',statistic=r.statistic,p_value=r.pvalue,effect=float(np.mean(b)-np.mean(a)),n_reference=len(a),n_test=len(b)))
   corrected=holm([r['p_value'] for r in out]) if cfg.get('adjustment','holm')=='holm' else [r['p_value'] for r in out]
   for r,p in zip(out,corrected):r['p_adjusted']=p
  rows.extend(dict(r,region=region,channel=channel,metric=metric) for r in out)
 if cfg.get('also_paired'):rows.extend(paired_table(obs,{**cfg,'metrics':cfg.get('metrics',['Int','TauPhase','TauModulation'])}))
 return pd.DataFrame(rows)

def expected_checks(actual,reference,cfg):
 checks=[]
 if actual.empty or reference is None:return pd.DataFrame()
 for i,r in reference.iterrows():
  test=str(r.get('test',r.get('statistical_test',r.get('primary_test','')))).lower();level=str(r.get('analysis_level',r.get('row_type',''))).lower()
  if any(k in test for k in ['shapiro','brown','anova','kruskal']) or level in ['normality','omnibus','omnibus_variance']:continue
  if test=='group summary':continue
  pcol=next((k for k in ['p_adjusted','p_adjusted_holm','p_raw','p_value','welch_p'] if k in r and pd.notna(r[k])),None)
  if not pcol:continue
  candidates=actual.copy()
  for target,fields in [('metric',['metric','metric_key']),('channel',['channel','marker_key']),('region',['region']),('condition',['condition_key','condition','test_condition_key','drug']),('reference',['reference_condition']),('analysis_unit',['analysis_unit'])]:
   field=next((f for f in fields if f in r and pd.notna(r[f]) and str(r[f])!=''),None)
   if field is None or target not in candidates:continue
   val=str(r[field]);val=METRIC_ALIASES.get(val,{'intensity':'Int','tau_phase':'TauPhase','tau_modulation':'TauModulation','phasor_g_intensity_weighted':'G','phasor_s_intensity_weighted':'S'}.get(val,val))
   if target=='analysis_unit' and val=='concentration mean':val='concentration'
   if target=='metric':val=val.split(' ')[-1];val=val.removesuffix('_mean')
   if target=='condition' and val in ['all']:continue
   matches=candidates[target].astype(str).str.lower().eq(val.lower())
   if matches.any():candidates=candidates[matches]
   elif target!='region':candidates=candidates.iloc[0:0]
  metric_text=str(r.get('metric',''))
  if metric_text.startswith(('TMRM ', 'NADPH ', 'NAD(P)H ')) and 'channel' in candidates:
   marker=metric_text.split(' ')[0].replace('NAD(P)H','NADPH')
   candidates=candidates[candidates.channel.astype(str).str.lower().eq(marker.lower())]
  kind='paired' if any(k in test for k in ['paired','one-sample','signed-rank']) or cfg['mode'] in ['paired','nile'] else 'posthoc'
  if 'spearman' in test:kind='descriptive_correlation'
  candidates=candidates[candidates.kind.eq(kind)]
  target='p_adjusted' if pcol in ['p_adjusted','p_adjusted_holm'] else 'p_value'
  expected=float(r[pcol]);vals=candidates[target].dropna().to_numpy(float) if target in candidates else np.array([])
  matches=np.isclose(vals,expected,rtol=1e-6,atol=1e-12)
  checks.append(dict(reference_row=i,metric=r.get('metric',r.get('metric_key','')),test=r.get('test',''),reference_p=expected,computed_p=float(vals[0]) if len(vals)==1 else np.nan,candidate_count=len(vals),status='pass' if len(vals)==1 and matches[0] else 'ambiguous' if len(vals)>1 else 'mismatch' if len(vals) else 'unmapped',p_column=pcol))
 return pd.DataFrame(checks)

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def verify_data(root,manifest):
 errors=[]
 for r in manifest:
  p=(root/r['file']).resolve()
  if not p.is_relative_to(root.resolve()):errors.append(r['file']+': unsafe path');continue
  if not p.exists():errors.append(r['file']+': missing')
  elif sha(p)!=r['sha256']:errors.append(r['file']+': checksum mismatch')
 return errors

def special(cfg,data,out,plots=True):
 frame=pd.read_csv(data/cfg['input']);mode=cfg['mode'];outputs=[]
 if mode=='spatial_cells':return run_spatial_cells(cfg,data,out,plots=plots)
 if mode=='valid_percentage':
  frame['percentage']=100*frame['count']/frame['total_count'];frame.to_csv(out/f'{cfg["panel"]}_recalculated_percentages.csv',index=False)
  fig,ax=plt.subplots(figsize=(6.3,3));labels=frame.get('condition_label',frame.get('photon_label',frame.index.astype(str))).astype(str)+' '+frame.metric_key.astype(str)+' '+frame.roi_class.astype(str);ax.bar(np.arange(len(frame)),frame.percentage);ax.set_xticks(np.arange(len(frame)),labels,rotation=90,fontsize=6);ax.set_ylabel('ROI (%)');fig.tight_layout()
 elif mode=='spectra':
  fig,ax=plt.subplots(figsize=(4.2,3))
  for label,g in frame.groupby('detection_window',sort=False):ax.plot(g.excitation_nm,g.intensity,'o-',ms=3,label=label.replace('_',' '))
  ax.set_xlabel('Excitation wavelength (nm)');ax.set_ylabel('Exported intensity');ax.legend(frameon=False,fontsize=6);fig.tight_layout()
  frame.to_csv(out/f'{cfg["panel"]}_scan_values.csv',index=False)
 elif mode=='profile':
  fig,axes=plt.subplots(1,3,figsize=(8,2.8))
  for ax,m in zip(axes,['Int','TauPhase','TauModulation']):
   for label,g in frame.groupby('Group',sort=False):
    ax.plot(g.BinCenter,g['Mean_'+m],label=label);ax.fill_between(g.BinCenter,g['Mean_'+m]-1.96*g['SEM_'+m],g['Mean_'+m]+1.96*g['SEM_'+m],alpha=.16)
   ax.set_xlabel('Normalized crypt-villus position');ax.set_title({'Int':'Intensity','TauPhase':'Tau phase','TauModulation':'Tau modulation'}[m]);ax.legend(fontsize=5,frameon=False)
  fig.tight_layout();frame.to_csv(out/f'{cfg["panel"]}_profile_summary.csv',index=False)
 elif mode=='scrna':
  if 'pathway' in frame:
   means=frame.groupby(['cell_type_label','pathway'],sort=False).sample_mean_score.mean().unstack();z=(means-means.mean())/means.std(ddof=0);means.to_csv(out/f'{cfg["panel"]}_equal_sample_means.csv');z.to_csv(out/f'{cfg["panel"]}_column_z_scores.csv')
   fig,ax=plt.subplots(figsize=(4.2,3));im=ax.imshow(z,aspect='auto',cmap='RdBu_r');ax.set_xticks(range(len(z.columns)),z.columns);ax.set_yticks(range(len(z.index)),z.index);fig.colorbar(im,ax=ax,label='Column z score');fig.tight_layout()
  else:
   col='complex_v_structural_nuclear_mean';fig,ax=plt.subplots(figsize=(4.2,3));types=list(frame.cell_type_label.unique())
   for sample,g in frame.groupby('normal_sample_label'):
    vals=g.set_index('cell_type_label')[col].reindex(types);ax.plot(range(len(types)),vals,'o-',ms=3,lw=.6,alpha=.6)
   ax.set_xticks(range(len(types)),types,rotation=25,ha='right');ax.set_ylabel('Complex V module score');fig.tight_layout();frame.groupby('cell_type_label')[col].agg(['count','mean','std']).to_csv(out/f'{cfg["panel"]}_sample_summary.csv')
 elif mode=='timelapse':
  rows=[]
  for (channel,ti),g in frame.groupby(['channel','time_index'],sort=False):
   for m in ['Int_mean','IntRel_mean','TauPhase_mean','TauModulation_mean']:
    if m in g:rows.append(dict(channel=channel,time_index=ti,elapsed_min=float(g.elapsed_min.mean()),metric=m,**summary(g[m])))
  summary_frame=pd.DataFrame(rows);summary_frame.to_csv(out/f'{cfg["panel"]}_timepoint_summary.csv',index=False)
  fig,axes=plt.subplots(2,2,figsize=(6.3,5),squeeze=False)
  for ax,m in zip(axes.flat,['Int_mean','IntRel_mean','TauPhase_mean','TauModulation_mean']):
   for channel,g in summary_frame[summary_frame.metric.eq(m)].groupby('channel'):
    ax.plot(g.elapsed_min,g['mean'],label=channel);ax.fill_between(g.elapsed_min,g.ci95_low,g.ci95_high,alpha=.16)
   ax.set_xlabel('Elapsed time (min)');ax.set_title({'Int_mean':'Intensity','IntRel_mean':'Intensity F/F0','TauPhase_mean':'Tau phase','TauModulation_mean':'Tau modulation'}[m]);ax.legend(frameon=False,fontsize=6)
  fig.tight_layout()
 else:raise ValueError('Unsupported workflow '+mode)
 if plots:
  for ext in ['png','svg']:fig.savefig(out/f'{cfg["panel"]}_reproduced.{ext}',dpi=300)
 plt.close(fig)
 extra_checks=verify_additional(cfg,data,out,plots=plots)
 return dict(panel=cfg["panel"],status='computed',statistics_checks=0,additional_checks=extra_checks)

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--data-root',type=Path,help='Extracted Zenodo data folder containing metadata/ and Fig*/')
 actions=parser.add_mutually_exclusive_group(required=True)
 actions.add_argument('--panel',help='Panel/workflow identifier or companion alias; see --list-panels')
 actions.add_argument('--all',action='store_true')
 actions.add_argument('--verify-only',action='store_true')
 actions.add_argument('--list-panels',action='store_true',help='List configured workflows without requiring the data download')
 parser.add_argument('--output',type=Path,default=Path('results'))
 parser.add_argument('--skip-hashes',action='store_true',help='Skip bulk experimental-file hashes for repeated local debugging; release indexes are still checked')
 parser.add_argument('--no-plots',action='store_true')
 args=parser.parse_args();here=Path(__file__).resolve().parent
 configs=[json.loads(f.read_text()) for f in sorted((here/'configs').glob('*.json'))]
 if args.list_panels:
  for cfg in configs:print(cfg['panel']+'\t'+cfg['mode']+'\t'+', '.join(cfg.get('aliases',[])))
  return
 if args.data_root is None:parser.error('--data-root is required unless --list-panels is used')
 if args.panel:
  requested=args.panel
  configs=[c for c in configs if requested==c['panel'] or requested in c.get('aliases',[]) or c['panel'].startswith(requested+'_')]
  if not configs:parser.error('Unknown panel; choose an identifier from --list-panels')
 manifest=json.loads((here/'manifests/data_files.json').read_text())
 release=json.loads((here/'manifests/data_release.json').read_text())
 errors=verify_data(args.data_root,release['metadata_files'])
 if errors:raise RuntimeError('Zenodo release metadata differs from this code package:\n'+'\n'.join(errors[:30]))
 if not args.skip_hashes:
  errors=verify_data(args.data_root,manifest)
  if errors:raise RuntimeError('\n'.join(errors[:30]))
 errors=validate_raw_inputs(args.data_root,manifest)
 if errors:raise RuntimeError('\n'.join(errors[:30]))
 if args.verify_only:print(f'Validated {len(manifest)} deposited files, matching release metadata and raw acquisition crosswalks');return
 args.output.mkdir(parents=True,exist_ok=True)
 reports=[]
 for cfg in configs:
  out=args.output/cfg["panel"];out.mkdir(exist_ok=True)
  try:
   if not cfg.get('observations'):reports.append(special(cfg,args.data_root,out,plots=not args.no_plots));continue
   obs=pd.read_csv(args.data_root/cfg['observations'],dtype={'condition':str,'sample':str,'region':str,'channel':str});obs=obs[np.isfinite(obs.value)]
   object_mask_checks=verify_object_masks(cfg,args.data_root)
   if len(object_mask_checks):object_mask_checks.to_csv(out/f'{cfg["panel"]}_object_mask_verification.csv',index=False)
   obs,object_checks=from_object_exports(obs,cfg,args.data_root)
   if len(object_checks):object_checks.to_csv(out/f'{cfg["panel"]}_object_verification.csv',index=False)
   obs,pixel_checks=from_pixels(obs,cfg,args.data_root)
   obs,mask_checks=from_measurement_masks(obs,cfg,args.data_root)
   if len(mask_checks):mask_checks.to_csv(out/f'{cfg["panel"]}_mask_verification.csv',index=False)
   if len(pixel_checks):pixel_checks.to_csv(out/f'{cfg["panel"]}_pixel_verification.csv',index=False)
   if cfg['mode']=='nile':
    ratios=[]
    sub=obs[obs.metric.eq('Int')&obs.channel.isin(['nile_red_green','nile_red_red'])]
    for (condition,sample,region),g in sub.groupby(['condition','sample','region']):
     vals=g.set_index('channel').value;ratios.append(dict(condition=condition,sample=sample,region=region,channel='nile_red_ratio',metric='green_red_ratio',value=vals['nile_red_green']/vals['nile_red_red'],input_row=-1))
    obs=pd.concat([obs,pd.DataFrame(ratios)],ignore_index=True)
   sums=[]
   for key,g in obs.groupby(['condition','region','channel','metric'],sort=False):sums.append(dict(zip(['condition','region','channel','metric'],key),**summary(g.value)))
   pd.DataFrame(sums).to_csv(out/f'{cfg["panel"]}_summary.csv',index=False)
   actual=analyse(obs,cfg);actual.to_csv(out/f'{cfg["panel"]}_statistics.csv',index=False)
   reference=pd.read_csv(args.data_root/cfg['reference_statistics']) if cfg.get('reference_statistics') else None
   checks=expected_checks(actual,reference,cfg);checks.to_csv(out/f'{cfg["panel"]}_verification.csv',index=False)
   extra_outputs(obs,cfg,out,plots=not args.no_plots)
   additional_checks=verify_additional(cfg,args.data_root,out,plots=not args.no_plots)
   if not args.no_plots:plot_endpoints(obs,cfg,out,actual)
   reports.append(dict(panel=cfg["panel"],status='computed',observations=len(obs),pixel_checks=len(pixel_checks),object_checks=len(object_checks),object_mask_checks=len(object_mask_checks),mask_checks=len(mask_checks),additional_checks=additional_checks,statistics_checks=len(checks),passed=int(checks.status.eq('pass').sum()) if len(checks) else 0,mismatched=int(checks.status.isin(['mismatch','ambiguous']).sum()) if len(checks) else 0,unmapped=int(checks.status.eq('unmapped').sum()) if len(checks) else 0))
   print(cfg["panel"],reports[-1],flush=True)
  except Exception as e:
   reports.append(dict(panel=cfg["panel"],status='error',error=str(e)));traceback.print_exc()
 pd.DataFrame(reports).to_csv(args.output/'reproduction_report.csv',index=False)
 (args.output/'reproduction_report.json').write_text(json.dumps(reports,indent=2))
 if any(r['status']=='error' or r.get('mismatched',0) or r.get('unmapped',0) for r in reports):sys.exit(2)
if __name__=='__main__':main()
