"""Replay the archived Figure 1j-l MATLAB spatial analysis from selected objects."""
from pathlib import Path
import json,re,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def project_cells(table,crypt_id,villus_id,num_bins=50):
 if table.Id.duplicated().any():raise ValueError('Duplicate object IDs')
 anchors=table.set_index('Id').loc[[crypt_id,villus_id],['CentroidX','CentroidY']].to_numpy(float)
 if not np.isfinite(anchors).all():raise ValueError('Non-finite anchor centroids')
 base,tip=anchors;v=tip-base;length=float(v@v)
 if length<=0:raise ValueError('Crypt and villus anchors coincide')
 raw=(table[['CentroidX','CentroidY']].to_numpy(float)-base)@v/length
 # Four archived empty object rows have NaN centroids and recorded position zero.
 position=np.nan_to_num(np.clip(raw,0,1),nan=0.0)
 bins=np.minimum(np.searchsorted(np.linspace(0,1,num_bins+1),position,side='right'),num_bins)
 return position,bins

def summarize_bins(table,metrics,num_bins=50,percentiles=(10,25,50,75,90)):
 rows=[]
 for metric in metrics:
  for b in range(1,num_bins+1):
   values=table.loc[table.Bin.eq(b),metric].dropna().to_numpy(float);n=len(values)
   row=dict(metric=metric,Bin=b,BinCenter=(b-.5)/num_bins,N=n,Mean=float(values.mean()) if n else np.nan,SEM=float(values.std(ddof=1)/np.sqrt(n)) if n>1 else (0.0 if n else np.nan))
   quantiles=np.percentile(values,percentiles,method='hazen') if n else np.repeat(np.nan,len(percentiles));row.update({f'P{q:02d}':float(v) for q,v in zip(percentiles,quantiles)});row['Median']=row['P50'];rows.append(row)
 return pd.DataFrame(rows)

def svg_profile_checks(summary,svg_file,transform,tolerance):
 paths=[n for n in ET.parse(svg_file).getroot().iter() if n.tag.endswith('}path')]
 lines=[p for p in paths if 'stroke:none' not in p.get('style','') and p.get('d','').startswith('M137.75')]
 bands=[p for p in paths if 'stroke:none' in p.get('style','') and p.get('d','').startswith(('M137.75','M881.75'))]
 if len(lines)!=1 or len(bands)!=1:raise ValueError('Archived SVG profile paths are ambiguous')
 def points(node):return np.array([(float(a),float(b)) for a,b in re.findall(r'[ML]([\d.eE+-]+)\s+([\d.eE+-]+)',node.get('d'))])
 line=points(lines[0]);band=points(bands[0]);rows=[];scale=transform['scale'];offset=transform['offset'];s=summary.set_index('Bin')
 for x,y in line:
  b=int(round((x-137.75)/15.5))+1;row=s.loc[b];mean_error=abs(scale*row.Mean+offset-y);ys=band[np.isclose(band[:,0],x),1]
  if not len(ys):raise ValueError('Archived SVG error band missing a bin')
  bounds=np.sort([scale*(row.Mean-row.SEM)+offset,scale*(row.Mean+row.SEM)+offset]);sem_error=float(np.max(abs(bounds-[ys.min(),ys.max()])))
  rows.append(dict(check='archived_svg_mean_and_SEM',metric=row.metric,Bin=b,mean_error_svg_pixels=mean_error,sem_error_svg_pixels=sem_error,passed=bool(mean_error<tolerance and sem_error<tolerance)))
 if len(rows)!=50:raise ValueError('Expected 50 archived SVG bins')
 return rows

def run_spatial_cells(cfg,data,out,plots=True):
 params=json.loads((Path(__file__).resolve().parents[1]/cfg['parameter_file']).read_text());index=pd.read_csv(data/cfg['input']);tables=[];checks=[];per_sample=[]
 def compare(kind,expected,actual,**keys):
  ok=bool(np.isclose(expected,actual,rtol=1e-9,atol=1e-10,equal_nan=True));checks.append(dict(check=kind,reference=expected,recomputed=actual,passed=ok,**keys))
 for r in index.itertuples(index=False):
  t=pd.read_csv(data/r.object_file);archived=pd.read_csv(data/r.archived_projection)
  if archived.Id.duplicated().any() or set(t.Id)!=set(archived.Id):raise ValueError('Archived projection/object IDs do not match')
  position,bins=project_cells(t,int(r.crypt_cell_id),int(r.villus_cell_id),params['num_bins']);t['Proj01']=position;t['Bin']=bins;t['sample']=r.sample;t['source_object_file']=r.object_file
  ref=archived.set_index('Id').loc[t.Id]
  for oid,actual,expected,value,ref_value in zip(t.Id,position,ref.Proj01,t.MeanTauModulation,ref.MeanTauModulation):
   passed=bool(abs(actual-expected)<=params['projection_tolerance']);checks.append(dict(check='archived_cell_projection',sample=r.sample,Id=int(oid),reference=expected,recomputed=actual,passed=passed));compare('archived_cell_modulation',ref_value,value,sample=r.sample,Id=int(oid))
  tables.append(t);per_sample.append(summarize_bins(t,params['metrics'],params['num_bins'],params['percentiles']).assign(sample=r.sample))
 cells=pd.concat(tables,ignore_index=True);pooled=summarize_bins(cells,params['metrics'],params['num_bins'],params['percentiles']);samples=pd.concat(per_sample,ignore_index=True)
 ref=pd.read_csv(data/cfg['intensity_reference']);sample_lookup={v:k for k,v in cfg['source_sample_names'].items()};ref['sample']=ref.Sample.map(sample_lookup);actual=samples[samples.metric.eq('MeanInt')].set_index(['sample','Bin'])
 for r in ref.to_dict('records'):
  a=actual.loc[(r['sample'],r['Bin'])]
  for metric in ['BinCenter','N','Mean','SEM','P10','P25','P50','P75','P90','Median']:compare('archived_sample_intensity_summary',r[metric],a[metric],sample=r['sample'],Bin=r['Bin'],metric=metric)
 ref=pd.read_csv(data/cfg['phase_bins_reference']);ref['sample']=ref.Sample.map(sample_lookup);actual=cells.set_index(['sample','Id'])
 for r in ref.to_dict('records'):
  a=actual.loc[(r['sample'],r['Id'])];compare('archived_phase_bin',r['Bin'],a.Bin,sample=r['sample'],Id=r['Id']);compare('archived_phase_value',r['MeanTauPhase'],a.MeanTauPhase,sample=r['sample'],Id=r['Id'])
 for metric,path in params['archived_profiles'].items():checks.extend(svg_profile_checks(pooled[pooled.metric.eq(metric)],data/path,params['svg_axis_transforms'][metric],params['svg_tolerance_pixels']))
 verification=pd.DataFrame(checks);verification.to_csv(out/f'{cfg["panel"]}_spatial_verification.csv',index=False)
 if not verification.passed.all():raise ValueError('Spatial reproduction mismatch: '+str(verification[~verification.passed].iloc[0].to_dict()))
 cells.to_csv(out/f'{cfg["panel"]}_cells_with_projection.csv',index=False);pooled.to_csv(out/f'{cfg["panel"]}_pooled_bin_summary.csv',index=False);samples.to_csv(out/f'{cfg["panel"]}_per_sample_bin_summary.csv',index=False)
 if plots:
  labels={'MeanInt':('NADPH intensity','NAD(P)H photon/pixel'),'MeanTauPhase':('Tau phase','NAD(P)H τ(φ) ns'),'MeanTauModulation':('Tau modulation','NAD(P)H τ(m) ns')}
  with plt.rc_context({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':7,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.7,'axes.spines.top':False,'axes.spines.right':False}):
   fig,axes=plt.subplots(1,3,figsize=(7.2,2.7),layout='constrained')
   for ax,metric in zip(axes,params['metrics']):
    s=pooled[pooled.metric.eq(metric)];ax.plot(s.BinCenter,s.Mean,color='black',lw=1.1);ax.fill_between(s.BinCenter,s.Mean-s.SEM,s.Mean+s.SEM,color='.5',alpha=.2,linewidth=0);ax.set(xlim=(0,1),xlabel='Normalized crypt–villus position',ylabel=labels[metric][1]);ax.set_title(labels[metric][0],loc='left',fontsize=9.2,fontweight='bold')
   fig.suptitle('Pooled cell means ± SEM; five acquisitions',fontsize=8)
   for ext in ['png','svg']:fig.savefig(out/f'{cfg["panel"]}_spatial_profiles.{ext}',dpi=300)
   plt.close(fig)
 return dict(panel=cfg['panel'],status='computed',observations=len(cells),statistics_checks=0,spatial_checks=len(verification),spatial_projection_checks=int(verification.check.eq('archived_cell_projection').sum()),spatial_svg_checks=int(verification.check.eq('archived_svg_mean_and_SEM').sum()),mismatched=0,unmapped=0)
