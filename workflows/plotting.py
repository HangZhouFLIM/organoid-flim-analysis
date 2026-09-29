"""Readable figure displays with the original statistical units and comparison results."""
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

NAMES = {'NADPH':'NAD(P)H','nadph':'NAD(P)H','nadph_mito':'NAD(P)H','tmrm':'TMRM','ATP Red 1':'ATP Red 1','nile_red_red':'Nile Red (red)','nile_red_green':'Nile Red (green)','nile_red_ratio':'Nile Red','lgr5_gfp':'Lgr5-GFP','muc2_yfp':'Muc2-YFP','cellrox':'CellROX','rostar':'RoStar'}
METRICS = {'Int':'Intensity','TauPhase':'Tau phase','TauModulation':'Tau modulation','TauNormalized':'Tau normalized','tissue_normalized_skeleton_density_um_per_um2':'Skeleton density','median_apparent_width_um':'Apparent width','green_red_ratio':'Green/red ratio'}

def label(channel, metric):
 channel=NAMES.get(channel,channel)
 if metric=='tissue_normalized_skeleton_density_um_per_um2':return 'Skeleton length / tissue area (µm/µm²)'
 if metric=='median_apparent_width_um':return 'Median apparent width (µm)'
 if metric=='green_red_ratio':return 'Green/red intensity ratio'
 if metric=='Int':return f'{channel} photon/pixel' if channel in ['NAD(P)H','TMRM'] else f'{channel} intensity (exported units)'
 if metric.startswith('Tau'):
  symbol={'TauPhase':'φ','TauModulation':'m','TauNormalized':'N'}[metric]
  return f'{channel} τ({symbol}) ns'
 return metric

def stars(p):
 return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'ns'

def draw_points(ax,g,group_column,order,calibration=False):
 xmap={k:float(k) if calibration else i*.5 for i,k in enumerate(order)}
 bounds=[]
 for i,k in enumerate(order):
  vals=g.loc[g[group_column].eq(k),'value'].to_numpy(float);n=len(vals);mean=float(np.mean(vals));sem=float(stats.sem(vals)) if n>1 else 0.
  err=sem if calibration else sem*stats.t.ppf(.975,n-1) if n>1 else 0.
  jitter=np.random.default_rng(100+i).uniform(-.06,.06,n)*(5 if calibration else 1)
  ax.scatter(xmap[k]+jitter,vals,s=9,alpha=.6,edgecolors='none',color=plt.get_cmap('tab10')(i%10),zorder=2)
  ax.errorbar(xmap[k],mean,yerr=err,fmt='o',ms=3,color='black',capsize=2,zorder=3)
  bounds.extend([float(np.min(vals)),float(np.max(vals)),mean-err,mean+err])
 ax.set_xticks(list(xmap.values()),[str(k) for k in order],rotation=35 if not calibration else 0,ha='right' if not calibration else 'center')
 if not calibration:ax.set_xlim(min(xmap.values())-.35,max(xmap.values())+.35)
 return xmap,min(bounds),max(bounds)

def bracket(ax,x1,x2,y,height,p):
 ax.plot([x1,x1,x2,x2],[y,y+height,y+height,y],color='black',lw=.55)
 ax.text((x1+x2)/2,y+height,stars(float(p)),ha='center',va='bottom',fontsize=7)

def plot_endpoints(obs,cfg,out,actual):
 plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':7,'axes.linewidth':.7,'svg.fonttype':'none'})
 paired=cfg['mode'] in ['paired','nile'];calibration=cfg['mode']=='calibration'
 fields=['condition','channel','metric'] if paired else ['region','channel','metric']
 allowed=cfg.get('plot_metrics',cfg.get('metrics',['Int','TauPhase','TauModulation']))
 facets=list(obs.groupby(fields,sort=False));display=[f for f in facets if f[0][2] in allowed]
 if not display:display=facets
 for page,start in enumerate(range(0,len(display),6),1):
  chunk=display[start:start+6];fig,axes=plt.subplots(3,2,figsize=(6.3,8.7),squeeze=False)
  for ax,(key,g) in zip(axes.flat,chunk):
   group_column='region' if paired else 'condition'
   order=[r for r in cfg.get('pair_regions',['crypt','villus']) if r in set(g.region)] if paired else list(dict.fromkeys(g.condition))
   if not order:order=list(dict.fromkeys(g[group_column]))
   if calibration:order=sorted(order,key=float)
   xmap,lo,hi=draw_points(ax,g,group_column,order,calibration)
   span=max(hi-lo,abs(hi)*.05,1e-6)
   tests=actual
   for field,value in [('channel',key[1]),('metric',key[2]),(fields[0],key[0])]:
    if field in tests:tests=tests[tests[field].eq(value)]
   if paired:
    wide=g.pivot(index='sample',columns='region',values='value').reindex(columns=order).dropna()
    ax.plot([xmap[k] for k in order],wide.to_numpy().T,color='.7',lw=.5,zorder=1)
    tests=tests[tests.kind.eq('paired')] if 'kind' in tests else tests.iloc[:0]
    if len(tests)==1 and len(order)==2:
     t=tests.iloc[0];bracket(ax,xmap[order[0]],xmap[order[1]],hi+.10*span,.025*span,t.get('p_adjusted',t.p_value));hi+=.25*span
   elif not calibration and 'kind' in tests:
    tests=tests[tests.kind.eq('posthoc')];level=0
    for _,t in tests.iterrows():
     if t.reference not in xmap or t.condition not in xmap:continue
     bracket(ax,xmap[t.reference],xmap[t.condition],hi+(.10+.18*level)*span,.025*span,t.get('p_adjusted',t.p_value));level+=1
    hi+=.18*span*level
   ax.set_ylim(lo-.10*span,hi+.12*span)
   channel=cfg.get('display_channel',key[1]);title=METRICS.get(key[2],key[2]);prefix=NAMES.get(channel,channel)
   if key[0]!='all':prefix=str(key[0])+' / '+prefix
   ax.set_title(prefix+'\n'+title,loc='left',fontsize=8)
   ax.set_ylabel(label(channel,key[2]));ax.set_box_aspect(1.25);ax.spines[['top','right']].set_visible(False)
   if calibration:ax.set_xlabel('ATP concentration (mM)')
  for ax in list(axes.flat)[len(chunk):]:ax.set_visible(False)
  fig.suptitle(cfg['panel'],x=.05,ha='left',fontsize=9.2,fontweight='bold');fig.tight_layout(rect=(0,0,1,.975));stem=f'{cfg["panel"]}_reproduced_{page:02d}'
  for ext in ['png','svg']:fig.savefig(out/f'{stem}.{ext}',dpi=300,bbox_inches='tight')
  plt.close(fig)
