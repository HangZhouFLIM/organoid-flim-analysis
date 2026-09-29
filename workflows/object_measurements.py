"""Read exact source object rows linked to the archived selected observations."""
from pathlib import Path
import numpy as np
import pandas as pd

def from_object_exports(obs, cfg, data_root):
 column=cfg.get('object_input_column')
 if not column:return obs,pd.DataFrame()
 source=pd.read_csv(data_root/cfg['input']);result=obs.copy();cache={};checks=[]
 for idx,row in result.iterrows():
  n=int(row.input_row)
  if n<0 or n>=len(source):raise ValueError('Missing object input row')
  origin=source.iloc[n];path=origin[column]
  if not isinstance(path,str) or not path or Path(path).is_absolute() or '..' in Path(path).parts:raise ValueError('Invalid source object path')
  if path not in cache:cache[path]=pd.read_csv(data_root/path)
  frame=cache[path];j=int(origin.source_object_row)
  if j<0 or j>=len(frame):raise ValueError('Object row outside the deposited export')
  obj=frame.iloc[j]
  if float(obj.Id)!=float(origin.source_object_id):raise ValueError('Object identifier does not match the deposited row')
  value=float(obj['Mean'+row.metric])
  passed=bool(np.isfinite(value) and np.isclose(value,row.value,rtol=0,atol=1e-8) and np.isclose(value,float(origin.value),rtol=0,atol=1e-8))
  checks.append(dict(input_file=path,object_row=j,object_id=int(float(obj.Id)),acquisition=origin.source_acquisition,metric=row.metric,reference=float(row.value),recomputed=value,passed=passed))
  result.at[idx,'value']=value
 checks=pd.DataFrame(checks)
 if len(checks) and not checks.passed.all():raise ValueError('Object-to-observation mismatch: '+str(checks[~checks.passed].iloc[0].to_dict()))
 return result,checks

def verify_object_masks(cfg, data_root):
 """Rebuild label/exclusion recipes and check archived cell areas and centroids."""
 if not cfg.get('object_mask_index'):return pd.DataFrame()
 from PIL import Image
 from scipy import ndimage
 index=pd.read_csv(data_root/cfg['object_mask_index']).fillna('');checks=[]
 for _,record in index.iterrows():
  labels=np.asarray(Image.open(data_root/record.source_labels))
  if record.recipe=='labels_where_exclusion_zero':
   exclusion=np.asarray(Image.open(data_root/record.exclusion_mask))
   if labels.shape!=exclusion.shape:raise ValueError('Cell label/exclusion shapes differ')
   labels=np.where(exclusion==0,labels,0)
  elif record.recipe!='identity':raise ValueError('Unknown object mask recipe')
  accepted=np.asarray(Image.open(data_root/record['mask']))
  if not np.array_equal(labels,accepted):raise ValueError('Rebuilt cell mask differs from its deposited mask')
  frame=pd.read_csv(data_root/record.object_file);frame=frame[pd.to_numeric(frame.Id,errors='coerce').notna()];ids=pd.to_numeric(frame.Id).to_numpy(int);counts=np.bincount(labels.ravel());counts=np.pad(counts,(0,max(0,int(ids.max())+1-len(counts))));centers=np.asarray(ndimage.center_of_mass(np.ones_like(labels),labels,ids));base=int(record.centroid_coordinate_base)
  for i,(_,row) in enumerate(frame.iterrows()):
   area=int(counts[ids[i]]);cy,cx=centers[i]+base;passed=bool(area==float(row.Area) and np.isclose(cx,float(row.CentroidX),rtol=0,atol=5.01e-5) and np.isclose(cy,float(row.CentroidY),rtol=0,atol=5.01e-5))
   checks.append(dict(acquisition=record['sample'],object_id=int(ids[i]),mask=record['mask'],area=area,centroid_x=cx,centroid_y=cy,passed=passed))
 checks=pd.DataFrame(checks)
 if not checks.passed.all():raise ValueError('Cell-mask geometry differs from the archived objects')
 return checks
