"""Recompute image-level measurements from selected pixel exports when available."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import ndimage

def from_pixels(obs,cfg,data_root):
 column=cfg.get('pixel_input_column')
 if not column:return obs,pd.DataFrame()
 source=pd.read_csv(data_root/cfg['input']);result=obs.copy();cache={};checks=[]
 for idx,row in result.iterrows():
  if int(row.input_row)<0:continue
  path=source.iloc[int(row.input_row)].get(column,'')
  if not isinstance(path,str) or 'pixel_exports/' not in path:continue
  if path not in cache:
   frame=pd.read_csv(data_root/path)
   if cfg.get('pixel_recipe')=='calibration_common_positive_finite':
    req=['Int','G','S','TauPhase','TauModulation'];valid=np.isfinite(frame[req].to_numpy(float)).all(axis=1)&(frame.Int>0);frame=frame[valid]
   cache[path]=frame
  frame=cache[path];metric=row.metric
  source_row=source.iloc[int(row.input_row)]
  if 'component_id' in source_row and pd.notna(source_row.component_id):
   frame=frame[frame.Xcoord.between(source_row.component_bbox_min_x,source_row.component_bbox_max_x)&frame.Ycoord.between(source_row.component_bbox_min_y,source_row.component_bbox_max_y)]
   if len(frame)!=int(source_row.pixel_count):raise ValueError('Lumen ROI pixel count differs from its recorded annotation')
  values=pd.to_numeric(frame[metric],errors='coerce').to_numpy(float);values=values[np.isfinite(values)]
  value=float(np.average(frame[metric],weights=frame.Int)) if metric in ['G','S'] and cfg.get('pixel_recipe')=='calibration_common_positive_finite' else float(np.mean(values))
  if metric in cfg.get('measurement_rounding',{}):value=float(np.round(value,cfg['measurement_rounding'][metric]))
  passed=bool(np.isclose(value,row.value,rtol=1e-8,atol=1e-10))
  checks.append(dict(input_file=path,metric=metric,reference=float(row.value),recomputed=value,passed=passed))
  result.at[idx,'value']=value
 checks=pd.DataFrame(checks)
 if len(checks) and not checks.passed.all():
  bad=checks[~checks.passed]
  raise ValueError(f'{len(bad)} pixel-to-image discrepancies; first: '+str(bad.iloc[0].to_dict()))
 return result,checks

def from_measurement_masks(obs, cfg, data_root):
 """Recalculate compartment means with deposited masks and explicit coordinate origins."""
 if not cfg.get('measurement_mask_index'):
  return obs, pd.DataFrame()
 from PIL import Image
 index = pd.read_csv(data_root / cfg['measurement_mask_index'], dtype={'condition': str, 'sample': str})
 result = obs.copy()
 checks = []
 for _, record in index.iterrows():
  mask = np.asarray(Image.open(data_root / record['mask']))
  shape = (int(record.height), int(record.width))
  if mask.shape != shape or mask.dtype != np.uint8 or not np.isin(mask, [0, 255]).all():
   raise ValueError(f'Invalid binary mask shape or values: {record["mask"]}')
  frame = pd.read_csv(data_root / record.nadph)
  coords = frame[['Ycoord', 'Xcoord']].to_numpy(float)
  if not np.isfinite(coords).all() or not np.equal(coords, np.floor(coords)).all():
   raise ValueError('Non-finite or non-integer image coordinates')
  coords = coords.astype(int) - int(record.coordinate_base)
  if (coords < 0).any() or (coords[:, 0] >= shape[0]).any() or (coords[:, 1] >= shape[1]).any():
   raise ValueError('Coordinates exceed the recorded measurement-mask extent')
  if frame[['Ycoord', 'Xcoord']].duplicated().any():
   raise ValueError('Duplicate NAD(P)H coordinates are not permitted')
  selected = mask[coords[:, 0], coords[:, 1]] > 0
  if int(selected.sum()) != int((mask > 0).sum()):
   raise ValueError('The measurement mask contains pixels absent from the NAD(P)H export')
  rows = result[(result.condition == record.condition) & (result['sample'].astype(str) == str(record['sample']))]
  for idx, row in rows.iterrows():
   values = frame.loc[selected, row.metric].to_numpy(float)
   value = float(np.nanmean(values))
   passed = bool(np.isclose(value, row.value, rtol=1e-8, atol=1e-10))
   checks.append(dict(condition=row.condition, sample=row['sample'], mask=record['mask'], input_file=record.nadph, metric=row.metric, reference=float(row.value), recomputed=value, selected_pixels=int(selected.sum()), passed=passed))
   result.at[idx, 'value'] = value
 checks = pd.DataFrame(checks)
 if len(checks) and not checks.passed.all():
  raise ValueError('Mask-to-measurement mismatch: ' + str(checks[~checks.passed].iloc[0].to_dict()))
 return result, checks
