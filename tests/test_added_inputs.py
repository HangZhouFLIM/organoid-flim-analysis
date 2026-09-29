"""Tests for newly supplied pixel precision, object lineage and exclusion masks."""
import unittest,tempfile,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from PIL import Image
from workflows.object_measurements import from_object_exports,verify_object_masks
from workflows.pixel_measurements import from_pixels
class AddedInputTests(unittest.TestCase):
 def test_pixel_rounding_preserves_archived_precision(self):
  with tempfile.TemporaryDirectory() as tmp:
   r=Path(tmp);(r/'pixel_exports').mkdir();pd.DataFrame({'Int':[1.23454,1.23456]}).to_csv(r/'pixel_exports/p.csv',index=False);pd.DataFrame({'px':['pixel_exports/p.csv']}).to_csv(r/'i.csv',index=False);obs=pd.DataFrame([{'input_row':0,'metric':'Int','value':1.2346}]);cfg={'input':'i.csv','pixel_input_column':'px','measurement_rounding':{'Int':4}}
   _,checks=from_pixels(obs,cfg,r);self.assertTrue(checks.passed.all())
   cfg.pop('measurement_rounding')
   with self.assertRaisesRegex(ValueError,'pixel-to-image'):from_pixels(obs,cfg,r)
 def test_wrong_object_identity_is_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   r=Path(tmp);pd.DataFrame({'Id':[7],'MeanInt':[2.]}).to_csv(r/'o.csv',index=False);inp=pd.DataFrame([dict(source_object_file='o.csv',source_object_row=0,source_object_id=7,source_acquisition='sample1',value=2.)]);inp.to_csv(r/'i.csv',index=False);obs=pd.DataFrame([dict(input_row=0,metric='Int',value=2.)]);cfg={'input':'i.csv','object_input_column':'source_object_file'}
   _,checks=from_object_exports(obs,cfg,r);self.assertTrue(checks.passed.all());inp.loc[0,'source_object_id']=8;inp.to_csv(r/'i.csv',index=False)
   with self.assertRaisesRegex(ValueError,'identifier'):from_object_exports(obs,cfg,r)
 def test_exclusion_mask_is_applied_and_geometry_checked(self):
  with tempfile.TemporaryDirectory() as tmp:
   r=Path(tmp);lab=np.array([[1,1],[2,2]],np.uint16);exc=np.array([[0,1],[0,0]],np.uint8);final=np.where(exc==0,lab,0).astype(np.uint16)
   for name,arr in [('labels',lab),('exclude',exc),('final',final)]:Image.fromarray(arr).save(r/(name+'.png'))
   pd.DataFrame({'Id':[1,2],'Area':[1,2],'CentroidX':[1.,1.5],'CentroidY':[1.,2.]}).to_csv(r/'o.csv',index=False)
   pd.DataFrame([dict(sample='s1',source_labels='labels.png',exclusion_mask='exclude.png',mask='final.png',recipe='labels_where_exclusion_zero',object_file='o.csv',centroid_coordinate_base=1)]).to_csv(r/'idx.csv',index=False);cfg={'object_mask_index':'idx.csv'}
   self.assertTrue(verify_object_masks(cfg,r).passed.all());Image.fromarray(lab).save(r/'final.png')
   with self.assertRaisesRegex(ValueError,'Rebuilt cell mask'):verify_object_masks(cfg,r)
if __name__=='__main__':unittest.main()
