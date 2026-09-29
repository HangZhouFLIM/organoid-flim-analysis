"""Integrity and coordinate-alignment regression tests for the public reader."""
import sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from PIL import Image
from reproduce import holm, sha, verify_data, expected_checks
from workflows.pixel_measurements import from_measurement_masks

class IntegrityTests(unittest.TestCase):
 def test_changed_input_is_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp); p=root/'input.csv'; p.write_text('value\n1\n')
   manifest=[{'file':'input.csv','sha256':sha(p)}]
   self.assertEqual(verify_data(root,manifest),[])
   p.write_text('value\n2\n')
   self.assertIn('checksum mismatch',verify_data(root,manifest)[0])
 def test_outside_path_is_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   self.assertIn('unsafe path',verify_data(Path(tmp),[{'file':'../external.csv','sha256':''}])[0])
 def test_holm_preserves_input_order_and_caps(self):
  np.testing.assert_allclose(holm([.04,.001,.03,.8]),[.09,.004,.09,.8])

class CoordinateTests(unittest.TestCase):
 def compute(self,root,base,wrong_origin=False,duplicate=False):
  mask=np.array([[255,0],[0,255]],np.uint8)
  Image.fromarray(mask).save(root/'mask.tif')
  frame=pd.DataFrame({'Xcoord':[base,base+1,base+1],'Ycoord':[base,base,base+1],'Int':[2.,100.,4.]})
  if duplicate:frame=pd.concat([frame,frame.iloc[[0]]],ignore_index=True)
  frame.to_csv(root/'pixels.csv',index=False)
  pd.DataFrame([dict(condition='control',sample=1,mask='mask.tif',nadph='pixels.csv',coordinate_base=base+1 if wrong_origin else base,height=2,width=2)]).to_csv(root/'index.csv',index=False)
  obs=pd.DataFrame([dict(condition='control',sample='1',metric='Int',value=3.)])
  return from_measurement_masks(obs,{'measurement_mask_index':'index.csv'},root)
 def test_zero_and_one_based_coordinates_match(self):
  for base in [0,1]:
   with tempfile.TemporaryDirectory() as tmp:
    obs,checks=self.compute(Path(tmp),base)
    self.assertEqual(obs.iloc[0].value,3.)
    self.assertEqual(checks.iloc[0].selected_pixels,2)
    self.assertTrue(checks.passed.all())
 def test_wrong_origin_is_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   with self.assertRaisesRegex(ValueError,'Coordinates exceed'):
    self.compute(Path(tmp),0,wrong_origin=True)
 def test_duplicate_coordinates_are_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   with self.assertRaisesRegex(ValueError,'Duplicate'):
    self.compute(Path(tmp),1,duplicate=True)
class ReferenceMatchingTests(unittest.TestCase):
 def test_equal_p_values_do_not_resolve_ambiguous_comparisons(self):
  actual=pd.DataFrame([dict(kind='posthoc',metric='Int',channel=c,condition='treated',p_adjusted=.04) for c in ['NADPH','TMRM']])
  reference=pd.DataFrame([dict(metric='Int',condition='treated',test='Dunn',p_adjusted=.04)])
  self.assertEqual(expected_checks(actual,reference,{'mode':'endpoint'}).iloc[0].status,'ambiguous')
  reference.loc[0,'metric']='TMRM Int'
  self.assertEqual(expected_checks(actual,reference,{'mode':'endpoint'}).iloc[0].status,'pass')
if __name__=='__main__':unittest.main()
