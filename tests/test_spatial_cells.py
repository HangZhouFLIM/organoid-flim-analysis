from pathlib import Path
import unittest,sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from workflows.spatial_cells import project_cells,summarize_bins
class SpatialTests(unittest.TestCase):
 def test_clamped_projection_and_final_bin_match_recorded_convention(self):
  t=pd.DataFrame({'Id':[1,2,3,4,5,6],'CentroidX':[0,.2,10,-1,12,np.nan],'CentroidY':[0,0,0,0,0,np.nan]})
  position,bins=project_cells(t,1,3);np.testing.assert_allclose(position,[0,.02,1,0,1,0]);np.testing.assert_array_equal(bins,[1,2,50,1,50,1])
 def test_matlab_percentiles_and_singleton_sem(self):
  t=pd.DataFrame({'Bin':[1,1,1,1,2,2],'MeanInt':[1,2,3,4,np.nan,9]});s=summarize_bins(t,['MeanInt']).set_index('Bin');self.assertEqual(s.loc[1,'P25'],1.5);self.assertEqual(s.loc[1,'P75'],3.5);self.assertAlmostEqual(s.loc[1,'SEM'],np.std([1,2,3,4],ddof=1)/2);self.assertEqual(s.loc[2,'N'],1);self.assertEqual(s.loc[2,'SEM'],0);self.assertTrue(np.isnan(s.loc[3,'Mean']))
 def test_ambiguous_object_ids_and_degenerate_axis_are_rejected(self):
  t=pd.DataFrame({'Id':[1,1],'CentroidX':[0,1],'CentroidY':[0,0]})
  with self.assertRaisesRegex(ValueError,'Duplicate'):project_cells(t,1,1)
  t.Id=[1,2];t.CentroidX=[0,0]
  with self.assertRaisesRegex(ValueError,'coincide'):project_cells(t,1,2)
if __name__=='__main__':unittest.main()
