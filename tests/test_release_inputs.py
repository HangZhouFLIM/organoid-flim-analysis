"""Release-index and raw-acquisition integrity checks for the portable reader."""
import csv,json,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from workflows.raw_inputs import validate_raw_inputs

def write_csv(path,rows):
 with path.open('w',newline='',encoding='utf-8') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

class RawAcquisitionTests(unittest.TestCase):
 def fixture(self,root):
  (root/'metadata').mkdir();manifest=[];groups=[];links=[]
  xml='FigS08/raw_for_GSLab/MetaData/series.xml';(root/xml).parent.mkdir(parents=True);(root/xml).write_text('<Acquisition/>');manifest.append(dict(file=xml,role='raw_for_GSLab',panels='FigS08g-i'))
  for frame in [0,1]:
   prefix=f'FigS08/raw_for_GSLab/series_t{frame}';channels=[]
   for channel in [0,1,2]:
    name=f'{prefix}_ch{channel}.tif';(root/name).write_bytes(b'fixture');channels.append(name);manifest.append(dict(file=name,role='raw_for_GSLab',panels='FigS08g-i'))
   pixel=f'FigS08/pixel_{frame}.csv';(root/pixel).write_text('Int\n1\n');manifest.append(dict(file=pixel,role='pixel_exports',panels='FigS08g-i'))
   group=dict(acquisition_prefix=prefix,channel_images=';'.join(channels),xml_metadata=xml);groups.append(group);links.append(dict(pixel_file=pixel,panels='FigS08g-i',**group))
  write_csv(root/'metadata/GSLab_input_index.csv',groups);write_csv(root/'metadata/pixel_to_raw_acquisition.csv',links)
  return manifest,groups,links
 def test_time_frames_share_original_series_xml(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);manifest,_,_=self.fixture(root);self.assertEqual(validate_raw_inputs(root,manifest),[])
 def test_incomplete_channel_triplet_is_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);manifest,groups,_=self.fixture(root);groups[0]['channel_images']=';'.join(groups[0]['channel_images'].split(';')[:2]);write_csv(root/'metadata/GSLab_input_index.csv',groups)
   self.assertTrue(any('incomplete' in s for s in validate_raw_inputs(root,manifest)))
 def test_crosswalk_cannot_substitute_another_frame(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);manifest,_,links=self.fixture(root);links[0]['channel_images']=links[1]['channel_images'];write_csv(root/'metadata/pixel_to_raw_acquisition.csv',links)
   self.assertTrue(any('crosswalk mismatch' in s for s in validate_raw_inputs(root,manifest)))

class ReaderCommandTests(unittest.TestCase):
 def run_reader(self,folder,*args):
  return subprocess.run([sys.executable,'-B',str(ROOT/'reproduce.py'),*args],cwd=folder,capture_output=True,text=True,encoding='utf-8')
 def test_panel_list_needs_no_data_and_no_output_folder(self):
  with tempfile.TemporaryDirectory() as tmp:
   result=self.run_reader(tmp,'--list-panels');self.assertEqual(result.returncode,0,result.stderr);panels=[line.split('\t')[0] for line in result.stdout.splitlines()];self.assertEqual(len(panels),43);self.assertIn('Fig01b-d',panels);self.assertEqual(list(Path(tmp).iterdir()),[])
 def test_invalid_panel_fails_before_reading_data(self):
  with tempfile.TemporaryDirectory() as tmp:
   result=self.run_reader(tmp,'--data-root',str(Path(tmp)/'absent'),'--panel','FigMissing');self.assertEqual(result.returncode,2);self.assertIn('Unknown panel',result.stderr);self.assertEqual(list(Path(tmp).iterdir()),[])
 def test_changed_release_indexes_are_rejected_before_analysis(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);release=json.loads((ROOT/'manifests/data_release.json').read_text())
   for row in release['metadata_files']:
    path=root/row['file'];path.parent.mkdir(parents=True,exist_ok=True);path.write_text('altered index\n')
   result=self.run_reader(tmp,'--data-root',str(root),'--verify-only');self.assertNotEqual(result.returncode,0);self.assertIn('Zenodo release metadata differs',result.stderr);self.assertIn('checksum mismatch',result.stderr);self.assertFalse((root/'results').exists())

if __name__=='__main__':unittest.main()

