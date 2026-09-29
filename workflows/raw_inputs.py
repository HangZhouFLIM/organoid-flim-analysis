"""Validate deposited GSLab channel groups and pixel-to-acquisition links."""
from pathlib import Path
import csv,re,xml.etree.ElementTree as ET

def read_rows(path):
 with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def validate_raw_inputs(root,manifest):
 root=Path(root).resolve();files={r['file']:r for r in manifest};errors=[];groups={};seen=set();xml_checked=set()
 def check_file(name,role):
  p=(root/name).resolve()
  if not p.is_relative_to(root) or not p.is_file() or name not in files or files[name]['role']!=role:
   errors.append(f'{name}: missing, unlisted or incorrect processing role');return False
  return True
 for row in read_rows(root/'metadata/GSLab_input_index.csv'):
  prefix=row['acquisition_prefix'];channels=row['channel_images'].split(';');metadata=row['xml_metadata'].split(';')
  if prefix in groups:errors.append(prefix+': duplicate acquisition group')
  groups[prefix]=(set(channels),set(metadata));indices=[]
  for name in channels:
   match=re.fullmatch(re.escape(prefix)+r'_ch(\d+)\.tiff?',name,re.I)
   if not match:errors.append(name+': channel prefix mismatch')
   else:indices.append(int(match[1]))
   check_file(name,'raw_for_GSLab');seen.add(name)
  if not indices or len(indices)!=len(set(indices)) or any(any(i not in indices for i in range(k//3*3,k//3*3+3)) for k in indices):errors.append(prefix+': incomplete or duplicate channel triplet')
  for name in metadata:
   if check_file(name,'raw_for_GSLab') and name not in xml_checked:
    try:
     try:ET.parse(root/name)
     except ET.ParseError:ET.fromstring((root/name).read_bytes().decode('cp1252'))
    except (ET.ParseError,UnicodeError) as exc:errors.append(name+': invalid XML '+str(exc))
    xml_checked.add(name)
 expected={r['file'] for r in manifest if r['role']=='raw_for_GSLab' and r['file'].lower().endswith(('.tif','.tiff'))}
 if seen!=expected:errors.append('Raw channel index does not cover exactly the deposited channel images')
 expected_xml={r['file'] for r in manifest if r['role']=='raw_for_GSLab' and r['file'].lower().endswith('.xml')}
 if xml_checked!=expected_xml:errors.append('Raw acquisition index does not cover exactly the deposited XML files')
 pixel_seen=set()
 for row in read_rows(root/'metadata/pixel_to_raw_acquisition.csv'):
  name=row['pixel_file'];check_file(name,'pixel_exports')
  if name in pixel_seen:errors.append(name+': duplicate pixel-to-raw assignment')
  pixel_seen.add(name);group=groups.get(row['acquisition_prefix']);channels=set(row['channel_images'].split(';'));metadata=set(row['xml_metadata'].split(';'))
  if group is None or channels!=group[0] or not metadata.issubset(group[1]):errors.append(name+': raw acquisition crosswalk mismatch')
  for raw in channels|metadata:
   if raw in files and not set(row['panels'].split(';')).issubset(set(files[raw]['panels'].split(';'))):errors.append(name+': panel not linked to raw input '+raw)
 return errors
