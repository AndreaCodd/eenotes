#!/usr/bin/python3
import json, os 
import shutil

COPY=[ '00_Intro.ipynb', '01_PythonBasics.ipynb', '02_Arrays.ipynb', '11_FEM.ipynb', '21_Gravity.ipynb', '22_Magnetic.ipynb',
   '22_Magnetotellurics.ipynb', '23_WavesFrequency.ipynb', '26_WavesTime.ipynb', '41_InversionAsLSQ.ipynb', '42_inversionGrav.ipynb',
   '43_inversionMT.ipynb', '44_inversionJoint.ipynb', 'mytools.py' ]
DCOPY = [ 'IMAGE_CS1.png', 'IMAGE_Dyke.jpg', 'IMAGE_Indicator_function.png', 'IMAGE_filter_function.png', 'IMAGE_fitting_grid.png',
  'Image_minima.png', 'MTImage1.png', 'SonicImage1.png', 'gravity_tasmania.csv', 'marmousi.bin' ]


SRCDIR='../Notebooks'
RELEASEDIR='../../material/Notebooks'

ANSWER= "SOLUTION"

def rmAnswers(src, dst):
    text=json.load(open(src, 'r'))
    for cell in text['cells']:
        if 'tags' in cell['metadata']:
            if cell['metadata']['tags'].count(ANSWER)>0:
                if 'source' in cell:
                    cell['source']=[]
                if 'outputs' in cell:
                    cell['outputs']=[]
 
    json.dump(text, open(dst, 'w'), indent=3) 


for fn in COPY + [ os.path.join("Data", f) for f in DCOPY]:
    frm=os.path.join(SRCDIR, fn)
    to=os.path.join(RELEASEDIR, fn)
    if os.path.splitext(fn)[1] == ".ipynb":
       rmAnswers(frm, to)
       print("cleaned ", frm," copied to", to)
    else:
       shutil.copy(frm, to)
       print(frm," copied to", to)
