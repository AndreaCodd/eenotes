#!/usr/bin/python3
import json, os 
import shutil

COPY1=['01_PythonBasics.ipynb', '02_Arrays.ipynb', 'Data/gravity_tasmania.csv', 'Data/marmousi.bin' ]
COPY=[ '03_Gravity.ipynb', '04_Magnetic.ipynb', '05_Magnetotellurics.ipynb', '06_WavesFrequency.ipynb', 'Data/IMAGE_CS1.png', 'Data/IMAGE_Dyke.jpg', 
        'Data/MTImage1.png', 'esys-escript.pdf', 'Data/SonicImage1.png']


SRCDIR='../Notebooks'
RELEASEDIR='../material/Notebooks'

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


for fn in COPY:
    frm=os.path.join(SRCDIR, fn)
    to=os.path.join(RELEASEDIR, fn)
    if os.path.splitext(fn)[1] == ".ipynb":
       rmAnswers(frm, to)
       print("cleaned ", frm," copied to", to)
    else:
       shutil.copy(frm, to)
       print(frm," copied to", to)
