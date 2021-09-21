import os
import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

observable = [
    ['m_dilep', '\"dilepton mass\"'],
    #['n_bjets', '\"b-jets multiplicity\"'],
    #['pt_lead', '\"Leading lepton pt\"']
]

year = [
    #'2016',
    '2017'
]

for y in year:
    for o in observable:
        for s in systematic_list:
            cmd = 'python ./bin/systematics_observable.py '+o[0]+' '+y+' '+s+' '+o[1]
            os.system(cmd)