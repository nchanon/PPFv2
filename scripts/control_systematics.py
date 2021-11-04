import os
import sys
sys.path.append('./')

import argparse

from tools.style_manager import *
from tools.sample_manager import *

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
year = args.year

observable = [
    ['m_dilep', '\"Dilepton mass (GeV)\"'],
    #['n_bjets', '\"b-jets multiplicity\"'],
    #['pt_lead', '\"Leading lepton pt\"']
]

#year = [
#    '2016',
    #'2017'
#]


#large_syst_list = []
large_syst_list = systematic_list
large_syst_list.append('Total')
large_syst_list.append('CP5')
large_syst_list.append('hdamp')
large_syst_list.append('erd')
large_syst_list.append('QCD')
large_syst_list.append('GluonMove')
large_syst_list.append('mtop')

#large_syst_list = ['syst_prefiring']
#large_syst_list = ['Total']
#large_syst_list = ['CP5', 'hdamp', 'mtop']
#large_syst_list = ['syst_pt_top']
#large_syst_list = ['mtop']
#large_syst_list =  ['erd','QCD','GluonMove']

for y in [year]:
    for o in observable:
        for s in large_syst_list: #['Total']:
            cmd = 'python ./bin/systematics_observable_new.py '+o[0]+' '+y+' '+s+' '+o[1]
            print cmd
            os.system(cmd)
