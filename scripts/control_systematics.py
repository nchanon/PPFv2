import os
import sys
sys.path.append('./')

import argparse

from tools.style_manager import *
from tools.sample_manager import *

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')
parser.add_argument('timed', help='timed or inclusive')
parser.add_argument('timebin', help='display the time bin')

args = parser.parse_args()
year = args.year
timed = args.timed
timebin = args.timebin

observable = [
    #['m_dilep', '\"Dilepton mass (GeV)\"'],
    #['pt_emu', '\"Dilepton p_{T} (GeV)\"'],
    ['n_bjets', '\"b-jets multiplicity\"'],
    #['pt_lead', '\"Leading lepton pt\"']
]

#year = [
#    '2016',
    #'2017'
#]


large_syst_list = []

#Complete list
large_syst_list = systematic_list

i=0
for jec_syst in jec_list[year]:
    if i%2==0:
        large_syst_list.append(jec_syst[:-3])
    i=i+1

#large_syst_list.append('Total')
#large_syst_list.append('Absolute')
#large_syst_list.append('Absolute_'+year)
#large_syst_list.append('FlavorQCD')
#large_syst_list.append('BBEC1')
#large_syst_list.append('BBEC1_'+year)
#large_syst_list.append('RelativeBal')
#large_syst_list.append('RelativeSample_'+year)
#large_syst_list.append('CP5')
large_syst_list.append('hdamp')
large_syst_list.append('erd')
large_syst_list.append('QCD')
large_syst_list.append('GluonMove')
large_syst_list.append('mtop')
#Adding timed uncertainties
if(timed=="timed"):
    large_syst_list.append('emu_trig_'+year)
    large_syst_list.append('lumi_stability_'+year)
    large_syst_list.append('lumi_linearity_'+year)

#Tests
#large_syst_list = ['syst_qcdscale', 'syst_pdfas']
#large_syst_list = ['syst_prefiring']
#large_syst_list = ['Total']
#large_syst_list = ['CP5', 'hdamp', 'mtop']
#large_syst_list = ['syst_pt_top']
#large_syst_list = ['mtop']
#large_syst_list =  ['erd','QCD','GluonMove']


for y in [year]:
    for o in observable:
        for s in large_syst_list: #['Total']:
            cmd = 'python ./bin/systematics_observable.py '+o[0]+' '+y+' '+timed+' '+s+' '+o[1]+' '+timebin
            print cmd
            os.system(cmd)
