#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString

import tools.tdrstyle as tdr
tdr.setTDRStyle()

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
year = args.year

observables = [
    'm_dilep',
#    'n_bjets',
#    'pt_lead',
#    'pt_sublead'
]

def complete_integrale(th1):
    return th1.Integral(0, int(th1.GetXaxis().GetXmax())+1)


AN_2016 = [
    ['witness', 'data_obs', 164396.0],
    ['witness','signal', 158230.9],
    ['witness','ttx', 1313.2],
    ['witness','singletop', 8110.2],
    ['witness','dibosons', 656.0],
    ['witness','wjets', 232.1],
    ['witness','zjets', 1558.0]
]

AN_2017 = [
    ['witness', 'data_obs', 207484.0],
    ['witness','signal', 196433.0],
    ['witness','ttx', 1564.5],
    ['witness','singletop', 9778.3],
    ['witness','dibosons', 608.1],
    ['witness','wjets', 372.1],
    ['witness','zjets', 1870.2]
]

rootfile_name = []

for l in observables:
    rootfile_name.append('./results/'+year+'/flattree/'+l+'.root')

################################################################################
## Core stuff
################################################################################

from tabulate import tabulate
count = []

for i in range(len(rootfile_name)):
    file = TFile(rootfile_name[i])
    
    for l in ['data_obs']+ttbar_list:
        foo = file.Get(l)
        count.append( [ observables[i], l, complete_integrale(foo) ] )

if year=='2016':
    print tabulate(AN_2016)
elif year=='2017':
    print tabulate(AN_2017)
print tabulate(count)