#!/usr/bin/env python

import os
import sys
sys.path.append('./')


import argparse

from tools.sample_manager import *
from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

nbin = 24

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
observable = args.observable
year = args.year

TH1.SetDefaultSumw2(1)


################################################################################
## Code body
################################################################################

lumi_file = TFile('./inputs/timed/LumiUncertainties_'+year+'.root')
lumi_histo = lumi_file.Get('hInstLumi_DataScaleFactor')

#lumiCorrEvent = 0 #No lumi correction
lumiCorrEvent = 1 #Per-sidereal time bin lumi correction from lumi studies
#lumiCorrEvent = 2 #Per-event lumi correction using tree lumi

if (lumiCorrEvent==0 or lumiCorrEvent==1):
    data_file = TFile('./results/'+year+'/flattree/'+observable+'_data_timed24.root')
if (lumiCorrEvent==2):
    data_file = TFile('./results/'+year+'/flattree/'+observable+'_lumicorrected_data_timed24.root')

results = TH1F('data_time', 'data_time', nbin, 0, nbin)
for i in range(nbin):
    hist = data_file.Get('data_obs_bin'+str(i))
    content = hist.GetEntries()
    if (lumiCorrEvent==1): 
	content = hist.GetEntries()/lumi_histo.GetBinContent(i+1)
    print content - hist.GetEntries()
    results.SetBinContent(i+1, content)

if (lumiCorrEvent==0):
    output = TFile('./results/'+year+'/flattree/sideraltime_data.root', 'RECREATE')
if (lumiCorrEvent==1):
    output = TFile('./results/'+year+'/flattree/sideraltime_globallumicorr_data.root', 'RECREATE')
if (lumiCorrEvent==2):
    output = TFile('./results/'+year+'/flattree/sideraltime_lumicorrected_data.root', 'RECREATE')


results.Write()
