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
## function
################################################################################



################################################################################
## Code body
################################################################################

# time dependant
lumi_syst_up = []
lumi_syst_down = []

histograms = []

data_file = TFile('./results/'+year+'/flattree/'+observable+'.root')
for l in data_file.GetListOfKeys():
    hist = data_file.Get(l.GetName())
    histograms.append(hist)

data_file2 = TFile('./results/'+year+'/flattree/'+observable+'_timed.root')
for l in data_file2.GetListOfKeys():
    hist = data_file2.Get(l.GetName())
    histograms.append(hist)

out = './combine/'+year+'/inclusive/inputs/'
output = TFile(out+observable+'.root', "RECREATE")
for h in histograms:
    h.Write()
output.Close()