#!/usr/bin/env python

import os
import sys
sys.path.append('./')


import argparse

from tools.sample_manager import *
from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

#nbin = 24

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
#lumi_syst_up = []
#lumi_syst_down = []

histograms = []


data_file = TFile('./results/'+year+'/flattree/'+observable+'_data.root')
for l in data_file.GetListOfKeys():
    hist = data_file.Get(l.GetName())
    histograms.append(hist)

mc_file = TFile('./results/'+year+'/flattree/'+observable+'_inclusive.root')
for l in mc_file.GetListOfKeys():
    hist = mc_file.Get(l.GetName())
    hname = hist.GetName()
    newname = hname
    if (hname.find("syst_em_trigUp")!=-1):
	newname = hname[:-2] + '_'+year+'Up'
    if (hname.find("syst_em_trigDown")!=-1):
        newname = hname[:-4] + '_'+year+'Down'
    hist.SetName(newname)
    hist.SetTitle(newname) 
    histograms.append(hist)

mc_alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco_inclusive.root')
for l in mc_alt_file.GetListOfKeys():
    hist = mc_alt_file.Get(l.GetName())
    histograms.append(hist)

mc_jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec_inclusive.root')
for l in mc_jec_file.GetListOfKeys():
    hist = mc_jec_file.Get(l.GetName())
    hname = hist.GetName()
    newname = hname
    if(hname.find('TotalUp')!= -1):
        newname = hname[:-7]+'jecUp'
    elif(hname.find('TotalDown')!= -1):
        newname = hname[:-9]+'jecDown'
    hist.SetName(newname)
    hist.SetTitle(newname)
    histograms.append(hist)


out = './combine/'+year+'/inclusive/inputs/'
output = TFile(out+observable+'_inclusive.root', "RECREATE")
for h in histograms:
    h.Write()
output.Close()

print 'Produced '+output.GetName()
