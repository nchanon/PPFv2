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
parser.add_argument('syst', help='display your systematics')


args = parser.parse_args()
observable = args.observable
year = args.year
syst = args.syst

TH1.SetDefaultSumw2(1)
canvas = TCanvas()

################################################################################
## function
################################################################################

mc_file = TFile('./results/'+year+'/flattree/'+observable+'.root')

h_nominal = mc_file.Get('signal')
nbin = h_nominal.GetNbinsX()
int_nom = h_nominal.Integral()
print int_nom

if syst == 'jec':
    infile = TFile('./results/'+year+'/flattree/'+observable+'_jec.root')
    h_up   =  infile.Get('signal_TotalUp')
    h_down =  infile.Get('signal_TotalDown')

else:
    infile = TFile('./results/'+year+'/flattree/'+observable+'_color_reco.root')
    h_up   =  infile.Get('signal_'+syst+'Up')
    h_down =  infile.Get('signal_'+syst+'Down')

h_up.SetLineColor(2)
h_down.SetLineColor(3)


h = h_nominal
h1 = h_nominal

h_up.Divide(h)
h_down.Divide(h1)
h.Divide(h_nominal)

h_nominal.SetStats(0)
h_nominal.GetXaxis().SetTitle('ratio '+syst+'/nominal')



legend_args = (0.645, 0.79, 0.985, 0.91, '', 'NDC')
legend = TLegend(*legend_args)
legend.AddEntry(h_up, syst+" up", "f")
legend.AddEntry(h_down,  syst+" down", "f")


h_nominal.Draw('HIST')
h_nominal.GetYaxis().SetRangeUser(0.92,1.08)
h_up.Draw('SAME HIST')
h_down.Draw('SAME HIST')
legend.Draw('SAME')

output = './results/'+year+'/other/'+observable+'_'+syst

canvas.SaveAs(output+'_'+year+'_check.pdf')


raw_input('quit')
