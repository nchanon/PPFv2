#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.directory_manager import *
from tools.sample_manager import *
from tools.style_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack 
from ROOT import TLegend, TApplication, TRatioPlot, TPad

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
observable = args.observable
year = args.year

bin_values = observable_values(observable)[0]

TH1.SetDefaultSumw2(1)

newfile = TFile(results_path(year,'combine','combine_'+observable+year+'.root'),'RECREATE')
newfile.Close()

################################################################################
## function
################################################################################

# DATA
canvas_data = TCanvas('data_obs', 'data_obs', 1000, 800)
hist_data = TH1F('data_obs','data_obs', bin_values[0], bin_values[1], bin_values[2])
rootfile = []
for i in range(len(sample_list['DATA'][year])):
    namefile = sample_list['DATA'][year][i]+'_'+observable+'.root'
    rootfile.append(TFile(results_path(year,'TH1/DATA',namefile)))
    foo = rootfile[i].Get(observable)
    hist_data.Add(foo)
del rootfile
rootfile =  TFile(results_path(year,'combine','combine_'+observable+year+'.root'),  'UPDATE')
hist_data.Draw()
hist_data.Write()
rootfile.Close()

# MC 
for l in ttbar_list:
    canvas = TCanvas(l, l)
    namefile = l+'_'+observable
    rfile = TFile(results_path(year,'groups/MC',namefile+'.root'))
    hist = rfile.Get(namefile)
    rootfile = TFile(results_path(year,'combine','combine_'+observable+year+'.root'), 'UPDATE')
    hist.Draw()
    hist.Write()
    rootfile.Close()


up_down = ['Up', 'Down']

#SYST 
for i, s in enumerate(up_down):
    for syst in systematic_list:
        for l in ttbar_list:
            canvas = TCanvas(namefile, namefile)
            namefile = l+'_'+observable+'_'+syst+s
            rfile = TFile(results_path(year,'groups/SYST',namefile+'.root'))
            hist = rfile.Get(namefile)
            rootfile =  TFile(results_path(year,'combine','combine_'+observable+year+'.root'), 'UPDATE')
            hist.Draw()
            hist.Write()
            rootfile.Close()
