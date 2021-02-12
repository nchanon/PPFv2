#!/usr/bin/env python

import os
import sys
sys.path.append('./')



import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

nbin = int(raw_input("number of bin : "))

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

histograms = []
data_number_of_event = []

def getTH1F(input, target):
    return input.Get(target)

data_file = TFile('./results/'+year+'/flattree/'+observable+'.root')
for l in data_file.GetListOfKeys():
    if TString(l.GetName()).Contains('grp'):
        hist = data_file.Get(l.GetName())
        hist.Scale(1./nbin)
        histograms.append(hist)

out = './combine/'+year+'/one_bin/inputs/'
for n in range(nbin):

    time_file = TFile('./results/'+year+'/flattree/'+observable+'_data_timed'+str(nbin)+'.root')
    histograms_timmed = time_file.Get('data_obs_bin'+str(n))
    data_number_of_event.append(histograms_timmed.GetEntries())

    output = TFile(out+observable+'_'+str(nbin)+'_'+str(n)+'.root', "RECREATE")
    histograms_timmed.Write()
    for h in histograms:
        h.Write()
    output.Close()


file_txt = ''
for i in data_number_of_event:
    file_txt += str(i)+'\n'

file = open('./combine/'+year+'/one_bin/inputs/'+observable+'_integrals_data_timmed.txt','w') 
file.write(file_txt) 
file.close() 