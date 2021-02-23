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

###########
# data part
###########

data_file = TFile('./results/'+year+'/flattree/'+observable+'_data_timed24.root')
data_integral = 0
hist_data_in = []

for i in range(nbin): 
    hist_data_in.append(data_file.Get('data_obs_bin'+str(i)))
    data_integral += hist_data_in[i].Integral()

# convenient variables
binning = hist_data_in[0].GetNbinsX()
min_bin = hist_data_in[0].GetXaxis().GetXmin()
max_bin = hist_data_in[0].GetXaxis().GetXmax()

canvas_data = TCanvas('data_obs', 'data_obs', 1000, 800)
hist_data = TH1F('data_obs','data_obs', binning*nbin,  0, nbin)
for iobs in range(nbin):
    for n in range(binning):
        hist_data.SetBinContent(n + iobs*binning + 1 , hist_data_in[iobs].GetBinContent(n))

###########
# MC part
###########

mc_file = TFile('./results/'+year+'/flattree/'+observable+'.root')
mc_integral = 0
canvas  = []
hist_mc = []

index = 0
for l in mc_file.GetListOfKeys():
    if not TString(l.GetName()).Contains('data_obs'):
        hist = mc_file.Get(l.GetName())
        hist_mc.append(TH1F("","", binning*nbin,  0, nbin))
        for iobs in range(nbin):
            for j in range(binning):
                hist_mc[index].SetBinContent(j + iobs*binning + 1 , hist.GetBinContent(j+1))
        del hist
        hist_mc[index].SetName(l.GetName())
        hist_mc[index].SetTitle(l.GetName())
        hist_mc[index].Scale(1./nbin)
        for g in ttbar_list:
            if TString(l.GetName()) == g:
                mc_integral += hist_mc[index].Integral()
        index += 1


###########
# SME part
###########

cmunu = 0.01
wilson = 'fXX_L'

sme_file = TFile('./results/'+year+'/flattree/sme.root')
sme_sig = sme_file.Get(wilson)

for g in hist_mc:
    if g.GetName() == 'signal':
        hist_sme = g.Clone()
hist_sme.SetName(wilson)
hist_sme.SetTitle(wilson)
for i in range(nbin):
    for j in range(binning):
        hist_sme.SetBinContent(j + i*binning + 1, 
                                    hist_sme.GetBinContent(j + i*binning + 1 )
                                    *(1+cmunu*sme_sig.GetBinContent(i+1))
                              )

print 'data = '+str(data_integral)
print 'mc   = '+str(mc_integral)

out = './combine/'+year+'/unrolled/inputs/'
output = TFile(out+observable+'.root', "RECREATE")
hist_data.Write()
hist_sme.Write()
for l in hist_mc:
    l.Write()
output.Close()

cmd = 'cp '+out+observable+'.root '+'./combine/'+year+'/interference/inputs/'
os.system(cmd)

file = open('./combine/'+year+'/'+observable+'_noe_data.txt','w') 
file.write(str(data_integral)) 
file.close() 