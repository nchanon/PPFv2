#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame

import tools.tdrstyle as tdr
tdr.setTDRStyle()

print systematic_list

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('systematic', help='display your systematic')
parser.add_argument('title', help='display your observable title')

args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title
systematic = args.systematic

nbin = 0
min_bin = 0
max_bin = 0
legend_coordinates = observable_values(observable)[1]
TH1.SetDefaultSumw2(1)
signal_integral = 0
mc_integral_i = []
mc_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('stack_'+observable,'stack_'+observable, 800, 800)
canvas.UseCurrentStyle()

rootfile_input = TFile('./results/'+year+'/flattree/'+observable+'.root')


################################################################################
## Create Histo 
################################################################################


###########
# data part
###########
hist_data = rootfile_input.Get('data_obs')
data_integral = hist_data.Integral()


# convenient variables
nbin    = hist_data.GetNbinsX()
min_bin = hist_data.GetXaxis().GetXmin()
max_bin = hist_data.GetXaxis().GetXmax()


###########
# mc part
###########

hist_mc = []
hist_mc_up = []
hist_mc_down = []


for l in rootfile_input.GetListOfKeys():
    for s in ttbar_list:
        if(l.GetName() == s):
            hist_mc.append(rootfile_input.Get(l.GetName()))
            mc_integral_i.append([rootfile_input.Get(l.GetName()).Integral(), l.GetName()])
        elif(TString(l.GetName()).Contains(s) and TString(l.GetName()).Contains(systematic) and TString(l.GetName()).Contains('Up')):
            name = s + '_' + systematic + 'Up'
            hist_mc_up.append(rootfile_input.Get(l.GetName()))
        elif(TString(l.GetName()).Contains(s) and TString(l.GetName()).Contains(systematic) and TString(l.GetName()).Contains('Down')):
            name = s + '_' + systematic + 'Down'
            hist_mc_down.append(rootfile_input.Get(l.GetName()))

for l in range(len(hist_mc)):
    print hist_mc[l].GetName(),  hist_mc_up[l].GetName(),  hist_mc_down[l].GetName() 



################################################################################
## Legend stuff
################################################################################

legend_args = (0.2, 0.79, 0.4
, 0.91, '', 'NDC')
legend = []
for index in range(len(hist_mc)):
    legend.append(TLegend(*legend_args))
    #legend[index].SetHeader(systematic)
    legend[index].AddEntry(hist_mc[index], hist_mc[index].GetName())
    legend[index].AddEntry(hist_mc_up[index], 'up')    
    legend[index].AddEntry(hist_mc_down[index], 'down')

    #legend_box(legend[index], legend_coordinates)

################################################################################
## Set Style
################################################################################

for index in range(len(hist_mc)):
    #style_histo(hist_mc_up[index], 7, 2, 0, 1000, 0)
    #style_histo(hist_mc_down[index], 6, 2, 0, 1000, 0)
    style_histo(hist_mc_up[index], 2, 1, 2, 3004, 0)
    style_histo(hist_mc_down[index], 4, 1, 4, 3005, 0)
    style_histo(hist_mc[index], 1, 1, 0, 3001, 1, 20)
    style_labels_counting(hist_mc[index], 'Uncertainty (in %)', title)


if(year=='2016'):
    tdr.cmsPrel(35900., 13.)
elif(year=='2017'):
    tdr.cmsPrel(41530., 13.)

################################################################################
## Save
################################################################################
import numpy as np

for h in range(len(hist_mc)):
    for i in range(hist_mc[h].GetNbinsX()):
        val = float(hist_mc[h].GetBinContent(i+1))
        if val == 0:
            val = 0.00001
        hist_mc[h].SetBinContent(i+1,0)
        hist_mc[h].SetBinError(i+1,hist_mc[h].GetBinError(i+1)/np.sqrt(val))

        up   = hist_mc_up[h].GetBinContent(i+1)-val
        down = hist_mc_down[h].GetBinContent(i+1)-val
        hist_mc_up[h].SetBinContent(i+1, float(up)/val*100)
        hist_mc_down[h].SetBinContent(i+1, float(down)/val*100)
    hist_mc[h].SetMaximum(hist_mc_up[h].GetMaximum()+hist_mc_up[h].GetMaximum()/10)
    hist_mc[h].SetMinimum(-hist_mc_up[h].GetMaximum()-hist_mc_up[h].GetMaximum()/10)

outputdir = './results/'+year+'/systematics/'

for index in range(len(hist_mc)):
    name = observable+'_'+hist_mc[index].GetName()+'_'+systematic
    canvas = TCanvas(name, name)
    hist_mc[index].Draw('')
    maxi = hist_mc[index].GetMaximum()
    hist_mc[index].GetYaxis().SetRangeUser(-maxi-2*maxi/3.,maxi+2*maxi/3.)
    hist_mc_up[index].Draw('HIST SAME')
    hist_mc_down[index].Draw('HIST SAME')
    legend[index].Draw('SAME')
    canvas.SaveAs(outputdir+name+'_'+year+'.png')
    canvas.SaveAs(outputdir+name+'_'+year+'.pdf')

#raw_input('exit')