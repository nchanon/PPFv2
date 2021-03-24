#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad


################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('title', help='display your observable title')
parser.add_argument('systematic', help='display your systematic')
parser.add_argument('subtitle',nargs='?', help='display your observable subtitle', default='')
parser.add_argument('option',nargs='?', help='display rootfile write option', default='RECREATE')


args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title
systematic = args.systematic
subtitle = args.subtitle
if(subtitle == ''):
    if year=='2016':
        subtitle = '35.9 fb^{-1} (13TeV)'
    elif year=='2017':
        subtitle = '41.53 fb^{-1} (13TeV)'
option = args.option

nbin = 0
min_bin = 0
max_bin = 0
legend_coordinates = observable_values(observable)[1]
TH1.SetDefaultSumw2(1)
canvas = TCanvas('stack_'+observable,'stack_'+observable, 1000, 800)

rootfile_input = TFile('./results/'+year+'/'+observable+'.root')

################################################################################
## Create Histo 
################################################################################


###########
# mc part
###########

#signal
hist_mc = rootfile_input.Get('grp_signal')
hist_mc.SetTitle(systematic)

# convenient variables
nbin    = hist_mc.GetNbinsX()
min_bin = hist_mc.GetXaxis().GetXmin()
max_bin = hist_mc.GetXaxis().GetXmax()

# background
hist_background  = TH1F("background", "background", nbin, min_bin, max_bin)
for l in rootfile_input.GetListOfKeys():
    if(TString(l.GetName()).Contains("grp") and not TString(l.GetName()).Contains("signal")  and not TString(l.GetName()).Contains("Up") and not TString(l.GetName()).Contains("Down")):
        hist_background.Add(rootfile_input.Get(l.GetName()))

hist_mc.Add(hist_background)
hist_mc.SetStats(0)

###########
# syst part
###########

hist_up = TH1F("up", "up", nbin, min_bin, max_bin)
hist_down = TH1F("down", "down", nbin, min_bin, max_bin)

for l in rootfile_input.GetListOfKeys():
    if(TString(l.GetName()).Contains("grp") and TString(l.GetName()).Contains("Up") and TString(l.GetName()).Contains(systematic)):
        hist_up.Add(rootfile_input.Get(l.GetName()))
        print l.GetName(), rootfile_input.Get(l.GetName()).Integral() 

    if(TString(l.GetName()).Contains("grp") and TString(l.GetName()).Contains("Down") and TString(l.GetName()).Contains(systematic)):
        hist_down.Add(rootfile_input.Get(l.GetName()))
        print l.GetName(), rootfile_input.Get(l.GetName()).Integral() 


rootfile_input.Close()

################################################################################
## Set Style
################################################################################

style_histo(hist_mc, 4, 2, 0, 1000, 1, 20)
style_histo(hist_up, 3, 2, 0, 1000, 0)
style_histo(hist_down, 6, 2, 0, 1000, 0)

style_labels_counting(hist_mc, 'Events', title)

################################################################################
## Legend stuff
################################################################################

legend = TLegend()
legend.Clear()
legend.SetHeader(subtitle, "C")
legend.AddEntry(hist_mc, 'MC')
legend.AddEntry(hist_up, 'up')
legend.AddEntry(hist_down, 'down')
legend_box(legend, legend_coordinates)

################################################################################
## Save
################################################################################

hist_mc.Draw()
hist_up.Draw("SAME")
hist_down.Draw("SAME")
legend.Draw("SAME")

print hist_mc.Integral()
print hist_up.Integral()
print hist_down.Integral()

canvas.SaveAs(observable+'_'+systematic+'.png')

exit = raw_input("Press key to quit : ") 
