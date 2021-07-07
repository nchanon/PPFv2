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

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
#parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('title', help='display your observable title')

args = parser.parse_args()
#observable = args.observable
year = args.year
title = args.title

def integral_complete(histo, max_bin):
    return histo.Integral()+histo.GetBinContent(int(max_bin+1))+histo.GetBinContent(0)


nbin = 0
min_bin = 0
max_bin = 0


legend_coordinates = [0.65, 0.75, 0.87, 0.87] 
TH1.SetDefaultSumw2(1)
signal_integral = 0
background_integral_i = []
background_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('pu','pu', 800, 800)
canvas.UseCurrentStyle()

if year == '2016':
    rootfile_data = TFile('./inputs/other/MyDataPileupHistogram.root')
    rootfile_mc = TFile('./inputs/other/pileup.root')
    puname = '#TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8#RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1#MINIAODSIM'
    rescale = 0.000657612

elif year == '2017':
    rootfile_data = TFile('./inputs/other/pudistributions_data_'+year+'.root')
    rootfile_mc = TFile('./inputs/other/pudistributions_mc_'+year+'.root')
    puname = '#TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2#MINIAODSIM'
    rescale = 0.000772359

################################################################################
## Create Histo 
################################################################################


###########
# data part
###########
hist_data  = rootfile_data.Get('pileup')
hist_data.Scale(1./hist_data.Integral())
hist_mc    = rootfile_mc.Get(puname)
#hist_mc.Scale(1./rescale)
hist_mc.Scale(1./(hist_mc.Integral()))
hist_data.SetAxisRange(0,100)

#hist_mc.Scale(1./hist_mc.Integral())
#hist_data.Scale(1./hist_data.Integral())
################################################################################
## Legend stuff
################################################################################

legend_args = (0.645, 0.79, 0.985, 0.85, '', 'NDC')
legend = TLegend(*legend_args)
legend.AddEntry(hist_data, "Data", "l")
legend.AddEntry(hist_mc, "MC", "l")
legend_box(legend, legend_coordinates)

################################################################################
## Draw stuff
################################################################################

hist_data.Draw("E HIST")
hist_mc.Draw(" E HIST SAME")
legend.Draw("SAME")

################################################################################
## Set Style
################################################################################

is_center=True

hist_mc.GetYaxis().SetTitle('Events')
hist_mc.GetXaxis().SetTitle(title)
hist_mc.GetYaxis().SetMaxDigits(4)
if(is_center):
    hist_mc.GetXaxis().CenterTitle()
    hist_mc.GetYaxis().CenterTitle()

# line_color, line_width, fill_style, marker_style
#style_histo(hist_data, 2, 5, 1, 3001, 1, 1)
#style_histo(hist_mc,   4, 5, 1, 3001, 1, 1)
style_histo(hist_data, 2, 1, 2, 3004, 0)
style_histo(hist_mc,   4, 1, 4, 3005, 0)

if(year=='2016'):
    tdr.cmsPrel(35900., 13.)
elif(year=='2017'):
    tdr.cmsPrel(41530., 13.)

################################################################################
## Save
################################################################################

resultname = './results/'+year+'/other/pileup_'+year

rootfile_output = TFile(resultname+'.root', "RECREATE")
canvas.Write()
canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
rootfile_output.Close()

raw_input('exit')