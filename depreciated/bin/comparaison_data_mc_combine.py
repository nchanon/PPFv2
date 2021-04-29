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
parser.add_argument('title', help='display your observable title')
parser.add_argument('file', help='input file')
parser.add_argument('systematic',nargs='?', help='display your systematic', default='')
parser.add_argument('subtitle',nargs='?', help='display your observable subtitle', default='')

args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title
systematic = args.systematic
subtitle = args.subtitle
inputfile = args.file
if(subtitle == ''):
    if year=='2016':
        subtitle = '35.9 fb^{-1} (13TeV)'
    elif year=='2017':
        subtitle = '41.53 fb^{-1} (13TeV)'

binning = observable_values(observable)[0]
legend_coordinates = observable_values(observable)[1]
TH1.SetDefaultSumw2(1)
mc_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('stack_'+observable,'stack_'+observable, 1000, 800)

################################################################################
## Create Histo 
################################################################################

ibin = inputfile[-6:-5]
hist_name = observable+'_bin_'+ibin

hist_signal =  TH1F(hist_name, hist_name, binning[0], binning[1], binning[2])
hist_background = TH1F("background","background", binning[0], binning[1], binning[2])
hist_data = TH1F("data","data", binning[0], binning[1], binning[2])

rootfile = TFile( inputfile )
hist_data = rootfile.Get("data_obs")
data_integral = hist_data.Integral()
for l in ttbar_list:
    tmp = rootfile.Get(l)
    hist_signal.Add(tmp)
    if(l.find('signal') != -1):
        continue
    else:
        hist_background.Add(tmp)
mc_integral = hist_signal.Integral()

################################################################################
## Draw stuff
################################################################################
'''
hist_mc = TH1F("","", nbin, min_bin, max_bin)
hist_mc.Add(hist_background)
hist_mc.Add(hist_signal)
#hist_mc.Draw("E HIST")
hist_background.Draw("E SAME")
hist_data.Draw("E SAME")
'''
hist_signal.Draw("E HIST")
hist_background.Draw("E HIST SAME")
hist_data.Draw("E SAME")

################################################################################
## Set Style
################################################################################

# line_color, line_width, fill_style, marker_style
style_histo(hist_signal, 2, 1, 2, 3004, 0)
style_histo(hist_background, 4, 1, 4, 3005, 0)
style_histo(hist_data, 1, 1, 0, 3001, 1, 20)
if(systematic != ''):
    style_histo(hist_systematics_up, 6, 2, 0, 1000, 0)
    style_histo(hist_systematics_down, 6, 2, 0, 1000, 0)

################################################################################
## Save
################################################################################

result_name = results_path(year,'stack/combine',hist_name)

newfile = TFile(result_name+'.root', 'RECREATE')
canvas.Update()
canvas.Write()
newfile.Close()
canvas.SaveAs(result_name+'.png')

print ''
print "Data : ", data_integral
print "MC : ", mc_integral
print "Ratio DATA/MC : ", percent(data_integral,mc_integral),'%'
print ''

#exit = raw_input("Press key to quit : ") 

