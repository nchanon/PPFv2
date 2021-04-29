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
parser.add_argument('subtitle',nargs='?', help='display your observable subtitle', default='')

args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title
subtitle = args.subtitle
if(subtitle == ''):
    if year=='2016':
        subtitle = '35.9 fb^{-1} (13TeV)'
    elif year=='2017':
        subtitle = '41.53 fb^{-1} (13TeV)'
bin_values = observable_values(observable)[0]
legend_coordinates = observable_values(observable)[1]


TH1.SetDefaultSumw2(1)
hist_signal = TH1F("","", bin_values[0], bin_values[1], bin_values[2])
hist_background = TH1F("","", bin_values[0], bin_values[1], bin_values[2])
hist_data = TH1F("","", bin_values[0], bin_values[1], bin_values[2])
canvas = TCanvas('stack_'+observable,'stack_'+observable, 1000, 800)
stack = THStack()

################################################################################
## Create Histo 
################################################################################

mc_integral = 0
data_integral = 0

# mc_signal
rootfile = []
for i in range(number_of_signal_samples[year]):
    namefile = sample_list['MC'][year][i]+'_'+observable+'.root'
    rootfile.append(TFile(results_path(year,'TH1/MC',namefile)))
    foo = rootfile[i].Get(observable)
    mc_integral += foo.Integral()
    hist_signal.Add(foo)
    print(sample_list['MC'][year][i])
del rootfile

# mc_background
rootfile = []
for i in range(number_of_signal_samples[year], len(sample_list['MC'][year])):
    namefile = sample_list['MC'][year][i]+'_'+observable+'.root'
    rootfile.append(TFile(results_path(year,'TH1/MC',namefile)))
    foo = rootfile[i-number_of_signal_samples[year]].Get(observable)
    mc_integral += foo.Integral()
    hist_background.Add(foo)
    print(sample_list['MC'][year][i])
del rootfile


# data
rootfile = []
for i in range(len(sample_list['DATA'][year])):
    namefile = sample_list['DATA'][year][i]+'_'+observable+'.root'
    rootfile.append(TFile(results_path(year,'TH1/DATA',namefile)))
    foo = rootfile[i].Get(observable)
    data_integral += foo.Integral()
    hist_data.Add(foo)
    print(sample_list['DATA'][year][i])
del rootfile


################################################################################
## Draw stuff
################################################################################

stack.Add(hist_background)
stack.Add(hist_signal)
stack.Draw("E HIST")
hist_data.Draw("E SAME")

################################################################################
## Ratio DATA/MC
################################################################################

ratio_limit = 0.4
ratio = TRatioPlot(stack, hist_data)

################################################################################
## Set Style
################################################################################

# line_color, line_width, fill_style, marker_style
style_histo(hist_signal, 2, 1, 2, 3004, 0)
style_histo(hist_background, 4, 1, 4, 3005, 0)
style_histo(hist_data, 1, 1, 0, 3001, 1, 20)

style_ratioplot(ratio, 'MC/DATA', ratio_limit)
style_labels_counting(stack, 'Events', title)

################################################################################
## Legend stuff
################################################################################

p = ratio.GetUpperPad()
legend = p.BuildLegend()
legend.Clear()
legend.SetHeader(subtitle, "C")
legend.AddEntry(hist_signal, "t#bar{t} signal", "f")
legend.AddEntry(hist_background, "non-t#bar{t}", "f")
legend.AddEntry(hist_data, "data")
legend_box(legend, legend_coordinates)

################################################################################
## Save
################################################################################

newfile = TFile(results_path(year,'stack',observable+'.root'), 'RECREATE')
canvas.Update()
canvas.Write()
newfile.Close()
canvas.SaveAs(results_path(year,'stack',observable+'.png'))

print ''
print "Data : ", data_integral
print "MC : ", mc_integral
print "Ratio MC/DATA : ", percent(mc_integral,data_integral),'%'
print ''

exit = raw_input("Press key to quit : ") 

