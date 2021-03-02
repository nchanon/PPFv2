#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

import tools.tdrstyle as tdr
tdr.setTDRStyle()

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('title', help='display your observable title')
parser.add_argument('systematic',nargs='?', help='display your systematic', default='')
parser.add_argument('option',nargs='?', help='display rootfile write option', default='RECREATE')

args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title
systematic = args.systematic
option = args.option


nbin = 0
min_bin = 0
max_bin = 0
legend_coordinates = observable_values(observable)[1]
TH1.SetDefaultSumw2(1)
signal_integral = 0
background_integral_i = []
background_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('stack_'+observable,'stack_'+observable, 600, 600)

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

# signal
hist_signal = rootfile_input.Get('signal')
signal_integral = hist_signal.Integral()


# background
background_integral = 0
hist_background  = TH1F("background", "background", nbin, min_bin, max_bin)
for l in rootfile_input.GetListOfKeys():
    for s in ttbar_list:
        if(l.GetName() == s and not l.GetName() == 'signal'):
            print l.GetName()
            hist_background.Add(rootfile_input.Get(l.GetName()))
            background_integral_i.append([rootfile_input.Get(l.GetName()).Integral(), l.GetName()])
            background_integral += rootfile_input.Get(l.GetName()).Integral()


################################################################################
## Legend stuff
################################################################################

legend_args = (0.645, 0.79, 0.985, 0.91, '', 'NDC')
legend = TLegend(*legend_args)
legend.AddEntry(hist_signal, "t#bar{t} signal", "f")
legend.AddEntry(hist_background, "non-t#bar{t}", "f")
legend.AddEntry(hist_data, "data")
legend_box(legend, legend_coordinates)

################################################################################
## Draw stuff
################################################################################

stack = THStack()
stack.Add(hist_background)
stack.Add(hist_signal)
stack.Draw("E HIST")
hist_data.Draw("E SAME")
legend.Draw("SAME")

################################################################################
## Set Style
################################################################################

# line_color, line_width, fill_style, marker_style
style_histo(hist_signal, 2, 1, 2, 3004, 0)
style_histo(hist_background, 4, 1, 4, 3005, 0)
style_histo(hist_data, 1, 1, 0, 3001, 1, 20)

tdr.cmsPrel(41530., 13.)
style_labels_counting(stack, 'Events', title)

################################################################################
## Save
################################################################################

resultname = './results/'+year+'/comparaison/'+observable

rootfile_output = TFile(resultname+'.root', option)
canvas.Write()
canvas.SaveAs(resultname+'.png')
rootfile_output.Close()


print 'Data integral       : ', data_integral
print 'Signal integral     : ', signal_integral
print 'Background integral : ', background_integral
print ''
for i in background_integral_i:
    print i[1],i[0]

if(systematic != ''):
    print "Systematic up       : ", syst_up_integral
    print "Systematic down     : ", syst_down_integral
print ''

exit = raw_input("Press key to quit : ") 
