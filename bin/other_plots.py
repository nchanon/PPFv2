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
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('title', help='display your observable title')

args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title

def integral_complete(histo, max_bin):
    return histo.Integral()+histo.GetBinContent(int(max_bin+1))+histo.GetBinContent(0)


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
data_integral = integral_complete(hist_data, max_bin)


# convenient variables
nbin    = hist_data.GetNbinsX()
min_bin = hist_data.GetXaxis().GetXmin()
max_bin = hist_data.GetXaxis().GetXmax()

###########
# mc part
###########

# signal
hist_signal = rootfile_input.Get('signal')
signal_integral = integral_complete(hist_signal, max_bin)

# background
background_integral = 0
hist_background  = TH1F("background", "background", nbin, min_bin, max_bin)
for l in rootfile_input.GetListOfKeys():
    for s in ttbar_list:
        if(l.GetName() == s and not l.GetName() == 'signal'):
            print l.GetName()
            hist_background.Add(rootfile_input.Get(l.GetName()))
            background_integral_i.append([integral_complete(rootfile_input.Get(l.GetName()), max_bin), l.GetName()])
            background_integral += integral_complete(rootfile_input.Get(l.GetName()), max_bin)



################################################################################
## Legend stuff
################################################################################

legend_args = (0.645, 0.79, 0.985, 0.91, '', 'NDC')
legend = TLegend(*legend_args)
legend.AddEntry(hist_signal, "t#bar{t}", "f")
legend.AddEntry(hist_background, "non-t#bar{t}", "f")
legend_box(legend, legend_coordinates)

################################################################################
## Draw stuff
################################################################################

hist_background.Scale(1./background_integral)
hist_signal.Scale(1./signal_integral)
hist_background.Draw("HIST")
hist_signal.Draw("HIST SAME")
legend.Draw("SAME")
#canvas.SetLogy()

################################################################################
## Set Style
################################################################################

# line_color, line_width, fill_style, marker_style
style_histo(hist_signal, 2, 1, 2, 3004, 0)
style_histo(hist_background, 4, 1, 4, 3005, 0)

style_labels_counting(hist_background, 'Events', title)
hist_signal.GetXaxis().SetLabelSize(0)
hist_signal.GetXaxis().SetTitleSize(0)

if(year=='2016'):
    tdr.cmsPrel(35900., 13.)
elif(year=='2017'):
    tdr.cmsPrel(41530., 13.)


################################################################################
## Save
################################################################################

resultname = './results/'+year+'/other/comp_'+observable+'_'+year

rootfile_output = TFile(resultname+'.root', "RECREATE")
canvas.Write()
canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
rootfile_output.Close()


print 'Data integral       : ', data_integral
print 'Signal integral     : ', signal_integral
print 'Background integral : ', background_integral
print ''
for i in background_integral_i:
    print i[1],i[0]

raw_input('exit')