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
parser.add_argument('systematic',nargs='?', help='display your systematic', default='')
parser.add_argument('subtitle',nargs='?', help='display your observable subtitle', default='')

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
canvas = TCanvas('stack_'+observable,'stack_'+observable, 1000, 800)

rootfile_input = TFile('./results/'+year+'/'+observable+'.root')

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
hist_signal = rootfile_input.Get('grp_signal')
signal_integral = hist_signal.Integral()


# background
hist_background  = TH1F("background", "background", nbin, min_bin, max_bin)
for l in rootfile_input.GetListOfKeys():
    if(TString(l.GetName()).Contains("grp") and not TString(l.GetName()).Contains("signal")  and not TString(l.GetName()).Contains("Up") and not TString(l.GetName()).Contains("Down")):
        print l.GetName()
        hist_background.Add(rootfile_input.Get(l.GetName()))
        background_integral_i.append([rootfile_input.Get(l.GetName()).Integral(), l.GetName()])
        background_integral += rootfile_input.Get(l.GetName()).Integral()

'''

## systematics
if(systematic != ''):
    hist_systematics_up = TH1F("","", nbin, min_bin, max_bin)
    hist_systematics_down = TH1F("","", nbin, min_bin, max_bin)

    rootfile = []
    for i, l in enumerate(ttbar_list):
        name = l+ '_'+observable+'_'+systematic+'Up'
        rootfile.append(TFile(results_path(year,'groups/SYST',name+'.root')))
        foo = rootfile[i].Get(name)
        syst_up_integral += foo.Integral()
        hist_systematics_up.Add(foo)
    del rootfile

    rootfile = []
    for i, l in enumerate(ttbar_list):
        name = l+ '_'+observable+'_'+systematic+'Down'
        rootfile.append(TFile(results_path(year,'groups/SYST',name+'.root')))
        foo = rootfile[i].Get(name)
        syst_down_integral += foo.Integral()
        hist_systematics_down.Add(foo)
    del rootfile
'''
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

stack = THStack()
stack.Add(hist_background)
stack.Add(hist_signal)
stack.Draw("E HIST")
hist_data.Draw("E SAME")
if(systematic != ''):
    hist_systematics_up.Draw(" SAME")
    hist_systematics_down.Draw(" SAME")


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
'''
if(systematic != ''):
    result_name = results_path(year,'stack',observable+'_'+systematic+'_groups.root')
else:
    result_name = results_path(year,'stack',observable+'_groups.root')
newfile = TFile(result_name, 'RECREATE')
canvas.Update()
canvas.Write()
newfile.Close()
canvas.SaveAs(result_name)
canvas.SaveAs("n_bjets_2017.png")

'''
print 'Data integral       : ', data_integral
print 'Signal integral     : ', signal_integral
print 'Background integral :  ', background_integral
print ''
for i in background_integral_i:
    print i[1],i[0]

if(systematic != ''):
    print "Systematic up : ", syst_up_integral
    print "Systematic down : ", syst_down_integral
print ''

exit = raw_input("Press key to quit : ") 
