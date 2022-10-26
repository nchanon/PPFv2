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
parser.add_argument('timebin', help='display the time bin')

args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title
timebin = int(args.timebin)

def integral_complete(histo, max_bin): 
     return histo.Integral()
#    return histo.Integral()+histo.GetBinContent(int(max_bin+1))+histo.GetBinContent(0)


nbin = 0
min_bin = 0
max_bin = 0
legend_coordinates = observable_values(observable)[1]
TH1.SetDefaultSumw2(1)
signal_integral = 0
background_integral_i = []
background_integral = 0
#data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
#canvas = TCanvas('compare_'+observable,'stack_'+observable, 800, 800)
#canvas.UseCurrentStyle()

#datafile_input = TFile('./results/'+year+'/flattree/'+observable+'_data.root')
stimebin="";
if (timebin==-1):
     stimebin = "_puold";
if (timebin==-2):
     stimebin = "_punew";
if (timebin==-3):
     stimebin = "_puinc";
if (timebin>=0):
     stimebin = "_put"+str(timebin);
rootfile_input = TFile('./results/'+year+'/flattree/'+observable+'_inclusive'+stimebin+'.root')

################################################################################
## Create Histo 
################################################################################


###########
# mc part
###########

maxi = []

# signal
hist_signal = rootfile_input.Get('signal')
a = hist_signal.Integral()
hist_signal.Scale(1/a)
maxi.append(hist_signal.GetMaximum())

# convenient variables
nbin    = hist_signal.GetNbinsX()
min_bin = hist_signal.GetXaxis().GetXmin()
max_bin = hist_signal.GetXaxis().GetXmax()

# backgrounds
hist_singletop = rootfile_input.Get('singletop')
a = hist_singletop.Integral()
hist_singletop.Scale(1/a)
maxi.append(hist_singletop.GetMaximum())

hist_ttx = rootfile_input.Get('ttx')
a = hist_ttx.Integral()
hist_ttx.Scale(1/a)
maxi.append(hist_ttx.GetMaximum())

hist_dibosons = rootfile_input.Get('dibosons')
a = hist_dibosons.Integral()
hist_dibosons.Scale(1/a)
maxi.append(hist_dibosons.GetMaximum())

hist_vjets = rootfile_input.Get('vjets')
a = hist_vjets.Integral()
hist_vjets.Scale(1/a)
maxi.append(hist_vjets.GetMaximum())


################################################################################
## Legend stuff
################################################################################

legend_args = (0.645, 0.65, 0.85, 0.92, '', 'NDC')
legend = TLegend(*legend_args)
legend.AddEntry(hist_signal, "t#bar{t} SM", "l")
legend.AddEntry(hist_singletop, "single top", "l")
legend.AddEntry(hist_vjets, "W/Z+jets", "l")
legend.AddEntry(hist_dibosons, "Dibosons", "l")
legend.AddEntry(hist_ttx, "t#bar{t}+X", "l")


################################################################################
## Draw stuff
################################################################################

canvas = TCanvas('compare_'+observable,'stack_'+observable, 800, 800)
canvas.UseCurrentStyle()

hist_signal.SetMinimum(0)
hist_signal.SetMaximum(max(maxi)*1.1)
#hist_signal.SetXTitle(title)
hist_signal.SetYTitle("a.u.")

style_labels_counting(hist_signal, "a.u.", title)

hist_signal.Draw("HIST")
hist_singletop.Draw("HISTsame")
hist_ttx.Draw("HISTsame")
hist_dibosons.Draw("HISTsame")
hist_vjets.Draw("HISTsame")

legend.Draw("SAME")

################################################################################
## Set Style
################################################################################

# line_color, line_width, fill_color, fill_style, marker_size, marker_style=1
style_histo(hist_signal, 2, 2, 0, 1004, 0)
style_histo(hist_singletop, 4, 2, 0, 1005, 0)
style_histo(hist_ttx, 8, 2, 0, 1005, 0)
style_histo(hist_dibosons, 42, 2, 0, 1005, 0)
style_histo(hist_vjets, 619, 2, 0, 1005, 0)

#style_labels_counting(stack, 'Events', title)
#stack.GetXaxis().SetLabelSize(0)
#stack.GetXaxis().SetTitleSize(0)

if(year=='2016'):
    tdr.cmsPrel(35900., 13.,simOnly=True,thisIsPrelim=True)
elif(year=='2017'):
   tdr.cmsPrel(41500., 13.,simOnly=True,thisIsPrelim=True)


################################################################################
## Save
################################################################################

resultname = './results/'+year+'/comparaison/mc_comp_'+observable+'_'+year

#rootfile_output = TFile(resultname+'.root', "RECREATE")
#canvas.Write()
canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
#rootfile_output.Close()


#print 'Signal integral     : ','%.2f'%signal_integral
#for i in background_integral_i:
#    print i[1], '      : ', '%.2f'%i[0]
#print 'Total Background integral : ', '%.2f'%background_integral
#print 'Total MC integral : ', '%.2f'%(signal_integral+background_integral)
#print 'Data integral       : ', data_integral
#print 'Data/MC agreement  : ', '%.1f'%(100*(signal_integral+background_integral-data_integral)/(signal_integral+background_integral)), '%'

#raw_input('exit')
