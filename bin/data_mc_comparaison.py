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
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('stack_'+observable,'stack_'+observable, 800, 800)
canvas.UseCurrentStyle()

datafile_input = TFile('./results/'+year+'/flattree/'+observable+'_data.root')
rootfile_input = TFile('./results/'+year+'/flattree/'+observable+'_forComp.root')

################################################################################
## Create Histo 
################################################################################


r = 0.3
epsilon = 0.1

pad1 = TPad("pad1", "pad1", 0, r-epsilon, 1, 1)
pad1.SetBottomMargin(epsilon)
canvas.cd()
#pad1.SetLogy()
pad1.Draw()
pad1.cd()


###########
# data part
###########
hist_data = datafile_input.Get('data_obs')
data_integral = integral_complete(hist_data, max_bin)


# convenient variables
nbin    = hist_data.GetNbinsX()
min_bin = hist_data.GetXaxis().GetXmin()
max_bin = hist_data.GetXaxis().GetXmax()

###########rootfile_input
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
legend.AddEntry(hist_signal, "t#bar{t} SM", "f")
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

style_labels_counting(stack, 'Events', title)
stack.GetXaxis().SetLabelSize(0)
stack.GetXaxis().SetTitleSize(0)

if(year=='2016'):
    tdr.cmsPrel(35900., 13.,simOnly=False,thisIsPrelim=True)
elif(year=='2017'):
   tdr.cmsPrel(41500., 13.,simOnly=False,thisIsPrelim=True)

################################################################################
## Ratio
################################################################################

pad2 = TPad("pad2", "pad2", 0, 0, 1, r*(1-epsilon))
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0.4)
pad2.SetFillStyle(0)
canvas.cd()
pad2.Draw()
pad2.cd()

ratio_coef = 0.3

h_one = TH1F("one", "one", 1, min_bin, max_bin)
h_one.SetBinContent(1, 1)
h_one.SetLineWidth(1)
h_one.SetLineColor(15)
h_num = hist_data.Clone()
h_denom = hist_signal+hist_background
h_num.Divide(h_denom)
h_num.GetXaxis().SetTitle("aksjd")
ratio = THStack()
ratio.Add(h_num)

ratio.SetMaximum(1+ratio_coef)
ratio.SetMinimum(1-ratio_coef)
ratio.Draw()
h_one.Draw("SAME")


style_labels_counting(ratio, 'Ratio data/mc', title)
ratio.GetYaxis().SetLabelSize(0.1)
ratio.GetYaxis().SetTitleSize(0.1)
ratio.GetYaxis().SetTitleOffset(0.5)

ratio.GetXaxis().SetLabelSize(0.15)
ratio.GetXaxis().SetTitleSize(0.17)
ratio.GetXaxis().SetLabelOffset(0.01)


################################################################################
## Save
################################################################################

resultname = './results/'+year+'/comparaison/'+observable+'_'+year

rootfile_output = TFile(resultname+'.root', "RECREATE")
canvas.Write()
canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
rootfile_output.Close()


print 'Signal integral     : ','%.2f'%signal_integral
for i in background_integral_i:
    print i[1], '      : ', '%.2f'%i[0]
print 'Total Background integral : ', '%.2f'%background_integral
print 'Total MC integral : ', '%.2f'%(signal_integral+background_integral)
print 'Data integral       : ', data_integral
print 'Data/MC agreement  : ', '%.1f'%(100*(signal_integral+background_integral-data_integral)/(signal_integral+background_integral)), '%'

#raw_input('exit')
