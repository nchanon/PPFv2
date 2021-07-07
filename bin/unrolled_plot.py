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




rootfile_mc = TFile('combine/2017/sme/inputs/m_dilep_cLXX.root')


################################################################################
## Create Histo 
################################################################################


###########
# data part
###########
hist  = rootfile_mc.Get('signal')

#hist_mc.Scale(1./hist_mc.Integral())
#hist_data.Scale(1./hist_data.Integral())
################################################################################
## Legend stuff
################################################################################

legend_args = (0.645, 0.79, 0.985, 0.85, '', 'NDC')
legend = TLegend(*legend_args)
legend.AddEntry(hist, "MC", "l")
legend_box(legend, legend_coordinates)

################################################################################
## Draw stuff
################################################################################

hist.Draw(" E HIST")
legend.Draw("SAME")

################################################################################
## Set Style
################################################################################

is_center=True

hist.GetYaxis().SetTitle('Events')
hist.GetXaxis().SetTitle('sidereal time')
hist.GetYaxis().SetMaxDigits(4)
if(is_center):
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()

# line_color, line_width, fill_style, marker_style
#style_histo(hist_data, 2, 5, 1, 3001, 1, 1)
#style_histo(hist_mc,   4, 5, 1, 3001, 1, 1)
style_histo(hist,   4, 1, 4, 3005, 0)

tdr.cmsPrel(41530., 13.)

################################################################################
## Save
################################################################################

resultname = './results/2017/other/unrolled'
'''
rootfile_output = TFile(resultname+'.root', "RECREATE")
canvas.Write()
canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
rootfile_output.Close()
'''
raw_input('exit')