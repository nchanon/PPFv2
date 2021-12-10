#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import os
import argparse
import numpy as np


from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame
from ROOT import TGraphAsymmErrors
from ROOT import TStyle, gStyle, TColor, TLatex

import tools.tdrstyle as tdr
tdr.setTDRStyle()

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
#parser.add_argument('observable', help='display your observable')
#parser.add_argument('year', help='year of samples')
parser.add_argument('wilson',nargs='?', help='wilson coefficent (cLXX, cRXX, ...)', default='cLXX')

args = parser.parse_args()
#observable = args.observable
#year = args.year
wilson = args.wilson

legend_coordinates = [0.65, 0.75, 0.87, 0.87] 
TH1.SetDefaultSumw2(1)
signal_integral = 0
background_integral_i = []
background_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('sme modulations','sme modulations', 700, 700)
canvas.UseCurrentStyle()

gStyle.SetPalette(55);

################################################################################
## Create Histo 
################################################################################

cmunu = 0.01

sme_file = TFile('./results/2016/flattree/sme.root')
hist = sme_file.Get(wilson+'_details')
hist.Scale(cmunu)

hist_massbin = []
for i in range(8):
    if (i>0):
        hist_massbin.append(sme_file.Get(wilson+'_details'+str(i)))
	hist_massbin[-1].Scale(cmunu)

################################################################################
## Legend stuff
################################################################################

smassbin = [20,60,100,140,180,220,260]

legend = TLegend(0.66,0.64,0.93,0.94)
legend.SetTextSize(0.023)

legend.AddEntry(hist, 'inclusive (m_{ll}>'+str(smassbin[0])+' GeV)', 'l')
for i in range(6):
    legend.AddEntry(hist_massbin[i], str(smassbin[i])+'<m_{ll}<'+str(smassbin[i+1])+' GeV', 'l')
legend.AddEntry(hist_massbin[6], 'm_{ll}>'+str(smassbin[6])+' GeV', 'l')

################################################################################
## Draw stuff
################################################################################

min = hist.GetMinimum()
max = hist.GetMaximum()
hist.SetLineColor(1)
hist.SetLineWidth(2)
hist.Draw("hist")

for i in range(7):
    if (hist_massbin[i].GetMinimum() < min): min = hist_massbin[i].GetMinimum()
    if (hist_massbin[i].GetMaximum() > max): max = hist_massbin[i].GetMaximum()
    hist_massbin[i].SetLineColor(gStyle.GetColorPalette((i+1)*255/8))
    hist_massbin[i].SetLineWidth(2)
    hist_massbin[i].Draw("histsame")

legend.Draw("SAME")

################################################################################
## Set Style
################################################################################

is_center=True

hist.GetYaxis().SetRangeUser(min*1.2,max*1.2)
hist.GetYaxis().SetTitle("f(t) = #sigma_{SME} / #sigma_{SM} - 1")
hist.GetYaxis().SetTitleSize(0.04)
hist.GetYaxis().SetLabelSize(0.04)

hist.GetXaxis().SetRangeUser(0,24)
hist.GetXaxis().SetTitle('sidereal time (h)')
hist.GetXaxis().SetTitleSize(0.04)
hist.GetXaxis().SetLabelSize(0.04)

if(is_center):
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()

tdr.cmsPrel(-1,13.)

latex = TLatex()
latex.SetTextSize(0.65*gStyle.GetPadTopMargin())
latex.SetNDC()

if (wilson=="cLXX"): modwilson = "c_{L,XX}=-c_{L,YY}=" + str(cmunu)
if (wilson=="cLXY"): modwilson = "c_{L,XY}=c_{L,YX}=" + str(cmunu)
if (wilson=="cLXZ"): modwilson = "c_{L,XZ}=c_{L,ZX}=" + str(cmunu)
if (wilson=="cLYZ"): modwilson = "c_{L,YZ}=c_{L,ZY}=" + str(cmunu)

if (wilson=="cRXX"): modwilson = "c_{R,XX}=-c_{R,YY}=" + str(cmunu)
if (wilson=="cRXY"): modwilson = "c_{R,XY}=c_{R,YX}=" + str(cmunu)
if (wilson=="cRXZ"): modwilson = "c_{R,XZ}=c_{R,ZX}=" + str(cmunu)
if (wilson=="cRYZ"): modwilson = "c_{R,YZ}=c_{R,ZY}=" + str(cmunu)

if (wilson=="cXX"): modwilson = "c_{XX}=-c_{YY}=" + str(cmunu)
if (wilson=="cXY"): modwilson = "c_{XY}=c_{YX}=" + str(cmunu)
if (wilson=="cXZ"): modwilson = "c_{XZ}=c_{ZX}=" + str(cmunu)
if (wilson=="cYZ"): modwilson = "c_{YZ}=c_{ZY}=" + str(cmunu)

if (wilson=="dXX"): modwilson = "d_{XX}=-d_{YY}=" + str(cmunu)
if (wilson=="dXY"): modwilson = "d_{XY}=d_{YX}=" + str(cmunu)
if (wilson=="dXZ"): modwilson = "d_{XZ}=d_{ZX}=" + str(cmunu)
if (wilson=="dYZ"): modwilson = "d_{YZ}=d_{ZY}=" + str(cmunu)

latex.DrawLatex(0.25,0.9,modwilson)

#if(year=='2016'):
#    tdr.cmsPrel(35900.,13.)
#elif(year=='2017'):
#    tdr.cmsPrel(41530.,13.)

################################################################################
## Save
################################################################################

resultname = './results/2016/other/sme_time_'+wilson

#rootfile_output = TFile(resultname+'.root', "RECREATE")
#canvas.Write()
#canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
#rootfile_output.Close()

raw_input('exit')
