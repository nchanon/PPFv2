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
parser.add_argument('observable', help='display your observable')
#parser.add_argument('year', help='year of samples')
parser.add_argument('wilson', help='wilson coefficent (cL, cR, ...)', default='cL')
parser.add_argument('singletop', help='year of samples', default='')

args = parser.parse_args()
observable = args.observable
#year = args.year
wilson = args.wilson
singletop = args.singletop

print observable+' '+wilson+' '+singletop

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

prefix=''
suffix=''

if singletop=='singletop':
    prefix='singletop_'
    suffix='_singletop'

cmunu = 0.01

sme_file = TFile('./results/2016/flattree/'+observable+'_sme.root')

skinbin="2"

hist_XX = sme_file.Get(prefix+wilson+'XX_details'+skinbin)
hist_XX.Scale(cmunu)
hist_XY = sme_file.Get(prefix+wilson+'XY_details'+skinbin)
hist_XY.Scale(cmunu)
hist_XZ = sme_file.Get(prefix+wilson+'XZ_details'+skinbin)
hist_XZ.Scale(cmunu)
hist_YZ = sme_file.Get(prefix+wilson+'YZ_details'+skinbin)
hist_YZ.Scale(cmunu)

'''
smassbin = [20,60,100,140,180,220,260]
sptemubin = [0,35,70,105,140,175,210,245]

if observable=='m_dilep':
    svar = 'm_{ll}'
    skinbin = smassbin
elif observable=='pt_emu':
    svar = 'p_{T,ll}'
    skinbin = sptemubin

hist_massbin = []
for i in range(8):
    if (observable=='m_dilep' and i>0):
        hist_massbin.append(sme_file.Get(prefix+wilson+'_details'+str(i)))
	hist_massbin[-1].Scale(cmunu)
    if (observable=='pt_emu' and i<7):
        hist_massbin.append(sme_file.Get(prefix+wilson+'_details'+str(i)))
        hist_massbin[-1].Scale(cmunu)
'''
################################################################################
## Legend stuff
################################################################################

if (wilson=="cL"): 
    modwilson1 = "c_{L,XX}=-c_{L,YY}=" + str(cmunu)
    modwilson2 = "c_{L,XY}=c_{L,YX}=" + str(cmunu)
    modwilson3 = "c_{L,XZ}=c_{L,ZX}=" + str(cmunu)
    modwilson4 = "c_{L,YZ}=c_{L,ZY}=" + str(cmunu)

if (wilson=="cR"): 
    modwilson1 = "c_{R,XX}=-c_{R,YY}=" + str(cmunu)
    modwilson2 = "c_{R,XY}=c_{R,YX}=" + str(cmunu)
    modwilson3= "c_{R,XZ}=c_{R,ZX}=" + str(cmunu)
    modwilson4 = "c_{R,YZ}=c_{R,ZY}=" + str(cmunu)

if (wilson=="c"): 
    modwilson1 = "c_{XX}=-c_{YY}=" + str(cmunu)
    modwilson2 = "c_{XY}=c_{YX}=" + str(cmunu)
    modwilson3 = "c_{XZ}=c_{ZX}=" + str(cmunu)
    modwilson4 = "c_{YZ}=c_{ZY}=" + str(cmunu)

if (wilson=="d"): 
    modwilson1 = "d_{XX}=-d_{YY}=" + str(cmunu)
    modwilson2 = "d_{XY}=d_{YX}=" + str(cmunu)
    modwilson3 = "d_{XZ}=d_{ZX}=" + str(cmunu)
    modwilson4 = "d_{YZ}=d_{ZY}=" + str(cmunu)

#smassbin = [20,60,100,140,180,220,260]
#sptemubin = [0,35,70,105,140,175,210,245]

legend = TLegend(0.66,0.8,0.93,0.94)
legend.SetTextSize(0.023)

legend.AddEntry(hist_XX, modwilson1, 'l')
legend.AddEntry(hist_XY, modwilson2, 'l')
legend.AddEntry(hist_XZ, modwilson3, 'l')
legend.AddEntry(hist_YZ, modwilson4, 'l')

#for i in range(6):
#    legend.AddEntry(hist_massbin[i], str(skinbin[i])+'<'+svar+'<'+str(skinbin[i+1])+' GeV', 'l')
#legend.AddEntry(hist_massbin[6], svar+'>'+str(skinbin[6])+' GeV', 'l')


################################################################################
## Draw stuff
################################################################################

min = hist_XX.GetMinimum()
max = hist_XX.GetMaximum()
hist_XX.SetLineColor(1)
hist_XX.SetLineWidth(2)
hist_XX.Draw("hist")
hist_XY.SetLineColor(2)
hist_XY.SetLineWidth(2)
hist_XY.Draw("histSAME")
hist_XZ.SetLineColor(4)
hist_XZ.SetLineWidth(2)
hist_XZ.Draw("histSAME")
hist_YZ.SetLineColor(8)
hist_YZ.SetLineWidth(2)
hist_YZ.Draw("histSAME")

'''
for i in range(7):
    if (hist_massbin[i].GetMinimum() < min): min = hist_massbin[i].GetMinimum()
    if (hist_massbin[i].GetMaximum() > max): max = hist_massbin[i].GetMaximum()
    hist_massbin[i].SetLineColor(gStyle.GetColorPalette((i+1)*255/8))
    hist_massbin[i].SetLineWidth(2)
    hist_massbin[i].Draw("histsame")
'''
legend.Draw("SAME")

################################################################################
## Set Style
################################################################################

is_center=True

hist_XX.GetYaxis().SetRangeUser(min*1.2,max*1.2)
hist_XX.GetYaxis().SetTitle("f(t) = #sigma_{SME} / #sigma_{SM} - 1")
hist_XX.GetYaxis().SetTitleSize(0.04)
hist_XX.GetYaxis().SetLabelSize(0.04)

hist_XX.GetXaxis().SetRangeUser(0,24)
hist_XX.GetXaxis().SetTitle('sidereal time (h)')
hist_XX.GetXaxis().SetTitleSize(0.04)
hist_XX.GetXaxis().SetLabelSize(0.04)

if(is_center):
    hist_XX.GetXaxis().CenterTitle()
    hist_XX.GetYaxis().CenterTitle()

tdr.cmsPrel(-1,13.)

latex = TLatex()
latex.SetTextSize(0.65*gStyle.GetPadTopMargin())
latex.SetNDC()
if skinbin!="":
    latex.DrawLatex(0.2,0.9,"Particle level, n_{bjets}="+skinbin);

#if(year=='2016'):
#    tdr.cmsPrel(35900.,13.)
#elif(year=='2017'):
#    tdr.cmsPrel(41530.,13.)

################################################################################
## Save
################################################################################

resultname = './results/2016/other/'+observable+'_sme_time_'+wilson+suffix+'_kinindep'
if skinbin!="":
    resultname += "_"+skinbin

#rootfile_output = TFile(resultname+'.root', "RECREATE")
#canvas.Write()
#canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
#rootfile_output.Close()


raw_input('exit')
