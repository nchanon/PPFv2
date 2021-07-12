#!/usr/bin/env python

import os
import sys
sys.path.append('./')


import argparse

from tools.sample_manager import *
from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

nbin = 24

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
observable = args.observable
year = args.year

TH1.SetDefaultSumw2(1)

################################################################################
## function
################################################################################

def max(a,b,c, nom):
    ratio = []
    for i in [a,b,c,nom]:
    #    if i>nom:
        ratio.append(i/nom)
        print ratio[-1]
    #    else:
    #        ratio.append(nom/i)
    max_rat = ratio[3]
    for i in range(3):
        if ratio[i]>max_rat:
            max_rat = ratio[i]
    return max_rat

def min(a,b,c, nom):
    ratio = []
    for i in [a,b,c, nom]:
    #    if i<nom:
        ratio.append(i/nom)
    #    else:
    #        ratio.append(nom/i)
    min_rat = ratio[3]
    for i in range(3):
        if ratio[i]<min_rat:
            min_rat = ratio[i]
    return min_rat

def rename(th1, name):
    th1.SetName(name)
    th1.SetTitle(name)


################################################################################
## Code body
################################################################################


mc_file = TFile('./results/'+year+'/flattree/'+observable+'.root')
mc_integral = 0
hist_mc = []

lumi_file = TFile('./inputs/timed/AllTimedSyst_'+year+'.root')
hist_weight = lumi_file.Get('h_SF_emu_sidereel_Full_UncBand')

h_nominal = mc_file.Get('signal')
nbin = h_nominal.GetNbinsX()
int_nom = h_nominal.Integral()
print int_nom
h_nominal.Draw()

####

alt_file = TFile('./results/'+year+'/flattree/'+observable+'_alt.root')
h_gluon = alt_file.Get('GluonMove')
h_erd = alt_file.Get('erdOn')
h_qcd = alt_file.Get('QCD')
h_hdampUp =  alt_file.Get('hdampUp')
h_hdampDown =  alt_file.Get('hdampDown')
h_CP5Up =  alt_file.Get('CP5Up')
h_CP5Down =  alt_file.Get('CP5Down')

################
## Color reco
################

h_gluon.SetLineColor(2)
h_erd.SetLineColor(5)
h_qcd.SetLineColor(3)

h_colorUp = h_nominal.Clone()
rename(h_colorUp,'signal_color_recoUp')

h_colorDown = h_nominal.Clone()
rename(h_colorDown, 'signal_color_recoDown')

for i in range(nbin):
    bin_max = max(h_gluon.GetBinContent(i+1),
                  h_erd.GetBinContent(i+1),
                  h_qcd.GetBinContent(i+1),
                  h_nominal.GetBinContent(i+1))
    bin_min = min(h_gluon.GetBinContent(i+1),
                  h_erd.GetBinContent(i+1),
                  h_qcd.GetBinContent(i+1),
                  h_nominal.GetBinContent(i+1))
    #h_colorUp.SetBinContent(i+1, h_nominal.GetBinContent(i+1)*bin_max)
    h_colorDown.SetBinContent(i+1, h_nominal.GetBinContent(i+1)*bin_min)
    h_colorUp.SetBinContent(i+1, h_nominal.GetBinContent(i+1)*(2-bin_min))

h_colorUp.SetLineColor(2)
h_colorDown.SetLineColor(3)

#h_colorUp.Draw('SAME')
#h_colorDown.Draw('SAME')


#int_gluon = int_nom/h_gluon.Integral()
#int_erd = int_nom/h_gluon.Integral()
#int_gluon = int_nom/h_gluon.Integral()


#h_gluon.Scale(int_gluon)
#h_erd.Scale(int_erd)
#h_qcd.Scale(int_gluon)

#h_gluon.Draw('SAME')
#h_erd.Draw('SAME')
#h_qcd.Draw('SAME')


################
## hdamp and CP5
################

h_CP5Up.SetLineColor(3)
h_CP5Down.SetLineColor(2)

h_hdampUp.SetLineColor(3)
h_hdampDown.SetLineColor(2)


int_cp5Up     = int_nom/h_CP5Up.Integral()
int_cp5Down   = int_nom/h_CP5Down.Intp/down variation (in tegral()
int_hdampUp   = int_nom/h_hdampUp.Integral()
int_hdampDown = int_nom/h_hdampDown.Integral()

h_hdampUp.Scale(int_hdampUp)
h_hdampDown.Scale(int_hdampDown)
h_CP5Up.Scale(int_cp5Up)
h_CP5Down.Scale(int_cp5Down)

rename(h_hdampUp,'signal_hdampUp')
rename(h_hdampDown,'signal_hdampDown')
rename(h_CP5Up,'signal_CP5Up')
rename(h_CP5Down,'signal_CP5Down')

h_CP5Up.Draw('SAME')
h_CP5Down.Draw('SAME')



################################################################################
## Code body
################################################################################

###########
# data part
###########


output = TFile('./results/'+year+'/flattree/'+observable+'_color_reco.root', "RECREATE")
h_colorUp.Write()
h_colorDown.Write()
h_hdampUp.Write()
h_hdampDown.Write()
h_CP5Up.Write()
h_CP5Down.Write()

output.Close()

raw_input('quit')