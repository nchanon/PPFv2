#!/usr/bin/env python

import os
import sys
sys.path.append('./')


import argparse

from tools.sample_manager import *
from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

#nbin = 24

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('timed', help='timed or inclusive', default='timed')
parser.add_argument('timebin', help='display the time bin')

args = parser.parse_args()
observable = args.observable
year = args.year
timed = args.timed
timebin = int(args.timebin)

doIntegrated = True

TH1.SetDefaultSumw2(1)

################################################################################
## functions
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

def get_hist_avg(hist, hist_nom):
    avg = hist.Integral() / hist_nom.Integral()
    #for i in range(nbin):
    #    avg += hist.GetBinContent(i+1)/hist_nom.GetBinContent(i+1)
    #avg = avg / nbin
    hist_avg = hist_nom.Clone() 
    hist_avg.Scale(avg)
    #for i in range(nbin): 
    #	hist_avg.SetBinContent(i+1,avg)
    return hist_avg

def rename(th1, name):
    th1.SetName(name)
    th1.SetTitle(name)



################################################################################
## Code body
################################################################################

if (timed=='timed'): stime=''
if (timed=='inclusive'): stime='_inclusive'

stimebin="";
if (timebin==-1):
     stimebin = "_puold";
if (timebin==-2):
     stimebin = "_punew";
if (timebin==-3):
     stimebin = "_puinc";
if (timebin>=0):
     stimebin = "_put"+str(timebin);

mc_file = TFile('./results/'+year+'/flattree/'+observable+stime+stimebin+'.root')

mc_integral = 0
hist_mc = []

#lumi_file = TFile('./inputs/timed/AllTimedSyst_'+year+'.root')
#hist_weight = lumi_file.Get('h_SF_emu_sidereel_Full_UncBand')

h_nominal = mc_file.Get('signal')
nbin = h_nominal.GetNbinsX()
int_nom = h_nominal.Integral()
print int_nom
h_nominal.Draw()

####

alt_file = TFile('./results/'+year+'/flattree/'+observable+'_alt'+stime+stimebin+'.root')

h_gluon = alt_file.Get('signal_GluonMove')
h_erd = alt_file.Get('signal_erdOn')
h_qcd = alt_file.Get('signal_QCD')
h_hdampUp =  alt_file.Get('signal_hdampUp')
h_hdampDown =  alt_file.Get('signal_hdampDown')
h_CP5Up =  alt_file.Get('signal_CP5Up')
h_CP5Down =  alt_file.Get('signal_CP5Down')
h_mtopUp = alt_file.Get('signal_mtopUp')
h_mtopDown = alt_file.Get('signal_mtopDown')

##################
## Color reco comb
##################

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
    h_colorUp.SetBinContent(i+1, h_nominal.GetBinContent(i+1)*bin_max)
    h_colorDown.SetBinContent(i+1, h_nominal.GetBinContent(i+1)*bin_min)
    #h_colorUp.SetBinContent(i+1, h_nominal.GetBinContent(i+1)*(2-bin_min))

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

if doIntegrated:
    h_hdampUpNew = get_hist_avg(h_hdampUp, h_nominal)
    h_hdampDownNew = get_hist_avg(h_hdampDown, h_nominal)
    h_CP5UpNew = get_hist_avg(h_CP5Up, h_nominal)
    h_CP5DownNew = get_hist_avg(h_CP5Down, h_nominal)
    rename(h_hdampUpNew,'signal_hdampUp')
    rename(h_hdampDownNew,'signal_hdampDown')
    rename(h_CP5UpNew,'signal_CP5Up')
    rename(h_CP5DownNew,'signal_CP5Down')
else:
    rename(h_hdampUp,'signal_hdampUp')
    rename(h_hdampDown,'signal_hdampDown')
    rename(h_CP5Up,'signal_CP5Up')
    rename(h_CP5Down,'signal_CP5Down')

#########################
## Color reco individuals
#########################

if doIntegrated:
    h_erdUpNew = get_hist_avg(h_erd, h_nominal)
    h_gluonUpNew = get_hist_avg(h_gluon, h_nominal)
    h_qcdUpNew = get_hist_avg(h_qcd, h_nominal)
    rename(h_erdUpNew, 'signal_erdOnUp')
    rename(h_gluonUpNew, 'signal_GluonMoveUp')
    rename(h_qcdUpNew, 'signal_QCDinspiredUp')
else:
    rename(h_erd, 'signal_erdOnUp')
    rename(h_gluon, 'signal_GluonMoveUp')
    rename(h_qcd, 'signal_QCDinspiredUp')

h_erdDown = h_nominal.Clone()
rename(h_erdDown, 'signal_erdOnDown')
h_gluonDown = h_nominal.Clone()
rename(h_gluonDown, 'signal_GluonMoveDown')
h_qcdDown = h_nominal.Clone()
rename(h_qcdDown, 'signal_QCDinspiredDown')


#h_CP5Up.Draw('SAME')
#h_CP5Down.Draw('SAME')


#########################
## Top mass 1 GeV
#########################

h_mtopUpNew = h_nominal.Clone()
h_mtopDownNew = h_nominal.Clone()

for i in range(nbin):
  val_up = h_nominal.GetBinContent(i+1) + (h_mtopUp.GetBinContent(i+1)-h_nominal.GetBinContent(i+1))/3.
  val_down = h_nominal.GetBinContent(i+1) + (h_mtopDown.GetBinContent(i+1)-h_nominal.GetBinContent(i+1))/3.
  h_mtopUpNew.SetBinContent(i+1,val_up)
  h_mtopDownNew.SetBinContent(i+1,val_down)

rename(h_mtopUpNew, 'signal_mtopUp')
rename(h_mtopDownNew, 'signal_mtopDown')


################################################################################
## Code body
################################################################################

##############
# storing part
##############


output = TFile('./results/'+year+'/flattree/'+observable+'_color_reco'+stime+stimebin+'.root', "RECREATE")
#h_nominal.Write()
if doIntegrated:
    h_hdampUpNew.Write()
    h_hdampDownNew.Write()
    h_CP5UpNew.Write()
    h_CP5DownNew.Write()
    h_erdUpNew.Write()
    h_gluonUpNew.Write()
    h_qcdUpNew.Write()
else:
    h_hdampUp.Write()
    h_hdampDown.Write()
    h_CP5Up.Write()
    h_CP5Down.Write()
    h_erd.Write()
    h_gluon.Write()
    h_qcd.Write()
 
h_erdDown.Write()
h_gluonDown.Write()
h_qcdDown.Write()
h_colorUp.Write()
h_colorDown.Write()
h_mtopUpNew.Write()
h_mtopDownNew.Write()

output.Close()

raw_input('quit')
