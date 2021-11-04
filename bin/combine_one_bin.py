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

def rename(th1, name):
    th1.SetName(name)
    th1.SetTitle(name)


################################################################################
## Code body
################################################################################

histograms = []
data_number_of_event = []


# Timed histo
lumi_syst_up = {}
lumi_syst_down = {}

lumi_file = TFile('./inputs/timed/AllTimedSyst_'+year+'.root')
lumi_histo = lumi_file.Get('hInstLumi_DataScaleFactor')
hist_weight = lumi_file.Get('h_SF_emu_sidereel_Full_UncBand')

for l in lumi_file.GetListOfKeys():
    if not TString(l.GetName()).Contains('DataScaleFactor'):
        if TString(l.GetName()).Contains('_Up'):
            if TString(l.GetName()).Contains('Inclusive'):
                lumi_syst_up.update({systematic_time_list[0]: lumi_file.Get(l.GetName())})
            elif TString(l.GetName()).Contains('Stability'):
                lumi_syst_up.update({systematic_time_list[1]: lumi_file.Get(l.GetName())})
            elif TString(l.GetName()).Contains('Linearity'):
                lumi_syst_up.update({systematic_time_list[2]: lumi_file.Get(l.GetName())})
        elif TString(l.GetName()).Contains('_Down'):
            if TString(l.GetName()).Contains('Inclusive'):
                lumi_syst_down.update({systematic_time_list[0]: lumi_file.Get(l.GetName())})
            elif TString(l.GetName()).Contains('Stability'):
                lumi_syst_down.update({systematic_time_list[1]: lumi_file.Get(l.GetName())})
            elif TString(l.GetName()).Contains('Linearity'):
                lumi_syst_down.update({systematic_time_list[2]: lumi_file.Get(l.GetName())})
    if TString(l.GetName()).Contains('SF_emu_sidereel_Full_UncBand'):
        lumi_syst_up.update({systematic_time_list[3]: lumi_file.Get(l.GetName())})
# Body
data_file = TFile('./results/'+year+'/flattree/'+observable+'.root')
for l in data_file.GetListOfKeys():
    if not TString(l.GetName()).Contains('data_obs'):
        hist = data_file.Get(l.GetName())
        hist.Scale(1./nbin)
        histograms.append(hist)


#alt_file = TFile('./results/'+year+'/flattree/'+observable+'_alt.root')
alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco.root')
for l in alt_file.GetListOfKeys():
    h = alt_file.Get(l.GetName())
    h.Scale(1./nbin)
    histograms.append(h)

jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec.root')
for l in jec_file.GetListOfKeys():
    hj = jec_file.Get(l.GetName())
    hname = l.GetName() 
    if(hname.find('TotalUp')!= -1):
        name = hname[:-7]+'jecUp'
    elif(hname.find('TotalDown')!= -1):
        name = hname[:-9]+'jecDown'
    rename(hj,name)
    hj.Scale(1./nbin)
    histograms.append(hj)

#print histograms[27].GetName(), histograms[27].GetTitle(), len(histograms)


out = './combine/'+year+'/one_bin/inputs/'
for n in range(nbin):


    
    hist_tim = []

    for g in ttbar_list:
        for s in systematic_time_list:
            hist_up = data_file.Get(g).Clone()
            hist_down = data_file.Get(g).Clone()
            if s == 'emu_trig':
                hist_up.Scale(1+lumi_syst_up[s].GetBinError(n+1))
                hist_up.SetName(g+'_'+s+'Up')
                hist_up.SetTitle(g+'_'+s+'Up')
                hist_down.Scale(1-lumi_syst_up[s].GetBinError(n+1))
                hist_down.SetName(g+'_'+s+'Down')
                hist_down.SetTitle(g+'_'+s+'Down')
            else:
                hist_up.Scale(lumi_syst_up[s].GetBinContent(n+1))
                hist_up.SetName(g+'_'+s+'Up')
                hist_up.SetTitle(g+'_'+s+'Up')
                hist_down.Scale(lumi_syst_down[s].GetBinContent(n+1))
                hist_down.SetName(g+'_'+s+'Down')
                hist_down.SetTitle(g+'_'+s+'Down')
            hist_up.Scale(hist_weight.GetBinContent(n+1))
            hist_down.Scale(hist_weight.GetBinContent(n+1))
            hist_tim.append(hist_up)
            hist_tim.append(hist_down)
    '''
        if n == 0:
            print g, s
            print histograms[28].GetName(), histograms[28].GetTitle(), len(histograms)
    '''

    time_file = TFile('./results/'+year+'/flattree/'+observable+'_data_timed'+str(nbin)+'.root')
    histograms_timmed = time_file.Get('data_obs_bin'+str(n))
    histograms_timmed.SetName('data_obs')
    histograms_timmed.SetTitle('data_obs')
    #histograms_timmed.Scale(lumi_histo.GetBinContent(n+1)*hist_weight.GetBinContent(n+1))
    histograms_timmed.Scale(hist_weight.GetBinContent(n+1))
    data_number_of_event.append(histograms_timmed.Integral())
    print 'Integral bin '+str(n)+' : '+str(data_number_of_event[n])

    output = TFile(out+observable+'_'+str(nbin)+'_'+str(n)+'.root', "RECREATE")
    histograms_timmed.Write()
    for h in histograms:
        h.Scale(hist_weight.GetBinContent(n+1))
        h.Write()
    for h in hist_tim:
        h.Write()
    output.Close()
    del hist_tim


file_txt = ''
for i in data_number_of_event:
    file_txt += str(i)+'\n'

file = open('./combine/'+year+'/'+observable+'_noe_data_timed.txt','w') 
file.write(file_txt) 
file.close() 
