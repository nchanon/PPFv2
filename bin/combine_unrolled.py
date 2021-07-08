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
parser.add_argument('wilson',nargs='?', help='wilson coefficent (cLXX, cCXX, ...)', default='cLXX')

args = parser.parse_args()
observable = args.observable
year = args.year
wilson = args.wilson

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

###########
# data part
###########

data_file = TFile('./results/'+year+'/flattree/'+observable+'_data_timed24.root')
data_integral = 0
hist_data_in = []

for i in range(nbin): 
    hist_data_in.append(data_file.Get('data_obs_bin'+str(i)))
    data_integral += hist_data_in[i].Integral()

# convenient variables
binning = hist_data_in[0].GetNbinsX()
min_bin = hist_data_in[0].GetXaxis().GetXmin()
max_bin = hist_data_in[0].GetXaxis().GetXmax()

canvas_data = TCanvas('data_obs', 'data_obs', 1000, 800)
hist_data = TH1F('data_obs','data_obs', binning*nbin,  0, nbin)
for iobs in range(nbin):
    for n in range(binning):
        hist_data.SetBinContent(n + iobs*binning + 1 , hist_data_in[iobs].GetBinContent(n+1))
        hist_data.SetBinError(n + iobs*binning + 1 , hist_data_in[iobs].GetBinError(n+1))

###########
# MC part
###########

mc_file = TFile('./results/'+year+'/flattree/'+observable+'.root')
mc_integral = 0
hist_mc = []

lumi_file = TFile('./inputs/timed/AllTimedSyst_'+year+'.root')
hist_weight = lumi_file.Get('h_SF_emu_sidereel_Full_UncBand')

index = 0
for l in mc_file.GetListOfKeys():
    if not TString(l.GetName()).Contains('data_obs'):
        hist = mc_file.Get(l.GetName())
        hist_mc.append(TH1F("","", binning*nbin,  0, nbin))
        for iobs in range(nbin):
            for j in range(binning):
                hist_mc[index].SetBinContent(j + iobs*binning + 1 , hist.GetBinContent(j+1)*hist_weight.GetBinContent(i+1))
                hist_mc[index].SetBinError(j + iobs*binning + 1 , hist.GetBinError(i+1))
        del hist
        hist_mc[index].SetName(l.GetName())
        hist_mc[index].SetTitle(l.GetName())
        hist_mc[index].Scale(1./nbin)
        for g in ttbar_list:
            if TString(l.GetName()) == g:
                mc_integral += hist_mc[index].Integral()
        index += 1


## timed systematics

lumi_syst_up = {}
lumi_syst_down = {}
#lumi_file = TFile('./inputs/timed/AllTimedSyst_'+year+'.root')

for l in lumi_file.GetListOfKeys():
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

for l in systematic_time_list:
    for g in hist_mc:
        for b in ttbar_list:
            if g.GetName() == b:
                hist_up = g.Clone()
                hist_up.SetName(b+'_'+l+'Up')
                hist_up.SetTitle(b+'_'+l+'Up')
                hist_down = g.Clone()
                hist_down.SetName(b+'_'+l+'Down')
                hist_down.SetTitle(b+'_'+l+'Down')
                for iobs in range(nbin):
                    for j in range(binning):
                        if l == 'emu_trig':
                            hist_up.SetBinContent(j + iobs*binning + 1 , 
                                                hist_up.GetBinContent(j+1)*(1+lumi_syst_up[l].GetBinError(iobs+1)))
                            hist_up.SetBinError(j + iobs*binning + 1 , 
                                                hist_up.GetBinError(j+1)*(lumi_syst_up[l].GetBinError(iobs+1)))
                            hist_down.SetBinContent(j + iobs*binning + 1 , 
                                                hist_down.GetBinContent(j+1)*(1-lumi_syst_up[l].GetBinError(iobs+1)))
                            hist_down.SetBinError(j + iobs*binning + 1 , 
                                                hist_down.GetBinError(j+1)*(lumi_syst_up[l].GetBinError(iobs+1)))
                        else:
                            hist_up.SetBinContent(j + iobs*binning + 1 , 
                                                hist_up.GetBinContent(j+1)*lumi_syst_up[l].GetBinContent(iobs+1))
                            hist_up.SetBinError(j + iobs*binning + 1 , 
                                                hist_up.GetBinError(j+1)*lumi_syst_up[l].GetBinError(iobs+1))
                            hist_down.SetBinContent(j + iobs*binning + 1 , 
                                                hist_down.GetBinContent(j+1)*lumi_syst_down[l].GetBinContent(iobs+1))
                            hist_down.SetBinError(j + iobs*binning + 1 , 
                                                hist_down.GetBinError(j+1)*lumi_syst_down[l].GetBinError(iobs+1))
                hist_mc.append(hist_up)
                hist_mc.append(hist_down)

###########
# SME part
###########

cmunu = 0.001
sme_file = TFile('./results/'+year+'/flattree/sme.root')
sme_sig = sme_file.Get(wilson)

hist_sme = []

for g in hist_mc:
    if g.GetName().find('signal') != -1:
        hist_sme.append(g.Clone())
        name = wilson+hist_sme[-1].GetName()[6:]
        hist_sme[-1].SetName(name)
        hist_sme[-1].SetTitle(name)

for n in range(len(hist_sme)):
    for i in range(nbin):
        for j in range(binning):
            hist_sme[n].SetBinContent(j + i*binning + 1, 
                                        hist_sme[n].GetBinContent(j + i*binning + 1 )*(1+cmunu*sme_sig.GetBinContent(i+1))
                                    )
            hist_sme[n].SetBinError(j + i*binning + 1, 
                                        hist_sme[n].GetBinError(j + i*binning + 1 )*(1+cmunu*sme_sig.GetBinError(i+1))
                                    )


hist_alt_jec = []

jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec.root')
for l in jec_file.GetListOfKeys():
    hh = jec_file.Get(l.GetName())
    hh.Scale(1./nbin)
    hname = l.GetName() 
    if(hname.find('TotalUp')!= -1):
        name = hname[:-7]+'jecUp'
    elif(hname.find('TotalDown')!= -1):
        name = hname[:-9]+'jecDown'
    h_jec = TH1F(name,name, binning*nbin,  0, nbin)
    for i in range(nbin):
        for j in range(binning):
            h_jec.SetBinContent(j + i*binning + 1, hh.GetBinContent(j+1)*hist_weight.GetBinContent(i+1)
                                    )
            h_jec.SetBinError(j + i*binning + 1, hh.GetBinError(j + 1))
    hist_alt_jec.append(h_jec)


alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco.root')
for l in alt_file.GetListOfKeys():
    hh = alt_file.Get(l.GetName())
    hh.Scale(1./nbin)
    name = l.GetName()
    h_alt = TH1F(name, name, binning*nbin,  0, nbin)
    for i in range(nbin):
        for j in range(binning):
            h_alt.SetBinContent(j + i*binning + 1, hh.GetBinContent(j+1)*hist_weight.GetBinContent(i+1)
                                    )
            h_alt.SetBinError(j + i*binning + 1, hh.GetBinError(j + 1))
    hist_alt_jec.append(h_alt)


print 'data = '+str(data_integral)
print 'mc   = '+str(mc_integral)

out = './combine/'+year+'/unrolled/inputs/'
output = TFile(out+observable+'.root', "RECREATE")
hist_data.Write()
for l in hist_sme:
    l.Write()
for l in hist_mc:
    l.Write()
for l in hist_alt_jec:
    l.Write()
output.Close()

#for i in range(l.GetNbinsX()-400):
#    print str(i)+'#########'+str(hist_mc[0].GetBinContent(i+1))

cmd = 'cp '+out+observable+'.root '+'./combine/'+year+'/sme/inputs/'+observable+'_'+wilson+'.root'
os.system(cmd)

file = open('./combine/'+year+'/'+observable+'_noe_data.txt','w') 
file.write(str(data_integral)) 
file.close() 