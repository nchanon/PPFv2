#!/usr/bin/env python

import os
import sys
sys.path.append('./')


import argparse

from tools.sample_manager import *
from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

#nbin = 24

doShapeOnly = True
#doShapeOnly = False

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('timebin', help='display the time bin')
args = parser.parse_args()
observable = args.observable
year = args.year
timebin = int(args.timebin)

TH1.SetDefaultSumw2(1)

################################################################################
## Functions
################################################################################

def applyNominalNorm(histo, hnom):
    area = histo.Integral()
    histo.Scale(hnom.Integral()/area)


################################################################################
## Code body
################################################################################

histograms = []

data_file = TFile('./results/'+year+'/flattree/'+observable+'_data.root')
for l in data_file.GetListOfKeys():
    hist = data_file.Get(l.GetName())
    histograms.append(hist)

h_nom = []
stimebin="";
if (timebin==-1):
     stimebin = "_puold";
if (timebin==-2):
     stimebin = "_punew";
if (timebin==-3):
     stimebin = "_puinc";
if (timebin>=0):
     stimebin = "_put"+str(timebin);

mc_file = TFile('./results/'+year+'/flattree/'+observable+'_inclusive'+stimebin+'.root')
for l in mc_file.GetListOfKeys():
    hist = mc_file.Get(l.GetName())
    histograms.append(hist)

    for proc in ttbar_list:
        if l.GetName()==proc:
            h_nom.append(mc_file.Get(l.GetName()))

    if doShapeOnly:
        curname = histograms[-1].GetName()
        for h in h_nom:
            proc = h.GetName()
            if TString(curname).Contains(proc):
		applyNominalNorm(histograms[-1], h)

    if TString(l.GetName()).Contains('syst_em_trig') or TString(l.GetName()).Contains('syst_b_uncorrelated') or TString(l.GetName()).Contains('syst_l_uncorrelated') or TString(l.GetName()).Contains('stat'):
        curname = histograms[-1].GetName()
        found = curname.find('Up')
        if (found==-1):
            found = curname.find('Down')
        newname = curname[:found] + '_' + year + curname[found:]
        histograms[-1].SetName(newname)

    if TString(l.GetName()).Contains('syst_qcdscale') or TString(l.GetName()).Contains('syst_ps_isr') or TString(l.GetName()).Contains('syst_ps_fsr'):
        curname = histograms[-1].GetName()
        found = curname.find('Up')
        if (found==-1):
            found = curname.find('Down')
        for h in h_nom:
            proc = h.GetName()
            if TString(curname).Contains(proc):
                newname = curname[:found] + '_' + proc + curname[found:]
                histograms[-1].SetName(newname)

    if TString(l.GetName()).Contains('syst_qcdscale') or TString(l.GetName()).Contains('syst_ps_isr') or TString(l.GetName()).Contains('syst_ps_fsr') or TString(l.GetName()).Contains('syst_pdfas') or TString(l.GetName()).Contains('syst_pt_top') or TString(l.GetName()).Contains('mtop'):
        curname = histograms[-1].GetName()
        for h in h_nom:
            proc = h.GetName()
            if TString(curname).Contains(proc):
                #area = histograms[-1].Integral()
                #histograms[-1].Scale(h.Integral()/area)
		if (doShapeOnly==False): 
		    applyNominalNorm(histograms[-1], h)



mc_alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco_inclusive'+stimebin+'.root')
for l in mc_alt_file.GetListOfKeys():
    hist = mc_alt_file.Get(l.GetName())
    histograms.append(hist)

    if TString(l.GetName()).Contains('mtop'):
        curname = histograms[-1].GetName()
        for h in h_nom:
            proc = h.GetName()
            if TString(curname).Contains(proc):
                #area = histograms[-1].Integral()
                #histograms[-1].Scale(h.Integral()/area)
                applyNominalNorm(histograms[-1], h)


mc_jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec_inclusive'+stimebin+'.root')
for l in mc_jec_file.GetListOfKeys():
    hist = mc_jec_file.Get(l.GetName())

    hname = hist.GetName()
    newname = hname
    if(hname.find('_up')!= -1):
        newname = hname[:-3]+'_jecUp'
    elif(hname.find('_down')!= -1):
        newname = hname[:-5]+'_jecDown'
    hist.SetName(newname)
    hist.SetTitle(newname)
    histograms.append(hist)

    if TString(hname).Contains("AbsoluteStat") or TString(hname).Contains("RelativeJEREC1") or TString(hname).Contains("RelativeJEREC2") or TString(hname).Contains("RelativePtEC1") or TString(hname).Contains("RelativePtEC2") or TString(hname).Contains("RelativeSample") or TString(hname).Contains("RelativeStatEC") or TString(hname).Contains("RelativeStatFSR") or TString(hname).Contains("RelativeStatHF") or TString(hname).Contains("TimePtEta"): #Careful for RelativeSample
        found = curname.find('_jecUp')
        if (found==-1):
            found = curname.find('_jecDown')
        newname = curname[:found] + '_' + year + curname[found:]
	histograms[-1].SetName(newname)
        histograms[-1].SetTitle(newname)

    if doShapeOnly:
        curname = histograms[-1].GetName()
        for h in h_nom:
            proc = h.GetName()
            if TString(curname).Contains(proc):
                applyNominalNorm(histograms[-1], h)


out = './combine/'+year+'/inclusive/inputs/'
output = TFile(out+observable+'_inclusive'+stimebin+'.root', "RECREATE")
for h in histograms:
    h.Write()
output.Close()

print 'Produced '+output.GetName()
