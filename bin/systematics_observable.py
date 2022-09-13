#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame
from ROOT import gROOT

gROOT.SetBatch()


import tools.tdrstyle as tdr
tdr.setTDRStyle()

#print "Systematics list"
#print ['Total'] #systematic_list

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('timed', help='timed or inclusive')
parser.add_argument('systematic', help='display your systematic')
parser.add_argument('title', help='display your observable title')
parser.add_argument('timebin', help='display the time bin')

args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title
systematic = args.systematic
timed = args.timed
timebin = int(args.timebin)
stimed = ''
if (timed=="inclusive"):
    stimed = '_inclusive'

#nbin = 0
#min_bin = 0
#max_bin = 0
legend_coordinates = observable_values(observable)[1]
TH1.SetDefaultSumw2(1)
signal_integral = 0
mc_integral_i = []
mc_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('stack_'+observable,'stack_'+observable, 800, 800)
canvas.UseCurrentStyle()


#data_input = TFile('./results/'+year+'/flattree/'+observable+'_data.root')
stimebin="";
if (timebin==-1):
     stimebin = "_puold";
if (timebin==-2):
     stimebin = "_punew";
if (timebin==-3):
     stimebin = "_puinc";
if (timebin>=0):
     stimebin = "_put"+str(timebin);
rootfile_input = TFile('./results/'+year+'/flattree/'+observable+stimed+stimebin+'.root')
rootfile_input_jec = TFile('./results/'+year+'/flattree/'+observable+'_jec'+stimed+stimebin+'.root')
rootfile_input_alt = TFile('./results/'+year+'/flattree/'+observable+'_alt'+stimed+stimebin+'.root')
rootfile_input_color_reco = TFile('./results/'+year+'/flattree/'+observable+'_color_reco'+stimed+stimebin+'.root')

rootfile_input_time = []
if timed=='timed':
    for i in range(24):
        rootfile_input_time.append(TFile('./combine/'+year+'/one_bin/inputs/'+observable+'_24_'+str(i)+'.root'))

################################################################################
## Create Histo 
################################################################################


###########
# data part
###########
#hist_data = data_input.Get('data_obs')
#data_integral = hist_data.Integral()


# convenient variables
#nbin    = hist_data.GetNbinsX()
#min_bin = hist_data.GetXaxis().GetXmin()
#max_bin = hist_data.GetXaxis().GetXmax()


###########
# mc part
###########

hist_mc = []
hist_mc_up = []
hist_mc_down = []

#print ttbar_list
rootfile_input_syst = rootfile_input
slist = ttbar_list

if (systematic == 'syst_pt_top' or systematic == 'CP5' or systematic == 'hdamp' or systematic == 'mtop' or systematic == 'erd' or systematic == 'QCD' or systematic == 'GluonMove'):
   slist = ['signal']

for l in rootfile_input.GetListOfKeys():
    for s in slist:
        #print l.GetName(),s
        if(l.GetName() == s):
            #print 'Nominal'
            hist_mc.append(rootfile_input.Get(l.GetName()))
            mc_integral_i.append([rootfile_input.Get(l.GetName()).Integral(), l.GetName()])

i=0
for jec_syst in jec_list[year]:
    if i%2==0 and systematic==jec_syst[:-3]:
        rootfile_input_syst = rootfile_input_jec
        slist = ttbar_list
    i=i+1

#if (systematic == 'Total' or systematic == 'Absolute' or systematic == 'Absolute_'+year or systematic == 'FlavorQCD' or systematic == 'BBEC1' or systematic == 'BBEC1_'+year or systematic=='RelativeBal' or systematic=='RelativeSample_'+year): #Total JEC
#    rootfile_input_syst = rootfile_input_jec
#    slist = ttbar_list

if (systematic == 'CP5' or systematic == 'hdamp' or systematic == 'mtop' or systematic == 'erd' or systematic == 'QCD' or systematic == 'GluonMove'):
    rootfile_input_syst = rootfile_input_alt
if systematic == 'mtop':
    rootfile_input_syst = rootfile_input_color_reco


for l in rootfile_input_syst.GetListOfKeys():
    for s in slist: #ttbar_list:
        #print l.GetName(),s
#        if(l.GetName() == s):
#	    print 'Nominal'
#            hist_mc.append(rootfile_input.Get(l.GetName()))
#            mc_integral_i.append([rootfile_input.Get(l.GetName()).Integral(), l.GetName()])
        if(TString(l.GetName()).Contains(s) and TString(l.GetName()).Contains(systematic) and (TString(l.GetName()).Contains('Up') or TString(l.GetName()).Contains('up'))):
	    print 'Up'
            name = s + '_' + systematic + 'Up'
            hist_mc_up.append(rootfile_input_syst.Get(l.GetName()))
        elif(TString(l.GetName()).Contains(s) and TString(l.GetName()).Contains(systematic) and (TString(l.GetName()).Contains('Down') or TString(l.GetName()).Contains('down'))):
	    print 'Down'
            name = s + '_' + systematic + 'Down'
            hist_mc_down.append(rootfile_input_syst.Get(l.GetName()))
	elif(TString(l.GetName()).Contains(s) and TString(l.GetName()).Contains(systematic)):
	    print 'Up'
	    name = s + '_' + systematic + 'Up'
            hist_mc_up.append(rootfile_input_syst.Get(l.GetName()))
	    hist_mc_down.append(hist_mc[-1].Clone())

obsname = observable

if timed=='timed' and (TString(systematic).Contains('emu_trig') or TString(systematic).Contains('lumi_stability') or TString(systematic).Contains('lumi_linearity') or TString(systematic).Contains('syst_pu')):
    del hist_mc[:]
    del hist_mc_up[:]
    del hist_mc_down[:]
    obsname = "sidereal"
    slist = ttbar_list
    for s in slist:
	h_time = TH1F(s, s, 24,0,24)
        h_time_up = TH1F(s + '_' + systematic + 'Up',s + '_' + systematic + 'Up',24,0,24)
        h_time_down = TH1F(s + '_' + systematic + 'Down',s + '_' + systematic + 'Down',24,0,24)
        for i in range(24):
	    for l in rootfile_input_time[i].GetListOfKeys():
		if l.GetName()==s:
		    print(str(i)+' Nom')
		    h_time.SetBinContent(i+1, rootfile_input_time[i].Get(l.GetName()).Integral())
		if(TString(l.GetName()).Contains(s) and TString(l.GetName()).Contains(systematic) and (TString(l.GetName()).Contains('Up'))):
		    print(str(i)+' Up')
		    h_time_up.SetBinContent(i+1, rootfile_input_time[i].Get(l.GetName()).Integral())
                    #name = s + '_' + systematic + 'Up'
		    #hist_mc_up.append(rootfile_input_time[i].Get(l.GetName()))
                if(TString(l.GetName()).Contains(s) and TString(l.GetName()).Contains(systematic) and (TString(l.GetName()).Contains('Down'))):
                    print(str(i)+' Down')
		    h_time_down.SetBinContent(i+1, rootfile_input_time[i].Get(l.GetName()).Integral())
                    #name = s + '_' + systematic + 'Down'
                    #hist_mc_down.append(rootfile_input_time[i].Get(l.GetName()))
        hist_mc.append(h_time)
	hist_mc_up.append(h_time_up)
	hist_mc_down.append(h_time_down)


print "Selected histos"
for l in range(len(hist_mc)):
    print hist_mc[l].GetName(),  hist_mc_up[l].GetName(),  hist_mc_down[l].GetName() 



################################################################################
## Legend stuff
################################################################################

legend_args = (0.2, 0.79, 0.55, 0.95, '', 'NDC')
legend = []
for index in range(len(hist_mc)):
    legend.append(TLegend(*legend_args))
    legend[index].SetHeader(systematic + ' ' + year)
    legend[index].AddEntry(hist_mc[index], hist_mc[index].GetName())
    legend[index].AddEntry(hist_mc_up[index], 'up')    
    legend[index].AddEntry(hist_mc_down[index], 'down')

    #legend_box(legend[index], legend_coordinates)

################################################################################
## Set Style
################################################################################

for index in range(len(hist_mc)):
    #style_histo(hist_mc_up[index], 7, 2, 0, 1000, 0)
    #style_histo(hist_mc_down[index], 6, 2, 0, 1000, 0)
    style_histo(hist_mc_up[index], 2, 1, 2, 3004, 0)
    style_histo(hist_mc_down[index], 4, 1, 4, 3005, 0)
    style_histo(hist_mc[index], 1, 1, 0, 3001, 1, 20)
    if (TString(systematic).Contains('emu_trig') or TString(systematic).Contains('lumi_stability') or TString(systematic).Contains('lumi_linearity')):
	style_labels_counting(hist_mc[index], 'Uncertainty (in %)', 'sidereal time (h)')
    else:
        style_labels_counting(hist_mc[index], 'Uncertainty (in %)', title)


if(year=='2016'):
    tdr.cmsPrel(35900., 13.)
elif(year=='2017'):
    tdr.cmsPrel(41530., 13.)

################################################################################
## Save
################################################################################
import numpy as np

edge = []

for h in range(len(hist_mc)):
    for i in range(hist_mc[h].GetNbinsX()):
        val = float(hist_mc[h].GetBinContent(i+1))
	val_err = float(hist_mc[h].GetBinError(i+1))
        hist_mc[h].SetBinContent(i+1,0)
        #hist_mc[h].SetBinError(i+1,hist_mc[h].GetBinError(i+1)/np.sqrt(val))
	if val!=0:
	    hist_mc[h].SetBinError(i+1,val_err/val*100)
	elif val==0:
	     hist_mc[h].SetBinError(i+1,0)
        up   = hist_mc_up[h].GetBinContent(i+1)-val
        down = hist_mc_down[h].GetBinContent(i+1)-val
	
	if val!=0:
            if up==0:
                hist_mc_up[h].SetBinContent(i+1, 0)
            else:
                hist_mc_up[h].SetBinContent(i+1, float(up)/val*100)
            if down==0:
                hist_mc_down[h].SetBinContent(i+1,0)
            else:
                hist_mc_down[h].SetBinContent(i+1, float(down)/val*100)
	elif val==0:
	    hist_mc_up[h].SetBinContent(i+1, 0)
	    hist_mc_down[h].SetBinContent(i+1,0)

    max_nom = hist_mc[h].GetMaximum()
    max_up = hist_mc_up[h].GetMaximum()
    max_down = hist_mc_down[h].GetMaximum()
    min_nom = hist_mc[h].GetMinimum()
    min_up = hist_mc_up[h].GetMinimum()
    min_down = hist_mc_down[h].GetMinimum()
    newmax = 0
    newmin = 0
    if (max_up > max_nom and max_up > max_down):
        newmax = max_up
    elif (max_down > max_nom and max_down > max_up):
	newmax = max_down
    elif (max_nom > max_up and max_nom > max_down):
	newmax = max_nom
    if (min_up < min_nom and min_up < min_down):
        newmin = min_up
    elif (min_down < min_nom and min_down < min_up):
        newmin = min_down
    elif (min_nom < min_up and min_nom < min_down):
        newmin = min_nom
    if (abs(newmax) >= abs(newmin)): edge.append(abs(newmax))
    elif (abs(newmin) > abs(newmax)): edge.append(abs(newmin))
    else: print 'Did not find histo max and min'

outputdir = './results/'+year+'/systematics/'

for index in range(len(hist_mc)):
    name = obsname+stimebin+'_'+hist_mc[index].GetName()+'_'+systematic
    canvas = TCanvas(name, name)
    hist_mc[index].SetAxisRange(-edge[index]*1.7, edge[index]*1.7, "Y")
    #hist_mc[index].Draw('')
    #maxi = hist_mc[index].GetMaximum()
    #hist_mc[index].GetYaxis().SetRangeUser(-maxi-2*maxi/3.,maxi+2*maxi/3.)
    hist_mc[index].Draw('HIST')
    hist_mc_up[index].Draw('HIST SAME')
    hist_mc_down[index].Draw('HIST SAME')
    legend[index].Draw('SAME')
    #canvas.SaveAs(outputdir+stimed+'_'+name+'_'+year+'.png')
    canvas.SaveAs(outputdir+stimed+'_'+name+'_'+year+'.pdf')

#raw_input('exit')
