#!/usr/bin/env python

import os
import sys
sys.path.append('./')


import argparse

from tools.sample_manager import *
from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

nbin = 24

doExpTimeNuisance = True

triggerOption = 0 #Full trigger syst uncertainties including Nvtx partition
#triggerOption = 1 #Trigger syst uncertainties without Nvtx partition
#triggerOption = 2 #Trigger syst uncertainties treated as uncorrelated in time

#puOption = "puinc" #inc pu
#puOption = "punew" #new pu
puOption = "puold" #old pu
#puOption = "putime" #pu per time bin

doScaleLumiTime = True

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

#histograms = []
data_number_of_event = []


################################################################################
## Timed histo
################################################################################

time_syst_up = {}
time_syst_down = {}

#lumisyst_file = TFile('./inputs/timed/AllTimedSyst_'+year+'.root')
lumisyst_file = TFile('./inputs/timed/LumiUncertainties_'+year+'.root')
#hist_lumi_corr = lumisyst_file.Get('hInstLumi_DataScaleFactor')
hist_lumi_corr = lumisyst_file.Get('hIntegratedLumi_sidereal_mcweight')

#triggersyst_file = TFile('./inputs/timed/AllTimedSyst_'+year+'.root')
triggersyst_file = TFile('./inputs/timed/TriggerSF_'+year+'.root')
triggersystNovtx_file = TFile('./inputs/timed/TriggerSF_'+year+'_noNvtx.root')

if triggerOption==0:
    hist_triggerSF = triggersyst_file.Get('h_SF_emu_sidereel_Full_UncBand')
if triggerOption==1:
    hist_triggerSF = triggersystNovtx_file.Get('h_SF_emu_sidereel_Full_UncBand')

for l in lumisyst_file.GetListOfKeys():
    if not TString(l.GetName()).Contains('DataScaleFactor'):
        if TString(l.GetName()).Contains('_Up'):
            if TString(l.GetName()).Contains('Inclusive'):
                time_syst_up.update({systematic_time_list[0]: lumisyst_file.Get(l.GetName())})
            elif TString(l.GetName()).Contains('Stability'):
                time_syst_up.update({systematic_time_list[1]: lumisyst_file.Get(l.GetName())})
            elif TString(l.GetName()).Contains('Linearity'):
                time_syst_up.update({systematic_time_list[2]: lumisyst_file.Get(l.GetName())})
        elif TString(l.GetName()).Contains('_Down'):
            if TString(l.GetName()).Contains('Inclusive'):
                time_syst_down.update({systematic_time_list[0]: lumisyst_file.Get(l.GetName())})
            elif TString(l.GetName()).Contains('Stability'):
                time_syst_down.update({systematic_time_list[1]: lumisyst_file.Get(l.GetName())})
            elif TString(l.GetName()).Contains('Linearity'):
                time_syst_down.update({systematic_time_list[2]: lumisyst_file.Get(l.GetName())})
    #if TString(l.GetName()).Contains('SF_emu_sidereel_Full_UncBand'):
    #    time_syst_up.update({systematic_time_list[3]: lumisyst_file.Get(l.GetName())})


################################################################################
## MC histo
################################################################################

mc_file_time = []

if puOption=="puold" or puOption=="punew" or puOption=="puinc":
    if triggerOption==0 or triggerOption==1: 
	mc_file = TFile('./results/'+year+'/flattree/'+observable+'_'+puOption+'.root')
    if triggerOption==2:
	mc_file = TFile('./results/'+year+'/flattree/'+observable+'_inclusive_'+puOption+'.root')
    mc_file_time.append(mc_file)
else:
    for n in range(24):
	if triggerOption==0 or triggerOption==1:
 	    mc_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_put'+str(n)+'.root'))
        if triggerOption==2:
	    mc_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_inclusive_put'+str(n)+'.root'))


histograms_time = []
h_nom_time = []
#h_nom = []

for mc_file in mc_file_time:

    print mc_file.GetName()

    h_nom = []
    histograms = []

    for l in mc_file.GetListOfKeys():

	if TString(l.GetName()).Contains('responseMatrix'):
	    continue

	if not TString(l.GetName()).Contains('data_obs'):
	    hist = mc_file.Get(l.GetName())
	    hist.Scale(1./nbin)
	    histograms.append(hist)
	    for proc in ttbar_list:
		if l.GetName()==proc:
		    h_nom.append(mc_file.Get(l.GetName()))

	if TString(l.GetName()).Contains('syst_b_uncorrelated') or TString(l.GetName()).Contains('syst_l_uncorrelated') or TString(l.GetName()).Contains('syst_em_trig'):
	    curname = histograms[-1].GetName()
	    found = curname.find('Up')
	    if (found==-1): 
		found = curname.find('Down')
	    newname = curname[:found] + '_' + year + curname[found:]
	    histograms[-1].SetName(newname)

	if TString(l.GetName()).Contains('syst_qcdscale') or TString(l.GetName()).Contains('syst_ps_isr'):
	    curname = histograms[-1].GetName()
	    found = curname.find('Up')
	    if (found==-1):
		found = curname.find('Down')
	    for h in h_nom:
		proc = h.GetName()
		if TString(curname).Contains(proc):
		    newname = curname[:found] + '_' + proc + curname[found:]
		    histograms[-1].SetName(newname)
		    #area = histograms[-1].Integral()
		    #print(curname+"Before nom="+str(h.Integral())+" qcdscale="+str(area))
		    #if TString(l.GetName()).Contains('syst_qcdscale'):
		    #    histograms[-1].Scale(h.Integral()/area)
		    #print(curname+"After nom="+str(h.Integral())+" qcdscale="+str(histograms[-1].Integral()))

	if TString(l.GetName()).Contains('syst_qcdscale') or TString(l.GetName()).Contains('syst_ps_isr') or TString(l.GetName()).Contains('syst_ps_fsr') or TString(l.GetName()).Contains('syst_pdfas') or TString(l.GetName()).Contains('syst_pt_top') or TString(l.GetName()).Contains('mtop'):
	    curname = histograms[-1].GetName()
	    for h in h_nom:
		proc = h.GetName()
		if TString(curname).Contains(proc):
		    print curname
		    area = histograms[-1].Integral()
		    histograms[-1].Scale(h.Integral()/area)
	
    h_nom_time.append(hnom)
    histograms_time.append(histograms)

################################################################################
## ALT histo
################################################################################

alt_file_time = []

if puOption=="puold" or puOption=="punew" or puOption=="puinc":
    if triggerOption==0 or triggerOption==1:
	alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco_'+puOption+'.root')
    if triggerOption==2:
	alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco_inclusive_'+puOption+'.root')
    alt_file_time.append(alt_file)
else:
    for n in range(24):
        if triggerOption==0 or triggerOption==1:
	    alt_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_color_reco_put'+str(n)+'.root'))
        if triggerOption==2:
	    alt_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_color_reco_inclusive_put'+str(n)+'.root'))


for n in len(range(alt_file_time)):

    alt_file = alt_file_time[n]
    for l in alt_file.GetListOfKeys():

	if TString(l.GetName()).Contains('responseMatrix'):
	    continue

	h = alt_file.Get(l.GetName())
	h.Scale(1./nbin)
	histograms_time[n].append(h)

	if TString(l.GetName()).Contains('mtop'):
	    curname = histograms_time[n][-1].GetName()
	    for h in h_nom:
		proc = h.GetName()
		if TString(curname).Contains(proc):
		    area = histograms_time[n][-1].Integral()
		    histograms_time[n][-1].Scale(h.Integral()/area)

################################################################################
## JEC histo
################################################################################

jec_file_time = []

if puOption=="puold" or puOption=="punew" or puOption=="puinc":
    if triggerOption==0 or triggerOption==1:
	jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec_'+puOption+'.root')
    if triggerOption==2:
	jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec_inclusive_'+puOption+'.root')
    jec_file_time.append(jec_file)
else:
    for n in range(24):
        if triggerOption==0 or triggerOption==1:
	    jec_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_jec_put'+str(n)+'.root'))
        if triggerOption==2:
	    jec_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_jec_inclusive_put'+str(n)+'.root'))

for n in len(range(jec_file_time)):

    jec_file = jec_file_time[n]
    for l in jec_file.GetListOfKeys():

	if TString(l.GetName()).Contains('responseMatrix'):
	    continue

	hj = jec_file.Get(l.GetName())
	hname = l.GetName()
	if TString(hname).Contains('Total') or TString(hname).Contains('Absolute') or TString(hname).Contains('FlavorQCD') or TString(hname).Contains('BBEC1') or TString(hname).Contains('RelativeBal') or TString(hname).Contains('RelativeSample'):
	    if(hname.find('_up')!= -1):
		newname = hname[:-3]+'_jecUp'
	    elif(hname.find('_down')!= -1):
		newname = hname[:-5]+'_jecDown' 
	#if(hname.find('TotalUp')!= -1):
	#    name = hname[:-7]+'jecUp'
	#elif(hname.find('TotalDown')!= -1):
	#    name = hname[:-9]+'jecDown'
	rename(hj,newname)
	hj.Scale(1./nbin)
	histograms_time[n].append(hj)

#print histograms[27].GetName(), histograms[27].GetTitle(), len(histograms)


################################################################################
## Produce histo with time uncertainties + store output
################################################################################


out = './combine/'+year+'/one_bin/inputs/'

print systematic_time_list

for n in range(nbin):

    if puOption=="puold" or puOption=="punew" or puOption=="puinc":
	timebin = 0
    else:
	timebin = n  

## Produce histo with time uncertainties    

    hist_tim = []

    for g in ttbar_list:
        for s in systematic_time_list:
            hist_up = mc_file_time[timebin].Get(g).Clone()
            hist_down = mc_file_time[timebin].Get(g).Clone()
	    hist_up.Scale(1./nbin) #Was forgotten. Does it matter in combine?
            hist_down.Scale(1./nbin) #Was forgotten. Does it matter in combine?
            if (s == 'emu_trig' or s=='lumi_stability' or s=='lumi_linearity'):
                newprefix = g+'_'+s+'_'+year
            else:
                newprefix = g+'_'+s
            if s == 'emu_trig':
                hist_up.SetName(newprefix+'Up')
                hist_up.SetTitle(newprefix+'Up')
                hist_down.SetName(newprefix+'Down')
                hist_down.SetTitle(newprefix+'Down')
		if triggerOption==0 or triggerOption==1:
                    hist_up.Scale(hist_triggerSF.GetBinContent(n+1)+hist_triggerSF.GetBinError(n+1))
                    hist_down.Scale(hist_triggerSF.GetBinContent(n+1)-hist_triggerSF.GetBinError(n+1))
            else:
                hist_up.Scale(time_syst_up[s].GetBinContent(n+1))
                hist_up.SetName(newprefix+'Up')
                hist_up.SetTitle(newprefix+'Up')
                hist_down.Scale(time_syst_down[s].GetBinContent(n+1))
                hist_down.SetName(newprefix+'Down')
                hist_down.SetTitle(newprefix+'Down')
		if triggerOption==0 or triggerOption==1:
                    hist_up.Scale(hist_triggerSF.GetBinContent(n+1))
                    hist_down.Scale(hist_triggerSF.GetBinContent(n+1))
            hist_tim.append(hist_up)
            hist_tim.append(hist_down)
    '''
        if n == 0:
            print g, s
            print histograms[28].GetName(), histograms[28].GetTitle(), len(histograms)
    '''

    ## Store data output
    data_file = TFile('./results/'+year+'/flattree/'+observable+'_data_timed'+str(nbin)+'.root')
    hist_data_time = data_file.Get('data_obs_bin'+str(n))
    hist_data_time.SetName('data_obs')
    hist_data_time.SetTitle('data_obs')
    data_number_of_event.append(hist_data_time.Integral())
    print 'Integral bin '+str(n)+' : '+str(data_number_of_event[n])

    output = TFile(out+observable+'_'+str(nbin)+'_'+str(n)+'.root', "RECREATE")
    hist_data_time.Write()

    ## Store mc output
    for h in histograms_time[timebin]:
	h_new = h.Clone()
	if doExpTimeNuisance:
	    if (TString(h.GetName()).Contains('syst_elec_reco') or TString(h.GetName()).Contains('syst_elec_id') or TString(h.GetName()).Contains('syst_muon_id') or TString(h.GetName()).Contains('syst_muon_iso') or TString(h.GetName()).Contains('syst_pu') or TString(h.GetName()).Contains('syst_b_correlated') or TString(h.GetName()).Contains('syst_b_uncorrelated') or TString(h.GetName()).Contains('syst_l_correlated') or TString(h.GetName()).Contains('syst_l_uncorrelated') or TString(h.GetName()).Contains('syst_prefiring') or TString(h.GetName()).Contains('jec') or TString(h.GetName()).Contains('syst_em_trig')):
                curname = h.GetName()
                found = curname.find('Up')
                if (found==-1):
                    found = curname.find('Down')
                newname = curname[:found] + '_t' + str(n) + curname[found:]
	        h_new.SetName(newname)
	if triggerOption==0 or triggerOption==1:
	    h_new.Scale(hist_triggerSF.GetBinContent(n+1)) #Is it correct?
	if doScaleLumiTime:
	    h_new.Scale(hist_lumi_corr.GetBinContent(n+1))
        h_new.Write()

    for h in hist_tim:
        #h.Scale(hist_triggerSF.GetBinContent(n+1)) #Is it correct?
        if doScaleLumiTime:
            h_new.Scale(hist_lumi_corr.GetBinContent(n+1))
        h.Write()

    output.Close()
    del hist_tim


file_txt = ''
for i in data_number_of_event:
    file_txt += str(i)+'\n'

file = open('./combine/'+year+'/'+observable+'_noe_data_timed.txt','w') 
file.write(file_txt) 
file.close() 
