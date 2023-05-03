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

#triggerOption = 0 #Full trigger syst uncertainties including Nvtx partition
triggerOption = 1 #Trigger syst uncertainties without Nvtx partition
#triggerOption = 2 #Trigger syst uncertainties treated as uncorrelated in time

#puOption = "puinc" #inc pu
#puOption = "punew" #new pu
#puOption = "puold" #old pu
puOption = "putime" #pu per time bin

doScaleLumiTime = True
#doScaleLumiTime = False

#doShapeOnly = True
doShapeOnly = False

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
## Functions
################################################################################

def rename(th1, name):
    th1.SetName(name)
    th1.SetTitle(name)


def applyNominalNorm(histo, hnom):
    print 'applyNominalNorm to '+histo.GetName()
    area = histo.Integral()
    histo.Scale(hnom.Integral()/area)


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
triggersystNovtx_file = TFile('./inputs/timed/TriggerSF_'+year+'_noNvtx_new.root') #new includes PU reweighting per time bin

if triggerOption==0:
    hist_triggerSF_stat = triggersyst_file.Get('h_SF_emu_sidereel_nominal')
    hist_triggerSF_syst = triggersyst_file.Get('h_SF_emu_sidereel_Full_UncBand')
if triggerOption==1:
    if (year=='2016'):
	hist_triggerSF_stat = triggersystNovtx_file.Get('h_SF_emu_sidereel_nominal')
	hist_triggerSF_syst = triggersystNovtx_file.Get('h_SF_emu_sidereel_Full_UncBand')
    if (year=='2017'):
        hist_triggerSF_stat = triggersystNovtx_file.Get('h_SF_emu_sidereal_nominal')
	hist_triggerSF_syst = triggersystNovtx_file.Get('h_SF_emu_sidereal_Full_UncBand')

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
    for n in range(nbin):
	if triggerOption==0 or triggerOption==1:
	    mc_file = TFile('./results/'+year+'/flattree/'+observable+'_put'+str(n)+'.root')
 	    #mc_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_put'+str(n)+'.root'))
        if triggerOption==2:
	    mc_file = TFile('./results/'+year+'/flattree/'+observable+'_inclusive_put'+str(n)+'.root')
	    #mc_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_inclusive_put'+str(n)+'.root'))
	mc_file_time.append(mc_file)


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
	    #hist = TH1F(l.GetName(),l.GetName(),hist_tmp.GetNbinsX(),hist_tmp.GetXaxis().GetXmin(),hist_tmp.GetXaxis().GetXmax())
	    #for ibin in range(hist_tmp.GetNbinsX()):
	    #    hist.SetBinContent(ibin+1,hist_tmp.GetBinContent(ibin+1))
	    #if (l.GetName()=='signal'):
	    print l.GetName()+' '+str(hist.Integral())
	    hist.Scale(1./nbin)
	    histograms.append(hist)
	    #histograms[-1].Scale(1./nbin)
	    for proc in ttbar_list:
		if l.GetName()==proc:
		    print l.GetName()+' '+str(hist.Integral())
		    h_nom.append(histograms[-1])
	    del hist

        if doShapeOnly==True:
            curname = histograms[-1].GetName()
            for h in h_nom:
                proc = h.GetName()
                if TString(curname).Contains(proc):
                    applyNominalNorm(histograms[-1], h)

	if TString(l.GetName()).Contains('syst_b') or TString(l.GetName()).Contains('syst_l') or TString(l.GetName()).Contains('syst_em_trig') or TString(l.GetName()).Contains('stat'):
	#if TString(l.GetName()).Contains('syst_b_uncorrelated') or TString(l.GetName()).Contains('syst_l_uncorrelated') or TString(l.GetName()).Contains('syst_em_trig') or TString(l.GetName()).Contains('stat'):
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

	#print h_nom
	#print 'Number of nominal histos: ' +str(len(h_nom))

	if TString(l.GetName()).Contains('syst_qcdscale') or TString(l.GetName()).Contains('syst_ps_isr') or TString(l.GetName()).Contains('syst_ps_fsr') or TString(l.GetName()).Contains('syst_pdfas') or TString(l.GetName()).Contains('syst_pt_top'): #or TString(l.GetName()).Contains('mtop'):
	    #curname = histograms[-1].GetName()
	    for h in h_nom:
		#proc = h.GetName()
		if TString(histograms[-1].GetName()).Contains(h.GetName()):
		    #print curname
		    if (doShapeOnly==False):
			applyNominalNorm(histograms[-1], h)
		    #area = histograms[-1].Integral()
		    #histograms[-1].Scale(h.Integral()/area)

	#Check
	for h in h_nom:
            if TString(histograms[-1].GetName()).Contains(h.GetName()):
		print histograms[-1].GetName()+' relative_uncertainty='+str(histograms[-1].Integral()/h.Integral())
	
    h_nom_time.append(h_nom)
    histograms_time.append(histograms)
    #mc_file.Close()


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
    for n in range(nbin):
        if triggerOption==0 or triggerOption==1:
	    alt_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_color_reco_put'+str(n)+'.root'))
        if triggerOption==2:
	    alt_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_color_reco_inclusive_put'+str(n)+'.root'))


for n in range(len(alt_file_time)):

    alt_file = alt_file_time[n]
    for l in alt_file.GetListOfKeys():

	if TString(l.GetName()).Contains('responseMatrix'):
	    continue

	h = alt_file.Get(l.GetName())
	h.Scale(1./nbin)
	histograms_time[n].append(h)

	if TString(l.GetName()).Contains('mtop'):
	    curname = histograms_time[n][-1].GetName()
	    for h in h_nom_time[n]: #Was: for h in h_nom (buggy)
		proc = h.GetName()
		if TString(curname).Contains(proc):
		    applyNominalNorm(histograms_time[n][-1], h)
		    #area = histograms_time[n][-1].Integral()
		    #histograms_time[n][-1].Scale(h.Integral()/area)

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
    for n in range(nbin):
        if triggerOption==0 or triggerOption==1:
	    jec_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_jec_put'+str(n)+'.root'))
        if triggerOption==2:
	    jec_file_time.append(TFile('./results/'+year+'/flattree/'+observable+'_jec_inclusive_put'+str(n)+'.root'))

for n in range(len(jec_file_time)):

    jec_file = jec_file_time[n]
    for l in jec_file.GetListOfKeys():

	if TString(l.GetName()).Contains('responseMatrix'):
	    continue

	hj = jec_file.Get(l.GetName())
	hname = l.GetName()
	#if TString(hname).Contains('Total') or TString(hname).Contains('Absolute') or TString(hname).Contains('FlavorQCD') or TString(hname).Contains('BBEC1') or TString(hname).Contains('RelativeBal') or TString(hname).Contains('RelativeSample'):
	if(hname.find('_up')!= -1):
	    newname = hname[:-3]+'_jecUp'
	elif(hname.find('_down')!= -1):
	    newname = hname[:-5]+'_jecDown' 
	rename(hj,newname)
	hj.Scale(1./nbin)
	histograms_time[n].append(hj)

        if TString(hname).Contains("AbsoluteStat") or TString(hname).Contains("RelativeJEREC1") or TString(hname).Contains("RelativeJEREC2") or TString(hname).Contains("RelativePtEC1") or TString(hname).Contains("RelativePtEC2") or TString(hname).Contains("RelativeStatEC") or TString(hname).Contains("RelativeStatFSR") or TString(hname).Contains("RelativeStatHF") or TString(hname).Contains("TimePtEta"): #Careful for RelativeSample
            curname = histograms_time[n][-1].GetName()
            found = curname.find('_jecUp')
            if (found==-1):
                found = curname.find('_jecDown')
            newname = curname[:found] + '_' + year + curname[found:]
            histograms_time[n][-1].SetName(newname)
            histograms_time[n][-1].SetTitle(newname)

        if doShapeOnly==True:
            curname = histograms_time[n][-1].GetName()
            for h in h_nom_time[n]:
                proc = h.GetName()
                if TString(curname).Contains(proc):
                    applyNominalNorm(histograms_time[n][-1], h)

################################################################################
## MC stat
################################################################################

for n in range(nbin):

    if puOption=="puold" or puOption=="punew" or puOption=="puinc":
        timebin = 0
    else:
        timebin = n

    for h_nom in h_nom_time[timebin]:
        for iobs in range(h_nom.GetNbinsX()):
            hist_up = h_nom.Clone()
            hist_down = h_nom.Clone()
            newname = h_nom.GetName()+'_MCstat_binobs'+str(iobs)+'_'+h_nom.GetName()+'_'+year
            hist_up.SetName(newname+'Up')
            hist_down.SetName(newname+'Down')
            binval = h_nom.GetBinContent(iobs + 1)
            binerr = h_nom.GetBinError(iobs + 1)
	    if binval!=0:
                #print newname+' itime='+str(n)+' rel_err='+str(binerr/binval)
                hist_up.SetBinContent(iobs + 1, binval*(1+binerr/binval))
                hist_down.SetBinContent(iobs + 1, binval/(1+binerr/binval))
            else:
                #print newname+' itime='+str(n)+' val=0, abs_err='+str(binerr)
                hist_up.SetBinContent(iobs + 1, binerr)
                hist_down.SetBinContent(iobs + 1, 0)
	    histograms_time[timebin].append(hist_up)
            histograms_time[timebin].append(hist_down)


################################################################################
## Produce histo with time uncertainties 
################################################################################

print systematic_time_list

integral_proc_central = []
integral_proc_syst_up = []
integral_proc_syst_down = []

systlist = []
for s in systematic_time_list:
    systlist.append(0)

for proc in ttbar_list:
    integral_proc_central.append(0)

integral_proc_syst_up = [[0] * len(systematic_time_list) for _ in range(len(ttbar_list))]
integral_proc_syst_down = [[0] * len(systematic_time_list) for _ in range(len(ttbar_list))]

#    integral_proc_syst_up.append(systlist)
#    integral_proc_syst_down.append(systlist)

hist_tim_time = []

for n in range(nbin):
    if puOption=="puold" or puOption=="punew" or puOption=="puinc":
	timebin = 0
    else:
	timebin = n  

## Produce histo with time uncertainties    

    if doShapeOnly:
	for h_nom in h_nom_time[timebin]:
	    for num in range(len(ttbar_list)):
		if h_nom.GetName()==ttbar_list[num]:
		    trg_factor = 1
		    if triggerOption==0 or triggerOption==1:
			trg_factor = hist_triggerSF_syst.GetBinContent(timebin+1)
		    integral_proc_central[num] = integral_proc_central[num] + h_nom.Integral() * trg_factor

    hist_tim = []
    for numproc in range(len(ttbar_list)):
        for numsyst in range(len(systematic_time_list)):
	    for h_nom in h_nom_time[timebin]:
		if h_nom.GetName()==ttbar_list[numproc]:
		    hist_up = h_nom.Clone()
		    hist_down = h_nom.Clone()
            #hist_up = mc_file_time[timebin].Get(ttbar_list[numproc]).Clone()
            #hist_down = mc_file_time[timebin].Get(ttbar_list[numproc]).Clone()
            if (systematic_time_list[numsyst].find('emu_trig')!=-1 or systematic_time_list[numsyst]=='lumi_stability' or systematic_time_list[numsyst]=='lumi_linearity'):
                newprefix = ttbar_list[numproc]+'_'+systematic_time_list[numsyst]+'_'+year
            else:
                newprefix = ttbar_list[numproc]+'_'+systematic_time_list[numsyst]
            if systematic_time_list[numsyst].find('emu_trig')!=-1:
                hist_up.SetName(newprefix+'Up')
                hist_up.SetTitle(newprefix+'Up')
                hist_down.SetName(newprefix+'Down')
                hist_down.SetTitle(newprefix+'Down')
		if triggerOption==0 or triggerOption==1:
		    if systematic_time_list[numsyst]=='emu_trig_syst':
                        hist_up.Scale(hist_triggerSF_syst.GetBinContent(n+1)+hist_triggerSF_syst.GetBinError(n+1))
                        hist_down.Scale(hist_triggerSF_syst.GetBinContent(n+1)-hist_triggerSF_syst.GetBinError(n+1))
		    if systematic_time_list[numsyst]=='emu_trig_stat':
                        hist_up.Scale(hist_triggerSF_syst.GetBinContent(n+1)+hist_triggerSF_stat.GetBinError(n+1))
                        hist_down.Scale(hist_triggerSF_syst.GetBinContent(n+1)-hist_triggerSF_stat.GetBinError(n+1))
            else:
                hist_up.Scale(time_syst_up[systematic_time_list[numsyst]].GetBinContent(n+1))
                hist_up.SetName(newprefix+'Up')
                hist_up.SetTitle(newprefix+'Up')
                hist_down.Scale(time_syst_down[systematic_time_list[numsyst]].GetBinContent(n+1))
                hist_down.SetName(newprefix+'Down')
                hist_down.SetTitle(newprefix+'Down')
		if triggerOption==0 or triggerOption==1:
                    hist_up.Scale(hist_triggerSF_syst.GetBinContent(n+1))
                    hist_down.Scale(hist_triggerSF_syst.GetBinContent(n+1))

    	    if doShapeOnly:
		integral_proc_syst_up[numproc][numsyst] = integral_proc_syst_up[numproc][numsyst] + hist_up.Integral()
		integral_proc_syst_down[numproc][numsyst] = integral_proc_syst_down[numproc][numsyst] + hist_down.Integral()
	        #print 'A UncUp integral over time bins '+ttbar_list[numproc]+' '+systematic_time_list[numsyst]+' '+str(integral_proc_syst_up[numproc][numsyst])
        	#print 'A UncDown integral over time bins '+ttbar_list[numproc]+' '+systematic_time_list[numsyst]+' '+str(integral_proc_syst_down[numproc][numsyst])

            hist_tim.append(hist_up)
            hist_tim.append(hist_down)

    hist_tim_time.append(hist_tim)

#if doShapeOnly:
    #for procnum in range(len(ttbar_list)):
	#print 'B Nominal integral over time bins, '+ttbar_list[procnum]+' '+ str(integral_proc_central[procnum])
	#for systnum in range(len(systematic_time_list)):
	    #print 'B UncUp integral over time bins '+ttbar_list[procnum]+' '+systematic_time_list[systnum]+' '+str(integral_proc_syst_up[procnum][systnum])
	    #print 'B UncDown integral over time bins '+ttbar_list[procnum]+' '+systematic_time_list[systnum]+' '+str(integral_proc_syst_down[procnum][systnum])


if doShapeOnly: #Normalize time systematics in time
    for n in range(nbin):
	if puOption=="puold" or puOption=="punew" or puOption=="puinc":
	    timebin = 0
	else:
	    timebin = n

	for h_tim in hist_tim_time[timebin]:
	    for procnum in range(len(ttbar_list)):
                for systnum in range(len(systematic_time_list)):
		    if TString(h_tim.GetName()).Contains(systematic_time_list[systnum]) and TString(h_tim.GetName()).Contains(ttbar_list[procnum]):
			if TString(h_tim.GetName()).Contains('Up'): 
			    #print 'Scaling '+h_tim.GetName()+' by '+str(integral_proc_central[procnum]/integral_proc_syst_up[procnum][systnum])
			    h_tim.Scale(integral_proc_central[procnum]/integral_proc_syst_up[procnum][systnum])
			if TString(h_tim.GetName()).Contains('Down'):
                            h_tim.Scale(integral_proc_central[procnum]/integral_proc_syst_down[procnum][systnum])
		

################################################################################
## Storing output
################################################################################


out = './combine/'+year+'/one_bin/inputs/'

for n in range(nbin):
    if puOption=="puold" or puOption=="punew" or puOption=="puinc":
        timebin = 0
    else:
        timebin = n

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
	    #or TString(h.GetName()).Contains('syst_muon_iso')
	    if (TString(h.GetName()).Contains('syst_elec_reco') or TString(h.GetName()).Contains('syst_elec_id') or TString(h.GetName()).Contains('syst_muon_id') or TString(h.GetName()).Contains('stat_muon_id') or TString(h.GetName()).Contains('stat_muon_iso') or TString(h.GetName()).Contains('syst_b_correlated') or TString(h.GetName()).Contains('syst_b_uncorrelated') or TString(h.GetName()).Contains('syst_l_correlated') or TString(h.GetName()).Contains('syst_l_uncorrelated') or TString(h.GetName()).Contains('syst_prefiring') or (TString(h.GetName()).Contains('jec') and TString(h.GetName()).Contains('FlavorPureBottom')==False and TString(h.GetName()).Contains('FlavorPureGluon')==False) or TString(h.GetName()).Contains('syst_em_trig')) or (puOption!="putime" and TString(h.GetName()).Contains('syst_pu')):
                curname = h.GetName()
                found = curname.find('Up')
                if (found==-1):
                    found = curname.find('Down')
                newname = curname[:found] + '_t' + str(n) + curname[found:]
	        h_new.SetName(newname)
	if triggerOption==0 or triggerOption==1:
	    h_new.Scale(hist_triggerSF_syst.GetBinContent(n+1)) #Is it correct?
	if doScaleLumiTime:
	    h_new.Scale(hist_lumi_corr.GetBinContent(n+1))
	#h_new.Scale(1./nbin)
        h_new.Write()

	#Check
        for hproc in h_nom_time[n]:
            if TString(h_new.GetName()).Contains(hproc.GetName()):
		if (abs(h_new.Integral()/hproc.Integral()-1)>1):
                    print 'Check: '+h_new.GetName()+' integral='+str(h_new.Integral())+' relative_uncertainty='+str(h_new.Integral()/hproc.Integral())

    for h in hist_tim_time[timebin]:
        #h.Scale(hist_triggerSF_syst.GetBinContent(n+1)) #Done before
	if doExpTimeNuisance and TString(h.GetName()).Contains('emu_trig_stat'):
	    curname = h.GetName()
	    found = curname.find('Up')
            if (found==-1):
                found = curname.find('Down')
            newname = curname[:found] + '_t' + str(n) + curname[found:]
            h.SetName(newname)
        if doScaleLumiTime:
            h.Scale(hist_lumi_corr.GetBinContent(n+1))
	#h.Scale(1./nbin)
        h.Write()

        #Check
        for hproc in h_nom_time[n]:
            if TString(h.GetName()).Contains(hproc.GetName()):
		if (abs(h.Integral()/hproc.Integral()-1)>1):
                    print 'Check: '+h.GetName()+' integral='+str(h.Integral())+' relative_uncertainty='+str(h.Integral()/hproc.Integral())

    output.Close()
    #del hist_tim


file_txt = ''
for i in data_number_of_event:
    file_txt += str(i)+'\n'

file = open('./combine/'+year+'/'+observable+'_noe_data_timed.txt','w') 
file.write(file_txt) 
file.close() 
