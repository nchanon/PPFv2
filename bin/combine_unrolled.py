#!/usr/bin/env python

import os
import sys
sys.path.append('./')


import argparse

from tools.sample_manager import *
from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

nbintime = 24

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

cmunu = 0.001
doSMEkindep = True
#doSMEkindep = False

doResponseMatrix = True
#doResponseMatrix = False

doDirectRecoSignal = True
#doDirectRecoSignal = False

doShapeOnly = False
#doShapeOnly = True

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
wilson_ = args.wilson

TH1.SetDefaultSumw2(1)

################################################################################
## function
################################################################################

def rename(th1, name):
    th1.SetName(name)
    th1.SetTitle(name)

def applyNominalNorm(histo, hnom):
    area = histo.Integral()
    histo.Scale(hnom.Integral()/area)

################################################################################
## Code body
################################################################################

###########
# data part
###########

data_file = TFile('./results/'+year+'/flattree/'+observable+'_data_timed24.root')
data_integral = 0
hist_data_in = []

for i in range(nbintime): 
    hist_data_in.append(data_file.Get('data_obs_bin'+str(i)))
    data_integral += hist_data_in[i].Integral()

# convenient variables
nbinkin = hist_data_in[0].GetNbinsX()
#min_bin = hist_data_in[0].GetXaxis().GetXmin()
#max_bin = hist_data_in[0].GetXaxis().GetXmax()

canvas_data = TCanvas('data_obs', 'data_obs', 1000, 800)
hist_data = TH1F('data_obs','data_obs', nbinkin*nbintime,  0, nbintime)
for itime in range(nbintime):
    for n in range(nbinkin):
        hist_data.SetBinContent(n + itime*nbinkin + 1 , hist_data_in[itime].GetBinContent(n+1))
        hist_data.SetBinError(n + itime*nbinkin + 1 , hist_data_in[itime].GetBinError(n+1))

###########
# MC part
###########

lumisyst_file = TFile('./inputs/timed/LumiUncertainties_'+year+'.root')
hist_lumi_corr = lumisyst_file.Get('hIntegratedLumi_sidereal_mcweight')

#lumi_file = TFile('./inputs/timed/AllTimedSyst_'+year+'.root')
triggersyst_file = TFile('./inputs/timed/TriggerSF_'+year+'.root')
triggersystNovtx_file = TFile('./inputs/timed/TriggerSF_'+year+'_noNvtx_new.root')

if triggerOption==0:
    hist_triggerSF_stat = triggersyst_file.Get('h_SF_emu_sidereel_nominal')
    hist_triggerSF_syst = triggersyst_file.Get('h_SF_emu_sidereel_Full_UncBand')
if triggerOption==1:
    if year=='2016':
        hist_triggerSF_stat = triggersystNovtx_file.Get('h_SF_emu_sidereel_nominal')
        hist_triggerSF_syst = triggersystNovtx_file.Get('h_SF_emu_sidereel_Full_UncBand')
    if year=='2017':
        hist_triggerSF_stat = triggersystNovtx_file.Get('h_SF_emu_sidereal_nominal')
        hist_triggerSF_syst = triggersystNovtx_file.Get('h_SF_emu_sidereal_Full_UncBand')


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

if puOption=="puold" or puOption=="punew" or puOption=="puinc":
    nputimebin = 1
else:
    nputimebin = nbintime

#Files for response matrices in triggerOption==2: do not use inclusive, or compute it also for inclusive?

histograms_time = []
h_nom_time = []
histograms_time_responseMatrix = []

#mc_integral = 0
#hist_mc = []

for mc_file in mc_file_time:

    print mc_file.GetName()

    h_nom = [] 
    #histograms = []
    mc_integral = 0
    hist_mc = []
    hist_responseMatrix = []

    triggerSF = 1
    #index = 0
    for l in mc_file.GetListOfKeys():

        for proc in ttbar_list:
            if l.GetName()==proc:
                h_nom.append(mc_file.Get(l.GetName()))
		mc_integral += h_nom[-1].Integral()	

        if TString(l.GetName()).Contains('responseMatrix'):
	    hist_responseMatrix.append(mc_file.Get(l.GetName()))
            continue

	if not TString(l.GetName()).Contains('data_obs'):
	    hist_mc.append(mc_file.Get(l.GetName()))
	    #hist_mc.append(TH1F("","", nbinkin*nbintime,  0, nbintime))
	    #for itime in range(nbintime):
		#for j in range(nbinkin):
		    #if triggerOption==0 or triggerOption==1:
			#triggerSF = hist_triggerSF.GetBinContent(itime+1)
		    #if  triggerOption==2:
			#triggerSF = 1
		    #hist_mc[index].SetBinContent(j + itime*nbinkin + 1 , hist.GetBinContent(j+1)*triggerSF) #i => itime
		    #hist_mc[index].SetBinError(j + itime*nbinkin + 1 , hist.GetBinError(j+1)*triggerSF) #i ou itime => j ! how to treat SF uncertainties (for nominal MC?)
	    #del hist
	    hist_mc[-1].SetName(l.GetName())
	    hist_mc[-1].SetTitle(l.GetName())
	    #hist_mc[index].Scale(1./nbintime)

	    if (doShapeOnly==True or (TString(hist_mc[-1].GetName()).Contains('syst_qcdscale') or TString(hist_mc[-1].GetName()).Contains('syst_ps_isr') or TString(hist_mc[-1].GetName()).Contains('syst_ps_fsr') or TString(hist_mc[-1].GetName()).Contains('syst_pdfas') or TString(hist_mc[-1].GetName()).Contains('syst_pt_top'))):
	        for h_proc in h_nom:
		    proc = h_proc.GetName()
                    if TString(hist_mc[-1].GetName()).Contains(proc):
                        applyNominalNorm(hist_mc[-1], h_proc)
	    #index += 1

    h_nom_time.append(h_nom)
    histograms_time.append(hist_mc)
    histograms_time_responseMatrix.append(hist_responseMatrix)


##############
# JEC part 
##############

#hist_alt_jec = []

#if triggerOption==0 or triggerOption==1:
#    jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec.root')
#if triggerOption==2:
#    jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec_inclusive.root')

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

for n in range(len(jec_file_time)):

    jec_file = jec_file_time[n]
    for l in jec_file.GetListOfKeys():

        if TString(l.GetName()).Contains('responseMatrix'):
            histograms_time_responseMatrix[n].append(jec_file.Get(l.GetName()))
            continue

	hh = jec_file.Get(l.GetName())
	#hh.Scale(1./nbintime)
	hname = l.GetName()
	#if TString(hname).Contains('Total') or TString(hname).Contains('Absolute') or TString(hname).Contains('FlavorQCD') or TString(hname).Contains('BBEC1') or TString(hname).Contains('RelativeBal') or TString(hname).Contains('RelativeSample'):
	if(hname.find('_up')!= -1):
	    newname = hname[:-3]+'_jecUp'
	elif(hname.find('_down')!= -1):
	    newname = hname[:-5]+'_jecDown'
	h_jec = hh
	h_jec.SetTitle(newname)
	h_jec.SetName(newname)
	#for i in range(nbintime):
	#    for j in range(nbinkin):
	#	if triggerOption==0 or triggerOption==1:
	#	    triggerSF = hist_triggerSF.GetBinContent(i+1)
	#	if  triggerOption==2:
	#	    triggerSF = 1
	#	h_jec.SetBinContent(j + i*nbinkin + 1, hh.GetBinContent(j+1)*triggerSF)
	#	h_jec.SetBinError(j + i*nbinkin + 1, hh.GetBinError(j + 1)*triggerSF)

	histograms_time[n].append(h_jec)

        if TString(hname).Contains("AbsoluteStat") or TString(hname).Contains("RelativeJEREC1") or TString(hname).Contains("RelativeJEREC2") or TString(hname).Contains("RelativePtEC1") or TString(hname).Contains("RelativePtEC2") or TString(hname).Contains("RelativeStatEC") or TString(hname).Contains("RelativeStatFSR") or TString(hname).Contains("RelativeStatHF") or TString(hname).Contains("TimePtEta"): #Careful for RelativeSample
            curname = histograms_time[n][-1].GetName()
            found = curname.find('_jecUp')
            if (found==-1):
                found = curname.find('_jecDown')
            newname = curname[:found] + '_' + year + curname[found:]
            histograms_time[n][-1].SetName(newname)
            histograms_time[n][-1].SetTitle(newname)

        if doShapeOnly==True:
            for h_proc in h_nom_time[n]:
                proc = h_proc.GetName()
                if TString(histograms_time[n][-1].GetName()).Contains(proc):
                    applyNominalNorm(histograms_time[n][-1], h_proc)


##############
# ALT part 
##############

#if triggerOption==0 or triggerOption==1:
#    alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco.root')
#if triggerOption==2:
#    alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco_inclusive.root')

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

for n in range(len(alt_file_time)):

    alt_file = alt_file_time[n]
    for l in alt_file.GetListOfKeys():
	#h_alt = alt_file.Get(l.GetName())
	#hh.Scale(1./nbintime)
	#hname = l.GetName()
	#h_alt = hh
	#h_alt = TH1F("", "", nbinkin*nbintime,  0, nbintime)
	#for i in range(nbintime):
	#    for j in range(nbinkin):
	#	if triggerOption==0 or triggerOption==1:
	#	    triggerSF = hist_triggerSF.GetBinContent(i+1)
	#	if  triggerOption==2:
	#	    triggerSF = 1
	#	h_alt.SetBinContent(j + i*nbinkin + 1, hh.GetBinContent(j+1)*triggerSF)
	#	h_alt.SetBinError(j + i*nbinkin + 1, hh.GetBinError(j + 1)*triggerSF)
	#h_alt.SetTitle(hname)
	#h_alt.SetName(hname)
	#hist_alt_jec.append(h_alt)
        histograms_time[n].append(alt_file.Get(l.GetName()))

	if TString(l.GetName()).Contains('mtop'):
            for h_proc in h_nom_time[n]:
                proc = h_proc.GetName()
                if TString(histograms_time[n][-1].GetName()).Contains(proc):
                    applyNominalNorm(histograms_time[n][-1], h_proc)

#hist_mc_all = hist_mc
#for i in range(len(hist_alt_jec)): hist_mc_all.append(hist_alt_jec[i])
#print hist_mc
#print hist_mc_all

#####################
# Make Unrolled hists
#####################

print 'Make Unrolled hists'

hist_unrolled = []
for h in histograms_time[0]:
    hist_unrolled.append(TH1F(h.GetName(),h.GetName(), nbinkin*nbintime,  0, nbintime))

for h_unrolled in hist_unrolled:
    for itime in range(nbintime):

        if triggerOption==0 or triggerOption==1:
            scale = hist_triggerSF_syst.GetBinContent(itime+1)
        elif triggerOption==2:
            scale = 1

        if doScaleLumiTime:
            scale = scale*hist_lumi_corr.GetBinContent(itime+1)

	#scale = scale / nbintime

        if nputimebin==1:
            iputime = 0
        else:
            iputime = itime

        for h in histograms_time[iputime]:
            if h_unrolled.GetName()==h.GetName():
                for j in range(nbinkin):
                    h_unrolled.SetBinContent(j + itime*nbinkin + 1 , h.GetBinContent(j+1)*scale) #i => itime
                    h_unrolled.SetBinError(j + itime*nbinkin + 1 , h.GetBinError(j+1)*scale) #i ou itime => j ! how to treat SF uncertainties (for nominal MC?)

    h_unrolled.Scale(1./nbintime)

#################
# Timed syst part
#################

lumi_syst_up = {}
lumi_syst_down = {}

for l in lumisyst_file.GetListOfKeys():
    if TString(l.GetName()).Contains('_Up'):
        if TString(l.GetName()).Contains('Inclusive'):
            lumi_syst_up.update({systematic_time_list[0]: lumisyst_file.Get(l.GetName())})
        elif TString(l.GetName()).Contains('Stability'):
            lumi_syst_up.update({systematic_time_list[1]: lumisyst_file.Get(l.GetName())})
        elif TString(l.GetName()).Contains('Linearity'):
            lumi_syst_up.update({systematic_time_list[2]: lumisyst_file.Get(l.GetName())})
    elif TString(l.GetName()).Contains('_Down'):
        if TString(l.GetName()).Contains('Inclusive'):
            lumi_syst_down.update({systematic_time_list[0]: lumisyst_file.Get(l.GetName())})
        elif TString(l.GetName()).Contains('Stability'):
            lumi_syst_down.update({systematic_time_list[1]: lumisyst_file.Get(l.GetName())})
        elif TString(l.GetName()).Contains('Linearity'):
            lumi_syst_down.update({systematic_time_list[2]: lumisyst_file.Get(l.GetName())})
    #if TString(l.GetName()).Contains('SF_emu_sidereel_Full_UncBand'):
    #    lumi_syst_up.update({systematic_time_list[3]: lumisyst_file.Get(l.GetName())})

print systematic_time_list
print lumi_syst_up
print lumi_syst_down

for l in systematic_time_list:
    for g in hist_unrolled:
	for b in ttbar_list:
	    if g.GetName() == b:
		print(b+' nom='+str(g.Integral()))
		if (l.find('emu_trig')!=-1 or l=='lumi_stability' or l=='lumi_linearity'):
		    newprefix = b+'_'+l+'_'+year
		else:
		    newprefix = b+'_'+l
		hist_up = g.Clone()
		hist_down = g.Clone()
		hist_up.SetName(newprefix+'Up')
		hist_up.SetTitle(newprefix+'Up')
		hist_down.SetName(newprefix+'Down')
		hist_down.SetTitle(newprefix+'Down')
		for itime in range(nbintime):
		    for j in range(nbinkin):
			if l.find('emu_trig')!=-1:
			    if triggerOption==0 or triggerOption==1:
				if l=='emu_trig_syst':
				    triggerSF_unc = hist_triggerSF_syst.GetBinError(itime+1)
                                if l=='emu_trig_stat':
                                    triggerSF_unc = hist_triggerSF_stat.GetBinError(itime+1)
			    if  triggerOption==2:
				triggerSF_unc = 0.
			    hist_up.SetBinContent(j + itime*nbinkin + 1 ,
						g.GetBinContent(j + itime*nbinkin + 1)*(1.+triggerSF_unc)) #hist_up => g
			    hist_up.SetBinError(j + itime*nbinkin + 1 ,
						g.GetBinError(j + itime*nbinkin + 1)*(1.+triggerSF_unc))
			    hist_down.SetBinContent(j + itime*nbinkin + 1 ,
						g.GetBinContent(j + itime*nbinkin + 1)*(1.-triggerSF_unc))
			    hist_down.SetBinError(j + itime*nbinkin + 1 ,
						g.GetBinError(j + itime*nbinkin + 1)*(1.-triggerSF_unc))
			else:
			    hist_up.SetBinContent(j + itime*nbinkin + 1 ,
						g.GetBinContent(j + itime*nbinkin + 1)*lumi_syst_up[l].GetBinContent(itime+1))
			    hist_up.SetBinError(j + itime*nbinkin + 1 ,
						g.GetBinError(j + itime*nbinkin + 1)*lumi_syst_up[l].GetBinError(itime+1))
			    hist_down.SetBinContent(j + itime*nbinkin + 1 ,
						g.GetBinContent(j + itime*nbinkin + 1)*lumi_syst_down[l].GetBinContent(itime+1))
			    hist_down.SetBinError(j + itime*nbinkin + 1 ,
						g.GetBinError(j + itime*nbinkin + 1)*lumi_syst_down[l].GetBinError(itime+1))
		print(l+' '+b+' up='+str(hist_up.Integral())+' down='+str(hist_down.Integral()))

		if doShapeOnly==True:
		    applyNominalNorm(hist_up, g)
                    applyNominalNorm(hist_down, g)

		hist_unrolled.append(hist_up)
		hist_unrolled.append(hist_down)

#################
# MC Stat uncert
#################


for g in hist_unrolled:
    for b in ttbar_list:
        if g.GetName() == b:
	    for iobs in range(nbinkin):
                hist_up = g.Clone()
                hist_down = g.Clone()
	        newname = b+'_MCstat_binobs'+str(iobs)+'_'+b+'_'+year
		hist_up.SetName(newname+'Up')
		hist_down.SetName(newname+'Down')
		for itime in range(nbintime):
		    binval = g.GetBinContent(iobs + itime*nbinkin + 1)
		    binerr = g.GetBinError(iobs + itime*nbinkin + 1)
		    if binval!=0:
		        print newname+' itime='+str(itime)+' rel_err='+str(binerr/binval)
		        hist_up.SetBinContent(iobs + itime*nbinkin + 1, binval*(1+binerr/binval))
                        hist_down.SetBinContent(iobs + itime*nbinkin + 1, binval/(1+binerr/binval))
		    else:
			hist_up.SetBinContent(iobs + itime*nbinkin + 1, binerr)
                        hist_down.SetBinContent(iobs + itime*nbinkin + 1, 0)
		hist_unrolled.append(hist_up)
		hist_unrolled.append(hist_down)

#exit()

###########
# SME part
###########

observable_forSME = observable

#if observable!="m_dilep" and observable!="pt_emu":
#    observable_forSME = "pt_emu" #Will use inclusive SME histos anyways

#if observable!="m_dilep" and observable!="pt_emu":
#    doSMEkindep = False

hist_sme = []
if doDirectRecoSignal==True:
    #sme_file_signal = TFile('./results/'+year+'/flattree/'+observable_forSME+'_sme.root')
    sme_file_signal = TFile('./results/'+year+'/flattree/sme_matrices_alt_puinc.root')
    sme_file_singletop = TFile('./results/'+year+'/flattree/'+observable_forSME+'_sme.root')
if doResponseMatrix==True and doDirectRecoSignal==False:
    #sme_file_signal = TFile('./inputs/pheno/signal_miniAOD_LO_Comb_n_bjets_13TeVCMSnanoGEN_particle.root')
    sme_file_signal = TFile('./results/'+year+'/flattree/'+observable_forSME+'_sme.root')
    sme_file_singletop = TFile('./results/'+year+'/flattree/'+observable_forSME+'_sme.root')

wilsonList = ["cLXX","cLXY","cLXZ","cLYZ","cRXX","cRXY","cRXZ","cRYZ","cXX","cXY","cXZ","cYZ","dXX","dXY","dXZ","dYZ"]
wilsonListSingleTop = ["cLXX","cLXY","cLXZ","cLYZ"]

if wilson_!="sme_all":
    wilsonList = [wilson_]

#Normalize all response matrices
nbingenkin = nbinkin+2

histograms_time_responseMatrix_norm = []
for putimebin in range(nputimebin):
    h_responseMatrix_norm = []
    for h2D in histograms_time_responseMatrix[putimebin]:
	if TString(h2D.GetName()).Contains('responseMatrix'):
	    #hist_responseMatrix.append(responseMatrix_file.Get(h.GetName()))
	    for ix in range(nbinkin):
		area = h2D.Integral(1+ix, 1+ix, 1, nbingenkin)
		for iy in range(nbingenkin):
		    bincontent = h2D.GetBinContent(ix+1,iy+1)
		    h2D.SetBinContent(ix+1, iy+1, bincontent/area)
		    #print h2D.GetName()+' ix='+str(ix)+' iy='+str(iy)+' content='+str(h2D.GetBinContent(ix+1,iy+1))
	h_responseMatrix_norm.append(h2D)
    histograms_time_responseMatrix_norm.append(h_responseMatrix_norm)

        #nbingenkin = nbinkin+2
        #Response matrix: nbingenkin=0: OOA. nbingenkin=1: 0 b-jets.
        #Modulations SME: nbingenkin=0: 0 b-jets. ...nbingenkin=5: OOA 
        #Reco bins: nbinkin=0: 1 b-jet...

for wilson in wilsonList:

	sme_sig = sme_file_signal.Get(wilson)
	sme_sig_kinbin = []
	sme_sig_kinbin_directreco = []
	sme_sig_kinbin_directreco_MCstatUp = []
	sme_sig_kinbin_directreco_MCstatDown = []
        hist_sme_MCstatUp = []
        hist_sme_MCstatDown = []

	for k in range(nbingenkin):
	    print(wilson+"_"+str(k))
	    #sme_sig_kinbin.append(sme_file_signal.Get(wilson+"_"+str(k)))
	    sme_sig_kinbin.append(sme_file_signal.Get("signal_miniAOD_LO_Comb_"+wilson+"_"+str(k)))
	for k in range(nbinkin):
	    if doDirectRecoSignal:
	        #sme_sig_kinbin_directreco.append(sme_file_signal.Get("reco_"+year+"_"+wilson+"_"+str(k+1)))
		sme_sig_kinbin_directreco.append(sme_file_signal.Get("reco_"+year+"_central_"+wilson+"_"+str(k+1)))
		sme_sig_kinbin_directreco_MCstatUp.append(sme_file_signal.Get("reco_"+year+"_central_"+wilson+"_MCstatUp_"+str(k+1)))
                sme_sig_kinbin_directreco_MCstatDown.append(sme_file_signal.Get("reco_"+year+"_central_"+wilson+"_MCstatDown_"+str(k+1)))

	for g in hist_unrolled:
	    if g.GetName().find('signal') != -1:
		hist_sme.append(g.Clone())
		name = wilson+hist_sme[-1].GetName()[6:]
		hist_sme[-1].SetName(name)
		hist_sme[-1].SetTitle(name)

		if g.GetName()=='signal':
                  for k in range(nbinkin):
                    if doSMEkindep:
                      hist_sme_MCstatUp.append(g.Clone())
                      hist_sme_MCstatUp[-1].SetName(name + '_MCstat_binobs'+str(k)+'_sme_'+year+'Up')
                      hist_sme_MCstatDown.append(g.Clone())
                      hist_sme_MCstatDown[-1].SetName(name + '_MCstat_binobs'+str(k)+'_sme_'+year+'Down')
		      for j in range(nbinkin):
                          for i in range(nbintime):
				if j==k:
				    coeff_stat_up = sme_sig_kinbin_directreco_MCstatUp[j].GetBinContent(i+1)
				    coeff_stat_down = sme_sig_kinbin_directreco_MCstatDown[j].GetBinContent(i+1)
				    #coeff_stat_up = sme_sig_kinbin_directreco[j].GetBinContent(i+1) * (1 + sme_sig_kinbin_directreco_MCstatUp[j].GetBinContent(i+1)/sme_sig_kinbin_directreco[j].GetBinContent(i+1))
				    #coeff_stat_down = sme_sig_kinbin_directreco[j].GetBinContent(i+1) * (1 - sme_sig_kinbin_directreco_MCstatDown[j].GetBinContent(i+1)/sme_sig_kinbin_directreco[j].GetBinContent(i+1))
				    #coeff_stat_up = (1+ sme_sig_kinbin_directreco[j].GetBinContent(i+1)) * (1 + sme_sig_kinbin_directreco_MCstatUp[j].GetBinContent(i+1)/sme_sig_kinbin_directreco[j].GetBinContent(i+1)) -1
                                    #coeff_stat_down = (1+ sme_sig_kinbin_directreco[j].GetBinContent(i+1)) * (1 + sme_sig_kinbin_directreco_MCstatDown[j].GetBinContent(i+1)/sme_sig_kinbin_directreco[j].GetBinContent(i+1)) -1
				    #coeff_stat_up = -1 + (1+ cmunu*sme_sig_kinbin_directreco[j].GetBinContent(i+1)) * abs(sme_sig_kinbin_directreco_MCstatUp[j].GetBinContent(i+1)/sme_sig_kinbin_directreco[j].GetBinContent(i+1))
				    #coeff_stat_down =  -1 + (1+ cmunu*sme_sig_kinbin_directreco[j].GetBinContent(i+1)) * abs(sme_sig_kinbin_directreco_MCstatDown[j].GetBinContent(i+1)/sme_sig_kinbin_directreco[j].GetBinContent(i+1))
				else: 
				    coeff_stat_up = sme_sig_kinbin_directreco[j].GetBinContent(i+1)
				    coeff_stat_down = coeff_stat_up
				hist_sme_MCstatUp[-1].SetBinContent(j + i*nbinkin + 1,
						g.GetBinContent(j + i*nbinkin + 1)*(1+cmunu*coeff_stat_up))
				hist_sme_MCstatUp[-1].SetBinError(j + i*nbinkin + 1,
						g.GetBinError(j + i*nbinkin + 1)*(1+cmunu*coeff_stat_up))
				hist_sme_MCstatDown[-1].SetBinContent(j + i*nbinkin + 1,
						g.GetBinContent(j + i*nbinkin + 1)*(1+cmunu*coeff_stat_down))
				hist_sme_MCstatDown[-1].SetBinError(j + i*nbinkin + 1,
						g.GetBinError(j + i*nbinkin + 1)*(1+cmunu*coeff_stat_down))

                for j in range(nbinkin):
                  for i in range(nbintime):
		    if doSMEkindep:
			if doDirectRecoSignal:
                            hist_sme[-1].SetBinContent(j + i*nbinkin + 1,
                                                g.GetBinContent(j + i*nbinkin + 1)*(1+cmunu*sme_sig_kinbin_directreco[j].GetBinContent(i+1)))
                            hist_sme[-1].SetBinError(j + i*nbinkin + 1,
                                                g.GetBinError(j + i*nbinkin + 1)*(1+cmunu*sme_sig_kinbin_directreco[j].GetBinContent(i+1)))
		        elif doResponseMatrix:
			    bincontent_i_j = 0
			    for jgen in range(nbingenkin):
				jreco = j
				propGen = 1
                                if jgen==0:
                                    jgenMat = nbingenkin-1
                                elif jgen==nbingenkin-1:
                                    jgenMat = 0
                                else:
                                    jgenMat = jgen
			        if nputimebin==1:
        			    iputime = 0
	        		else:
            			    iputime = i
				for h2D in histograms_time_responseMatrix_norm[iputime]:
				    #print h2D.GetName()+'  '+g.GetName()[:6]+'_responseMatrix'+g.GetName()[6:]
				    if h2D.GetName()==g.GetName()[:6]+'_responseMatrix'+g.GetName()[6:]:
					#print 'itime='+str(i)+' g='+g.GetName()+' responseMatrix='+h2D.GetName()
					propGen = h2D.GetBinContent(1+jreco, 1+jgenMat)
				if propGen==1:
				    for h2D in histograms_time_responseMatrix_norm[iputime]:
					if h2D.GetName()=='signal_responseMatrix':
	                                    #print 'itime='+str(i)+' g='+g.GetName()+' responseMatrix='+h2D.GetName()
					    propGen = h2D.GetBinContent(1+jreco, 1+jgenMat)
				#print  'itime='+str(i)+' jreco='+str(j)+' jgen='+str(jgen)+' g='+g.GetName()+'  propGen='+str(propGen)
				bincontent_i_j += cmunu*sme_sig_kinbin[jgen].GetBinContent(i+1)*propGen
			    #bincontent_i_j = (1 + bincontent_i_j) * g.GetBinContent(j + i*nbinkin + 1) 
			    #binerror_i_j = (1 + bincontent_i_j) * g.GetBinError(j + i*nbinkin + 1)
			    hist_sme[-1].SetBinContent(j + i*nbinkin + 1, g.GetBinContent(j + i*nbinkin + 1)*(1+bincontent_i_j))
                                                #g.GetBinContent(j + i*nbinkin + 1)*(1+cmunu*sme_sig_kinbin[jgen].GetBinContent(i+1)*propGen))
 			    hist_sme[-1].SetBinError(j + i*nbinkin + 1, g.GetBinError(j + i*nbinkin + 1)*(1+bincontent_i_j))
                                                #g.GetBinError(j + i*nbinkin + 1)*(1+cmunu*sme_sig_kinbin[jgen].GetBinContent(i+1)*propGen))
			else: #No response matrix
			    jgen = j+1 #Modulations SME: j=0: 0 b-jets, j=1: 1 b-jets... (nbingenkin=5: OOA  not used)
			    hist_sme[-1].SetBinContent(j + i*nbinkin + 1,
						g.GetBinContent(j + i*nbinkin + 1)*(1+cmunu*sme_sig_kinbin[jgen].GetBinContent(i+1)))
			    hist_sme[-1].SetBinError(j + i*nbinkin + 1,
						g.GetBinError(j + i*nbinkin + 1)*(1+cmunu*sme_sig_kinbin[jgen].GetBinContent(i+1)))
		    else:
			hist_sme[-1].SetBinContent(j + i*nbinkin + 1,
						g.GetBinContent(j + i*nbinkin + 1)*(1+cmunu*sme_sig.GetBinContent(i+1)))
			hist_sme[-1].SetBinError(j + i*nbinkin + 1,
						g.GetBinError(j + i*nbinkin + 1)*(1+cmunu*sme_sig.GetBinError(i+1)))

	for k in range(nbinkin):
            if doDirectRecoSignal:	
		hist_sme.append(hist_sme_MCstatUp[k])
		hist_sme.append(hist_sme_MCstatDown[k])

		#print(hist_sme[-1].GetName()+' g_Integral='+str(g.Integral())+' hist_sme_integral='+str(hist_sme[-1].Integral()))+' ratio='+str(hist_sme[-1].Integral()/g.Integral()) 


for wilson in wilsonListSingleTop:

	if wilson[2:]=='XX' or wilson[2:]=='XY':
	    cmunu_singletop = cmunu*3
	elif wilson[2:]=='XZ' or wilson[2:]=='YZ':
	    cmunu_singletop = cmunu*10

        sme_singletop = sme_file_singletop.Get('singletop_'+wilson)
        if TString(wilson).Contains('cR'):
            sme_singletop = sme_file_singletop.Get('singletop_cL'+wilson[2:])

        sme_singletop_kinbin = []
        for k in range(nbingenkin):
            print(wilson+"_"+str(k))
            if TString(wilson).Contains('cR'):
                sme_singletop_kinbin.append(sme_file_singletop.Get('singletop_cL'+wilson[2:]+"_"+str(k)))
            else:
                sme_singletop_kinbin.append(sme_file_singletop.Get('singletop_'+wilson+"_"+str(k)))

        for g in hist_unrolled:
	    if g.GetName()=='singletop': #How to treat single top sme decay uncertainty when multiple wilson in the datacard?
		name = g.GetName() + '_sme_decay_'+wilson[2:]
		h_singletop_sme_up = g.Clone()
		h_singletop_sme_up.SetName(name+'Up')
		h_singletop_sme_down = g.Clone()
		h_singletop_sme_down.SetName(name+'Down')
		for i in range(nbintime):
		  for j in range(nbinkin):
		    if doSMEkindep:
                        if doResponseMatrix:
                            bincontent_i_j = 0
                            for jgen in range(nbingenkin):
                                jreco = j
                                propGen = 1
                                if jgen==0:
                                    jgenMat = nbingenkin-1
                                elif jgen==nbingenkin-1:
                                    jgenMat = 0
                                else:
                                    jgenMat = jgen
                                if nputimebin==1:
                                    iputime = 0
                                else:
                                    iputime = i
                                for h2D in histograms_time_responseMatrix[iputime]:
                                    #print h2D.GetName()+'  '+g.GetName()[:6]+'_responseMatrix'+g.GetName()[6:]
                                    if h2D.GetName()==g.GetName()[:6]+'_responseMatrix'+g.GetName()[6:]:
                                        #print 'itime='+str(i)+' g='+g.GetName()+' responseMatrix='+h2D.GetName()
                                        propGen = h2D.GetBinContent(1+jreco, 1+jgenMat)
                                if propGen==1:
                                    for h2D in histograms_time_responseMatrix[iputime]:
                                        if h2D.GetName()=='singletop_responseMatrix':
                                            #print 'itime='+str(i)+' g='+g.GetName()+' responseMatrix='+h2D.GetName()
                                            propGen = h2D.GetBinContent(1+jreco, 1+jgenMat)
                                bincontent_i_j += cmunu_singletop*sme_singletop_kinbin[jgen].GetBinContent(i+1)*propGen
                            h_singletop_sme_up.SetBinContent(j + i*nbinkin + 1, g.GetBinContent(j + i*nbinkin + 1)*(1+bincontent_i_j))
                            h_singletop_sme_down.SetBinContent(j + i*nbinkin + 1, g.GetBinContent(j + i*nbinkin + 1)*(1-bincontent_i_j))
                                #h_singletop_sme_up.SetBinContent(j + i*nbinkin + 1,
                                #                g.GetBinContent(j + i*nbinkin + 1)*(1+cmunu*sme_singletop_kinbin[jgen].GetBinContent(i+1)*propGen))
                                #h_singletop_sme_down.SetBinContent(j + i*nbinkin + 1,
                                #                g.GetBinContent(j + i*nbinkin + 1)*(1-cmunu*sme_singletop_kinbin[jgen].GetBinContent(i+1)*propGen))
			else:
			    h_singletop_sme_up.SetBinContent(j + i*nbinkin + 1,
						g.GetBinContent(j + i*nbinkin + 1 )*(1+cmunu_singletop*sme_singletop_kinbin[jgen+1].GetBinContent(i+1)))
			    h_singletop_sme_down.SetBinContent(j + i*nbinkin + 1,
						g.GetBinContent(j + i*nbinkin + 1 )*(1-cmunu_singletop*sme_singletop_kinbin[jgen+1].GetBinContent(i+1)))
		    else:
			h_singletop_sme_up.SetBinContent(j + i*nbinkin + 1,
						g.GetBinContent(j + i*nbinkin + 1 )*(1+cmunu_singletop*sme_singletop.GetBinContent(i+1)))
			h_singletop_sme_down.SetBinContent(j + i*nbinkin + 1,
						g.GetBinContent(j + i*nbinkin + 1 )*(1-cmunu_singletop*sme_singletop.GetBinContent(i+1)))
		hist_sme.append(h_singletop_sme_up)
		hist_sme.append(h_singletop_sme_down)


for i in range(len(hist_sme)): hist_unrolled.append(hist_sme[i])


#################################################
# Making time-uncorrelated experimental nuisances
#################################################

hist_mc_timeuncorr = []
h_nom = []
ttbar_list_ext = ttbar_list

for wilson in wilsonList:
    ttbar_list_ext.append(wilson)

for h in hist_unrolled:

    for proc in ttbar_list_ext:
        if h.GetName()==proc:
            h_nom.append(h.Clone())

    #Uncorrelated between year
    if TString(h.GetName()).Contains('syst_b') or TString(h.GetName()).Contains('syst_l') or TString(h.GetName()).Contains('syst_em_trig') or TString(h.GetName()).Contains('stat_muon'):
    #if TString(h.GetName()).Contains('syst_b_uncorrelated') or TString(h.GetName()).Contains('syst_l_uncorrelated') or TString(h.GetName()).Contains('syst_em_trig') or TString(h.GetName()).Contains('stat_muon'):
        curname = h.GetName()
        found = curname.find('Up')
        if (found==-1):
            found = curname.find('Down')
        newname = curname[:found] + '_' + year + curname[found:]
        h.SetName(newname)

    #Uncorrelated between processes 
    if TString(h.GetName()).Contains('syst_qcdscale') or TString(h.GetName()).Contains('syst_ps_isr') or TString(h.GetName()).Contains('syst_ps_fsr'):
        curname = h.GetName()
        found = curname.find('Up')
        if (found==-1):
            found = curname.find('Down')
        for h_proc in h_nom:
            proc = h_proc.GetName()
            if TString(curname).Contains(proc):
                if (proc=="ttx" or proc=="vjets" or proc=="singletop" or proc=="dibosons"):
		    newname = curname[:found] + '_' + proc + curname[found:]
		else:
		    newname = curname[:found] + '_signal' + curname[found:]
                h.SetName(newname)

    #if TString(h.GetName()).Contains('syst_qcdscale') or TString(h.GetName()).Contains('syst_ps_isr') or TString(h.GetName()).Contains('syst_ps_fsr') or TString(h.GetName()).Contains('syst_pdfas') or TString(h.GetName()).Contains('syst_pt_top') or TString(h.GetName()).Contains('mtop'):
        #curname = h.GetName()
        #for h_proc in h_nom:
            #proc = h_proc.GetName()
            #if TString(h.GetName()).Contains(proc):  #Was buggy: should normalize per time bin independently (and globally for time syst only)
		#if (doShapeOnly==False): 
		    #applyNominalNorm(h, h_proc)
                #area = h.Integral()
                #h.Scale(h_proc.Integral()/area)

    if doExpTimeNuisance:
	#TString(h.GetName()).Contains('syst_muon_iso')
        if (TString(h.GetName()).Contains('syst_elec_reco') or TString(h.GetName()).Contains('syst_elec_id') or TString(h.GetName()).Contains('syst_muon_id') or TString(h.GetName()).Contains('stat_muon_id') or TString(h.GetName()).Contains('stat_muon_iso') or TString(h.GetName()).Contains('syst_pu') or TString(h.GetName()).Contains('syst_b_correlated') or TString(h.GetName()).Contains('syst_b_uncorrelated') or TString(h.GetName()).Contains('syst_l_correlated') or TString(h.GetName()).Contains('syst_l_uncorrelated') or TString(h.GetName()).Contains('syst_prefiring') or (TString(h.GetName()).Contains('jec') and TString(h.GetName()).Contains('FlavorPureBottom')==False and TString(h.GetName()).Contains('FlavorPureGluon')==False) or TString(h.GetName()).Contains('syst_em_trig') or (puOption!="putime" and TString(h.GetName()).Contains('syst_pu')) or TString(h.GetName()).Contains('emu_trig_stat')):
            curname = h.GetName()
            found = curname.find('Up')
            if (found==-1):
                found = curname.find('Down')
	    for h_proc in h_nom:
		if TString(h.GetName()).Contains(h_proc.GetName()):	
	    	    for i in range(nbintime):
                        hist_mc_timeuncorr.append(h_proc.Clone())
                	newname = curname[:found] + '_t' + str(i) + curname[found:]
                	hist_mc_timeuncorr[-1].SetName(newname)
                	hist_mc_timeuncorr[-1].SetTitle(newname)
			for j in range(nbinkin):
			    hist_mc_timeuncorr[-1].SetBinContent(j + i*nbinkin + 1, h.GetBinContent(j + i*nbinkin + 1))
    #h_new.Scale(hist_triggerSF.GetBinContent(n+1))
    


print 'data = '+str(data_integral)
print 'mc   = '+str(mc_integral)

out = './combine/'+year+'/unrolled/inputs/'
output = TFile(out+observable+'.root', "RECREATE")
hist_data.Write()

for l in hist_unrolled:
    l.Write()
for l in hist_mc_timeuncorr:
    l.Write()
output.Close()

#for i in range(l.GetNbinsX()-400):
#    print str(i)+'#########'+str(hist_mc[0].GetBinContent(i+1))

cmd = 'cp '+out+observable+'.root '+'./combine/'+year+'/sme/inputs/'+observable+'_'+wilson_+'.root'
os.system(cmd)

file = open('./combine/'+year+'/'+observable+'_noe_data.txt','w') 
file.write(str(data_integral)) 
file.close() 
