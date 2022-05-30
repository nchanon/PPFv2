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

for i in range(nbintime): 
    hist_data_in.append(data_file.Get('data_obs_bin'+str(i)))
    data_integral += hist_data_in[i].Integral()

# convenient variables
nbinkin = hist_data_in[0].GetNbinsX()
#min_bin = hist_data_in[0].GetXaxis().GetXmin()
#max_bin = hist_data_in[0].GetXaxis().GetXmax()

canvas_data = TCanvas('data_obs', 'data_obs', 1000, 800)
hist_data = TH1F('data_obs','data_obs', nbinkin*nbintime,  0, nbintime)
for iobs in range(nbintime):
    for n in range(nbinkin):
        hist_data.SetBinContent(n + iobs*nbinkin + 1 , hist_data_in[iobs].GetBinContent(n+1))
        hist_data.SetBinError(n + iobs*nbinkin + 1 , hist_data_in[iobs].GetBinError(n+1))

###########
# MC part
###########

mc_file = TFile('./results/'+year+'/flattree/'+observable+'.root')
mc_integral = 0
hist_mc = []

lumisyst_file = TFile('./inputs/timed/LumiUncertainties_'+year+'.root')
hist_lumi_corr = lumisyst_file.Get('hInstLumi_DataScaleFactor')

#lumi_file = TFile('./inputs/timed/AllTimedSyst_'+year+'.root')
triggersyst_file = TFile('./inputs/timed/TriggerSF_'+year+'.root')
hist_triggerSF = triggersyst_file.Get('h_SF_emu_sidereel_Full_UncBand')

index = 0
for l in mc_file.GetListOfKeys():
    if not TString(l.GetName()).Contains('data_obs'):
        hist = mc_file.Get(l.GetName())
        hist_mc.append(TH1F("","", nbinkin*nbintime,  0, nbintime))
        for iobs in range(nbintime):
            for j in range(nbinkin):
                hist_mc[index].SetBinContent(j + iobs*nbinkin + 1 , hist.GetBinContent(j+1)*hist_triggerSF.GetBinContent(iobs+1)) #i => iobs
                hist_mc[index].SetBinError(j + iobs*nbinkin + 1 , hist.GetBinError(iobs+1)) #i => iobs
        del hist
        hist_mc[index].SetName(l.GetName())
        hist_mc[index].SetTitle(l.GetName())
        hist_mc[index].Scale(1./nbintime)
        for g in ttbar_list:
            if TString(l.GetName()) == g:
                mc_integral += hist_mc[index].Integral()
        index += 1
#j + i*nbinkin + 1, hh.GetBinContent(j+1)*hist_triggerSF.GetBinContent(i+1)


#################
# Timed syst part
#################

## timed systematics

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
    for g in hist_mc:
        for b in ttbar_list:
            if g.GetName() == b:
		print(b+' nom='+str(g.Integral()))
                if (l == 'emu_trig' or l=='lumi_stability' or l=='lumi_linearity'): 
		    newprefix = b+'_'+l+'_'+year
		else:
		    newprefix = b+'_'+l
                hist_up = g.Clone()
                hist_down = g.Clone()
                hist_up.SetName(newprefix+'Up')
                hist_up.SetTitle(newprefix+'Up')
                hist_down.SetName(newprefix+'Down')
                hist_down.SetTitle(newprefix+'Down')
                for iobs in range(nbintime):
                    for j in range(nbinkin):
                        if l == 'emu_trig':
                            hist_up.SetBinContent(j + iobs*nbinkin + 1 , 
                                                g.GetBinContent(j + iobs*nbinkin + 1)*(1+hist_triggerSF.GetBinError(iobs+1))) #hist_up => g
                            hist_up.SetBinError(j + iobs*nbinkin + 1 , 
                                                g.GetBinError(j + iobs*nbinkin + 1)*(1+hist_triggerSF.GetBinError(iobs+1)))
                            hist_down.SetBinContent(j + iobs*nbinkin + 1 , 
                                                g.GetBinContent(j + iobs*nbinkin + 1)*(1-hist_triggerSF.GetBinError(iobs+1)))
                            hist_down.SetBinError(j + iobs*nbinkin + 1 , 
                                                g.GetBinError(j + iobs*nbinkin + 1)*(1-hist_triggerSF.GetBinError(iobs+1)))
                        else:
                            hist_up.SetBinContent(j + iobs*nbinkin + 1 , 
                                                g.GetBinContent(j + iobs*nbinkin + 1)*lumi_syst_up[l].GetBinContent(iobs+1))
                            hist_up.SetBinError(j + iobs*nbinkin + 1 , 
                                                g.GetBinError(j + iobs*nbinkin + 1)*lumi_syst_up[l].GetBinError(iobs+1))
                            hist_down.SetBinContent(j + iobs*nbinkin + 1 , 
                                                g.GetBinContent(j + iobs*nbinkin + 1)*lumi_syst_down[l].GetBinContent(iobs+1))
                            hist_down.SetBinError(j + iobs*nbinkin + 1 , 
                                                g.GetBinError(j + iobs*nbinkin + 1)*lumi_syst_down[l].GetBinError(iobs+1))
		print(l+' '+b+' up='+str(hist_up.Integral())+' down='+str(hist_down.Integral()))
                hist_mc.append(hist_up)
                hist_mc.append(hist_down)

###########
# SME part
###########

#cmunu = 0.001
#observable_forSME = observable
#if observable!="m_dilep" and observable!="pt_emu":
#    observable_forSME = "m_dilep" #Will use inclusive SME histos anyways
#sme_file = TFile('./results/'+year+'/flattree/'+observable_forSME+'_sme.root')
#sme_sig = sme_file.Get(wilson)

#sme_sig_kinbin = []
#for k in range(nbinkin+1):
#	print(wilson+"_"+str(k))
#	sme_sig_kinbin.append(sme_file.Get(wilson+"_"+str(k)))

#doSMEkindep = True
#if observable!="m_dilep" and observable!="pt_emu":
#    doSMEkindep = False

#hist_sme = []

#for g in hist_mc:
#    if g.GetName().find('signal') != -1:
#        hist_sme.append(g.Clone())
#        name = wilson+hist_sme[-1].GetName()[6:]
#        hist_sme[-1].SetName(name)
#        hist_sme[-1].SetTitle(name)

#for n in range(len(hist_sme)):
#    for i in range(nbintime):
#        for j in range(nbinkin):
#            hist_sme[n].SetBinContent(j + i*nbinkin + 1, 
#                                        hist_sme[n].GetBinContent(j + i*nbinkin + 1 )*(1+cmunu*sme_sig.GetBinContent(i+1))
#                                    )
#            hist_sme[n].SetBinError(j + i*nbinkin + 1, 
#                                        hist_sme[n].GetBinError(j + i*nbinkin + 1 )*(1+cmunu*sme_sig.GetBinError(i+1))
#                                    )


##############
# JEC,ALT part 
##############

hist_alt_jec = []

jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec.root')
for l in jec_file.GetListOfKeys():
    hh = jec_file.Get(l.GetName())
    hh.Scale(1./nbintime)
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
    h_jec = TH1F("","", nbinkin*nbintime,  0, nbintime)
    h_jec.SetTitle(newname)
    h_jec.SetName(newname)
    for i in range(nbintime):
        for j in range(nbinkin):
            h_jec.SetBinContent(j + i*nbinkin + 1, hh.GetBinContent(j+1)*hist_triggerSF.GetBinContent(i+1)
                                    )
            h_jec.SetBinError(j + i*nbinkin + 1, hh.GetBinError(j + 1))
    hist_alt_jec.append(h_jec)


alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco.root')
for l in alt_file.GetListOfKeys():
    hh = alt_file.Get(l.GetName())
    hh.Scale(1./nbintime)
    hname = l.GetName()
    h_alt = TH1F("", "", nbinkin*nbintime,  0, nbintime)
    for i in range(nbintime):
        for j in range(nbinkin):
            h_alt.SetBinContent(j + i*nbinkin + 1, hh.GetBinContent(j+1)*hist_triggerSF.GetBinContent(i+1)
                                    )
            h_alt.SetBinError(j + i*nbinkin + 1, hh.GetBinError(j + 1))
    h_alt.SetTitle(hname)
    h_alt.SetName(hname)
    hist_alt_jec.append(h_alt)

hist_mc_all = hist_mc
for i in range(len(hist_alt_jec)): hist_mc_all.append(hist_alt_jec[i])
#print hist_mc
#print hist_mc_all


###########
# SME part
###########

cmunu = 0.001
observable_forSME = observable
if observable!="m_dilep" and observable!="pt_emu":
    observable_forSME = "pt_emu" #Will use inclusive SME histos anyways

sme_file = TFile('./results/'+year+'/flattree/'+observable_forSME+'_sme.root')

sme_sig = sme_file.Get(wilson)
sme_singletop = sme_file.Get('singletop_'+wilson)
if TString(wilson).Contains('cR'):
    sme_singletop = sme_file.Get('singletop_cL'+wilson[2:])

sme_sig_kinbin = []
sme_singletop_kinbin = []
for k in range(nbinkin+1):
        print(wilson+"_"+str(k))
        sme_sig_kinbin.append(sme_file.Get(wilson+"_"+str(k)))
	if TString(wilson).Contains('cR'):
	    sme_singletop_kinbin.append(sme_file.Get('singletop_cL'+wilson[2:]+"_"+str(k)))
	else:
	    sme_singletop_kinbin.append(sme_file.Get('singletop_'+wilson+"_"+str(k)))


doSMEkindep = True
if observable!="m_dilep" and observable!="pt_emu":
    doSMEkindep = False

hist_sme = []

for g in hist_mc_all:
    if g.GetName().find('signal') != -1:
        hist_sme.append(g.Clone())
        name = wilson+hist_sme[-1].GetName()[6:]
        hist_sme[-1].SetName(name)
        hist_sme[-1].SetTitle(name)
        for i in range(nbintime):
          for j in range(nbinkin):
	    if doSMEkindep:
                hist_sme[-1].SetBinContent(j + i*nbinkin + 1,
					g.GetBinContent(j + i*nbinkin + 1 )*(1+cmunu*sme_sig_kinbin[j+1].GetBinContent(i+1)))
                hist_sme[-1].SetBinError(j + i*nbinkin + 1,
					g.GetBinError(j + i*nbinkin + 1 )*(1+cmunu*sme_sig_kinbin[j+1].GetBinContent(i+1)))
	    else:
                hist_sme[-1].SetBinContent(j + i*nbinkin + 1,
                                        g.GetBinContent(j + i*nbinkin + 1 )*(1+cmunu*sme_sig.GetBinContent(i+1)))
                hist_sme[-1].SetBinError(j + i*nbinkin + 1,
                                        g.GetBinError(j + i*nbinkin + 1 )*(1+cmunu*sme_sig.GetBinError(i+1)))

    	print(hist_sme[-1].GetName()+' g_Integral='+str(g.Integral())+' hist_sme_integral='+str(hist_sme[-1].Integral()))
    if g.GetName()=='singletop':
	#hist_sme.append(g.Clone())
	name = g.GetName() + '_sme_decay'
	h_singletop_sme_up = g.Clone()
	h_singletop_sme_up.SetName(name+'Up')
	h_singletop_sme_down = g.Clone()
	h_singletop_sme_down.SetName(name+'Down')
	#hist_sme[-1].SetName(name)
        #hist_sme[-1].SetTitle(name)
        for i in range(nbintime):
          for j in range(nbinkin):
            if doSMEkindep:
                h_singletop_sme_up.SetBinContent(j + i*nbinkin + 1,
                                        g.GetBinContent(j + i*nbinkin + 1 )*(1+cmunu*sme_singletop_kinbin[j+1].GetBinContent(i+1)))
                h_singletop_sme_down.SetBinContent(j + i*nbinkin + 1,
                                        g.GetBinContent(j + i*nbinkin + 1 )*(1-cmunu*sme_singletop_kinbin[j+1].GetBinContent(i+1)))
	    else:
                h_singletop_sme_up.SetBinContent(j + i*nbinkin + 1,
                                        g.GetBinContent(j + i*nbinkin + 1 )*(1+cmunu*sme_singletop.GetBinContent(i+1)))
	        h_singletop_sme_down.SetBinContent(j + i*nbinkin + 1,
                                        g.GetBinContent(j + i*nbinkin + 1 )*(1-cmunu*sme_singletop.GetBinContent(i+1)))
	hist_sme.append(h_singletop_sme_up)
	hist_sme.append(h_singletop_sme_down)


for i in range(len(hist_sme)): hist_mc_all.append(hist_sme[i])

#Making the time uncorrelated exp nuisances
hist_mc_timeuncorr = []
h_nom = []
ttbar_list_ext = ttbar_list
ttbar_list_ext.append(wilson)
for h in hist_mc_all:

    for proc in ttbar_list_ext:
        if h.GetName()==proc:
            h_nom.append(h.Clone())

    if TString(h.GetName()).Contains('syst_b_uncorrelated') or TString(h.GetName()).Contains('syst_l_uncorrelated'):
        curname = h.GetName()
        found = curname.find('Up')
        if (found==-1):
            found = curname.find('Down')
        newname = curname[:found] + '_' + year + curname[found:]
        h.SetName(newname)

    if TString(h.GetName()).Contains('syst_qcdscale') or TString(h.GetName()).Contains('syst_ps_isr'):
        curname = h.GetName()
        found = curname.find('Up')
        if (found==-1):
            found = curname.find('Down')
        for h_proc in h_nom:
            proc = h_proc.GetName()
            if TString(curname).Contains(proc):
                if (proc!=wilson):
		    newname = curname[:found] + '_' + proc + curname[found:]
		else:
		    newname = curname[:found] + '_signal' + curname[found:]
                h.SetName(newname)


    if TString(h.GetName()).Contains('syst_qcdscale') or TString(h.GetName()).Contains('syst_ps_isr') or TString(h.GetName()).Contains('syst_ps_fsr') or TString(h.GetName()).Contains('syst_pdfas') or TString(h.GetName()).Contains('syst_pt_top') or TString(h.GetName()).Contains('mtop'):
        curname = h.GetName()
        for h_proc in h_nom:
            proc = h_proc.GetName()
            if TString(h.GetName()).Contains(proc):
                area = h.Integral()
                h.Scale(h_proc.Integral()/area)

    if doExpTimeNuisance:
        if (TString(h.GetName()).Contains('syst_elec_reco') or TString(h.GetName()).Contains('syst_elec_id') or TString(h.GetName()).Contains('syst_muon_id') or TString(h.GetName()).Contains('syst_muon_iso') or TString(h.GetName()).Contains('syst_pu') or TString(h.GetName()).Contains('syst_b_correlated') or TString(h.GetName()).Contains('syst_b_uncorrelated') or TString(h.GetName()).Contains('syst_l_correlated') or TString(h.GetName()).Contains('syst_l_uncorrelated') or TString(h.GetName()).Contains('syst_prefiring') or TString(h.GetName()).Contains('jec')):
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

#for l in hist_sme:
#    l.Write()
for l in hist_mc_all:
    l.Write()
for l in hist_mc_timeuncorr:
    l.Write()
output.Close()

#for i in range(l.GetNbinsX()-400):
#    print str(i)+'#########'+str(hist_mc[0].GetBinContent(i+1))

cmd = 'cp '+out+observable+'.root '+'./combine/'+year+'/sme/inputs/'+observable+'_'+wilson+'.root'
os.system(cmd)

file = open('./combine/'+year+'/'+observable+'_noe_data.txt','w') 
file.write(str(data_integral)) 
file.close() 
