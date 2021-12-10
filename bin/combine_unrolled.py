#!/usr/bin/env python

import os
import sys
sys.path.append('./')


import argparse

from tools.sample_manager import *
from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

nbintime = 24

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
nbinmass = hist_data_in[0].GetNbinsX()
#min_bin = hist_data_in[0].GetXaxis().GetXmin()
#max_bin = hist_data_in[0].GetXaxis().GetXmax()

canvas_data = TCanvas('data_obs', 'data_obs', 1000, 800)
hist_data = TH1F('data_obs','data_obs', nbinmass*nbintime,  0, nbintime)
for iobs in range(nbintime):
    for n in range(nbinmass):
        hist_data.SetBinContent(n + iobs*nbinmass + 1 , hist_data_in[iobs].GetBinContent(n+1))
        hist_data.SetBinError(n + iobs*nbinmass + 1 , hist_data_in[iobs].GetBinError(n+1))

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
        hist_mc.append(TH1F("","", nbinmass*nbintime,  0, nbintime))
        for iobs in range(nbintime):
            for j in range(nbinmass):
                hist_mc[index].SetBinContent(j + iobs*nbinmass + 1 , hist.GetBinContent(j+1)*hist_weight.GetBinContent(iobs+1)) #i => iobs
                hist_mc[index].SetBinError(j + iobs*nbinmass + 1 , hist.GetBinError(iobs+1)) #i => iobs
        del hist
        hist_mc[index].SetName(l.GetName())
        hist_mc[index].SetTitle(l.GetName())
        hist_mc[index].Scale(1./nbintime)
        for g in ttbar_list:
            if TString(l.GetName()) == g:
                mc_integral += hist_mc[index].Integral()
        index += 1
#j + i*nbinmass + 1, hh.GetBinContent(j+1)*hist_weight.GetBinContent(i+1)

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
                    for j in range(nbinmass):
                        if l == 'emu_trig':
                            hist_up.SetBinContent(j + iobs*nbinmass + 1 , 
                                                g.GetBinContent(j + iobs*nbinmass + 1)*(1+lumi_syst_up[l].GetBinError(iobs+1))) #hist_up => g
                            hist_up.SetBinError(j + iobs*nbinmass + 1 , 
                                                g.GetBinError(j + iobs*nbinmass + 1)*(lumi_syst_up[l].GetBinError(iobs+1)))
                            hist_down.SetBinContent(j + iobs*nbinmass + 1 , 
                                                g.GetBinContent(j + iobs*nbinmass + 1)*(1-lumi_syst_up[l].GetBinError(iobs+1)))
                            hist_down.SetBinError(j + iobs*nbinmass + 1 , 
                                                g.GetBinError(j + iobs*nbinmass + 1)*(lumi_syst_up[l].GetBinError(iobs+1)))
                        else:
                            hist_up.SetBinContent(j + iobs*nbinmass + 1 , 
                                                g.GetBinContent(j + iobs*nbinmass + 1)*lumi_syst_up[l].GetBinContent(iobs+1))
                            hist_up.SetBinError(j + iobs*nbinmass + 1 , 
                                                g.GetBinError(j + iobs*nbinmass + 1)*lumi_syst_up[l].GetBinError(iobs+1))
                            hist_down.SetBinContent(j + iobs*nbinmass + 1 , 
                                                g.GetBinContent(j + iobs*nbinmass + 1)*lumi_syst_down[l].GetBinContent(iobs+1))
                            hist_down.SetBinError(j + iobs*nbinmass + 1 , 
                                                g.GetBinError(j + iobs*nbinmass + 1)*lumi_syst_down[l].GetBinError(iobs+1))
		print(l+' '+b+' up='+str(hist_up.Integral())+' down='+str(hist_down.Integral()))
                hist_mc.append(hist_up)
                hist_mc.append(hist_down)

###########
# SME part
###########

cmunu = 0.001
sme_file = TFile('./results/'+year+'/flattree/sme.root')
sme_sig = sme_file.Get(wilson)

sme_sig_massbin = []
for k in range(nbinmass+1):
	print(wilson+"_"+str(k))
	sme_sig_massbin.append(sme_file.Get(wilson+"_"+str(k)))

#hist_sme = []

#for g in hist_mc:
#    if g.GetName().find('signal') != -1:
#        hist_sme.append(g.Clone())
#        name = wilson+hist_sme[-1].GetName()[6:]
#        hist_sme[-1].SetName(name)
#        hist_sme[-1].SetTitle(name)

#for n in range(len(hist_sme)):
#    for i in range(nbintime):
#        for j in range(nbinmass):
#            hist_sme[n].SetBinContent(j + i*nbinmass + 1, 
#                                        hist_sme[n].GetBinContent(j + i*nbinmass + 1 )*(1+cmunu*sme_sig.GetBinContent(i+1))
#                                    )
#            hist_sme[n].SetBinError(j + i*nbinmass + 1, 
#                                        hist_sme[n].GetBinError(j + i*nbinmass + 1 )*(1+cmunu*sme_sig.GetBinError(i+1))
#                                    )


hist_alt_jec = []

jec_file = TFile('./results/'+year+'/flattree/'+observable+'_jec.root')
for l in jec_file.GetListOfKeys():
    hh = jec_file.Get(l.GetName())
    hh.Scale(1./nbintime)
    hname = l.GetName() 
    if(hname.find('TotalUp')!= -1):
        name = hname[:-7]+'jecUp'
    elif(hname.find('TotalDown')!= -1):
        name = hname[:-9]+'jecDown'
    h_jec = TH1F("","", nbinmass*nbintime,  0, nbintime)
    h_jec.SetTitle(name)
    h_jec.SetName(name)
    for i in range(nbintime):
        for j in range(nbinmass):
            h_jec.SetBinContent(j + i*nbinmass + 1, hh.GetBinContent(j+1)*hist_weight.GetBinContent(i+1)
                                    )
            h_jec.SetBinError(j + i*nbinmass + 1, hh.GetBinError(j + 1))
    hist_alt_jec.append(h_jec)


alt_file = TFile('./results/'+year+'/flattree/'+observable+'_color_reco.root')
for l in alt_file.GetListOfKeys():
    hh = alt_file.Get(l.GetName())
    hh.Scale(1./nbintime)
    hname = l.GetName()
    h_alt = TH1F("", "", nbinmass*nbintime,  0, nbintime)
    for i in range(nbintime):
        for j in range(nbinmass):
            h_alt.SetBinContent(j + i*nbinmass + 1, hh.GetBinContent(j+1)*hist_weight.GetBinContent(i+1)
                                    )
            h_alt.SetBinError(j + i*nbinmass + 1, hh.GetBinError(j + 1))
    h_alt.SetTitle(hname)
    h_alt.SetName(hname)
    hist_alt_jec.append(h_alt)


hist_sme = []

hist_mc_all = hist_mc
for i in range(len(hist_alt_jec)): hist_mc_all.append(hist_alt_jec[i])
#print hist_mc
#print hist_mc_all

for g in hist_mc_all:
    if g.GetName().find('signal') != -1:
        hist_sme.append(g.Clone())
        name = wilson+hist_sme[-1].GetName()[6:]
        hist_sme[-1].SetName(name)
        hist_sme[-1].SetTitle(name)
#for n in range(len(hist_sme)):
        for i in range(nbintime):
          for j in range(nbinmass):
            hist_sme[-1].SetBinContent(j + i*nbinmass + 1,
					g.GetBinContent(j + i*nbinmass + 1 )*(1+cmunu*sme_sig_massbin[j+1].GetBinContent(i+1))
                                        #g.GetBinContent(j + i*nbinmass + 1 )*(1+cmunu*sme_sig.GetBinContent(i+1))
                                    )
            hist_sme[-1].SetBinError(j + i*nbinmass + 1,
                                        #g.GetBinError(j + i*nbinmass + 1 )*(1+cmunu*sme_sig.GetBinError(i+1))
					g.GetBinError(j + i*nbinmass + 1 )*(1+cmunu*sme_sig_massbin[j+1].GetBinContent(i+1))
                                    )
    	print(hist_sme[-1].GetName()+' g_Integral='+str(g.Integral())+' hist_sme_integral='+str(hist_sme[-1].Integral()))




print 'data = '+str(data_integral)
print 'mc   = '+str(mc_integral)

out = './combine/'+year+'/unrolled/inputs/'
output = TFile(out+observable+'.root', "RECREATE")
hist_data.Write()
for l in hist_sme:
    l.Write()
for l in hist_mc_all:
    l.Write()
#for l in hist_mc:
#    l.Write()
#for l in hist_alt_jec:
#    l.Write()
output.Close()

#for i in range(l.GetNbinsX()-400):
#    print str(i)+'#########'+str(hist_mc[0].GetBinContent(i+1))

cmd = 'cp '+out+observable+'.root '+'./combine/'+year+'/sme/inputs/'+observable+'_'+wilson+'.root'
os.system(cmd)

file = open('./combine/'+year+'/'+observable+'_noe_data.txt','w') 
file.write(str(data_integral)) 
file.close() 
