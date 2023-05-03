#!/usr/bin/env python

import os
import sys
sys.path.append('./')

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString, TTree
from ROOT import TLegend, TApplication, TRatioPlot, TPad

nbin = 24

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('year', help='display your year')

args = parser.parse_args()
year = args.year

process = ['signal_dilep']
#process = ['signal_dilep', 'signal_semilep', 'signal_hadronic', 'singletop_ST_s', 'singletop_ST_t_top', 'singletop_ST_t_antitop', 'singletop_ST_tW_top', 'singletop_ST_tW_antitop']

#(n_jets>=2 && n_bjets>=1)*weight_puinc*weight_generator*weight_top*weight_prefiring*weight_sfe_id*weight_sfe_reco*weight_sfm_id*weight_sfm_iso*weight_sfb*weight_sf_em_trig*(1||(trg_muon_electron_mu23ele12_fired[0]==1)||(trg_muon_electron_mu8ele23_fired[0]==1)||(trg_muon_mu24_fired[0]==1)||(trg_muon_mutk24_fired[0]==1)||(trg_electron_ele27_fired[0]==1)||(trg_muon_electron_mu8ele23DZ_fired[0]==1)||(trg_muon_electron_mu23ele12DZ_fired[0]==1))

string_cut = "n_jets>=2 && n_bjets>=1"
string_weight_pu = ["weight_pu", "weight_punew", "weight_puinc"]
for i in range(24):
    string_weight_pu.append("weight_putime"+str(i))
print string_weight_pu

string_weight_notrigSF = "weight_generator*weight_top*weight_prefiring*weight_sfe_id*weight_sfe_reco*weight_sfm_id*weight_sfm_iso*weight_sfb"
string_trigger_accept = "(trg_muon_electron_mu23ele12_fired[0]==1)||(trg_muon_electron_mu8ele23_fired[0]==1)||(trg_muon_mu24_fired[0]==1)||(trg_muon_mutk24_fired[0]==1)||(trg_electron_ele27_fired[0]==1)||(trg_muon_electron_mu8ele23DZ_fired[0]==1)||(trg_muon_electron_mu23ele12DZ_fired[0]==1)"

string_weight_trig = "("+string_cut+")*("+string_weight_notrigSF+")*("+string_trigger_accept+")"
string_weight_notrig = "("+string_cut+")*("+string_weight_notrigSF+")"

effTrig_pu = []
for proc in process:

	print 'inputs/'+year+'/MC_save/MC_'+proc+'/NtupleProducer/tree.root'	   
        filein = TFile('inputs/'+year+'/MC_save/MC_'+proc+'_noHLT/NtupleProducer/tree.root')
	tree = filein.Get("events")

	for k in range(len(string_weight_pu)):
	    h = TH1F("h","h",10,0,10)
	    tree.Draw("n_bjets >> h","("+string_weight_notrig+")*"+string_weight_pu[k])
	    print h.Integral()

	    hTrig = TH1F("hTrig","hTrig",10,0,10)
            tree.Draw("n_bjets >> hTrig","("+string_weight_trig+")*"+string_weight_pu[k])
	    print hTrig.Integral()

	    print string_weight_pu[k]+" : MC trig eff="+str(hTrig.Integral()/h.Integral())
	    effTrig_pu.append(hTrig.Integral()/h.Integral())

print 'Summary:'
for k in range(len(string_weight_pu)):
    print string_weight_pu[k]+" : MC trig eff="+str(effTrig_pu[k])

raw_input()
exit()

        #event = filein.Get('events').GetEntriesFast()
        #print 'Reco '+proc+' -> '+str(event)

	#print 'inputs/'+year+'/GEN/'+proc+'_NanoGEN_2016_selection_particle.root' 
        #filein2 = TFile('inputs/'+year+'/GEN/'+proc+'_NanoGEN_'+year+'_selection_particle.root')
        #event2 = filein.Get('events').GetEntriesFast()
        #print 'GenMatched '+proc+' -> '+str(event2)

	#if event!=event2: print 'Number of events in GenMatched and Reco are different'
