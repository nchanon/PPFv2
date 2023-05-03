

import os
import sys
sys.path.append('./')

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad


f = TFile("../results/2016/flattree/sme_matrices_alt_puinc.root")

wilson = "cL_XX"

print 'RECO LEVEL'
for ib in range(5):
    #print str(ib)+": "+"signal_dilep_LO_amplitude_"+wilson+"_n_bjets_"+str(ib)+"_t0" 
    if ib>0:
	#h = f.Get("signal_dilep_LO_amplitude_d_YZ_n_bjets_4_t4")	
        h = f.Get("signal_dilep_LO_amplitude_"+wilson+"_n_bjets_"+str(ib)+"_t0")
        Mean = h.GetMean()
        MeanError = h.GetMeanError()
        relMeanError = abs(h.GetMeanError()/(h.GetMean()-1))
        print wilson+" n_bjets="+str(ib)+ " relMeanError="+str(round(relMeanError*100,2))+"%"

print 'GEN LEVEL'

fGen = TFile('../inputs/pheno/signal_miniAOD_LO_Comb_n_bjets_13TeVCMSnanoGEN_particle.root')

wilson = "cLXX"

for ib in range(6):
	hGen = fGen.Get("signal_miniAOD_LO_Comb_"+wilson+"_"+str(ib))
	hGenDown = fGen.Get("signal_miniAOD_LO_Comb_"+wilson+"_"+str(ib)+"_MCstatDown")
	hGenUp = fGen.Get("signal_miniAOD_LO_Comb_"+wilson+"_"+str(ib)+"_MCstatUp")
	Mean = hGen.GetBinContent(1)
	relMeanError = (hGenUp.GetBinContent(1)-Mean)/Mean
	if ib<5:
	    print wilson+" n_bjets="+str(ib)+ " relMeanError="+str(round(relMeanError*100,2))+"%"
	else:
	    print wilson+" OOA relMeanError="+str(round(relMeanError*100,2))+"%"
