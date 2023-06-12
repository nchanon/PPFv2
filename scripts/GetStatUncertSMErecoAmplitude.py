

import os
import sys
sys.path.append('./')

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad


f = TFile("../results/2016/flattree/sme_matrices_alt_puinc.root")

wilson = "cL_XX"

doRecoLevelStat = False

if doRecoLevelStat:
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

doGenLevelStat = False

if doGenLevelStat:
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


doPtReweightingCheck = True

if doPtReweightingCheck:
    print 'PT REWEIGHTING CHECK'
    fReweight2016 = TFile("../results/2016/flattree/sme_matrices_alt_puinc_paperV5.root")
    fNoReweight2016 = TFile("../results/2016/flattree/sme_matrices_alt_puinc_Approval.root")
    fReweight2017 = TFile("../results/2017/flattree/sme_matrices_alt_puinc_paperV5.root")
    fNoReweight2017 = TFile("../results/2017/flattree/sme_matrices_alt_puinc_Approval.root")


    for wilson in ["cLXX", "cLXY", "cLXZ", "cLYZ", "cRXX", "cRXY", "cRXZ", "cRYZ", "cXX", "cXY", "cXZ", "cYZ", "dXX", "dXY", "dXZ", "dYZ"]:
	for ib in range(4):

	    hReweight2016 = fReweight2016.Get("reco_2016_central_"+wilson+"_"+str(ib+1))
	    hNoReweight2016 = fNoReweight2016.Get("reco_2016_central_"+wilson+"_"+str(ib+1))
	    diff2016 = (hNoReweight2016.GetMaximum()-hReweight2016.GetMaximum())/hReweight2016.GetMaximum()*100
            hReweight2017 = fReweight2017.Get("reco_2017_central_"+wilson+"_"+str(ib+1))
            hNoReweight2017 = fNoReweight2017.Get("reco_2017_central_"+wilson+"_"+str(ib+1))
            diff2017 = (hNoReweight2017.GetMaximum()-hReweight2017.GetMaximum())/hReweight2017.GetMaximum()*100

	    #print wilson+" n_bjets="+str(ib+1)+" old="+str(hReweight.GetMaximum())+" new="+str(hNoReweight.GetMaximum())+" diff="+str(round(diff,2))+"%"
	    print wilson+" n_bjets="+str(ib+1)+" diff2016="+str(round(diff2016,2))+"% diff2017="+str(round(diff2017,2))+"%"


