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

process = ['signal_dilep', 'signal_semilep', 'signal_hadronic', 'singletop_ST_s', 'singletop_ST_t_top', 'singletop_ST_t_antitop', 'singletop_ST_tW_top', 'singletop_ST_tW_antitop']


for proc in process:

	#print 'inputs/'+year+'/MC/MC_'+proc+'/NtupleProducer/tree.root'	   
        filein = TFile('inputs/'+year+'/MC/MC_'+proc+'/NtupleProducer/tree.root')
        event = filein.Get('events').GetEntriesFast()
        print 'Reco '+proc+' -> '+str(event)

	#print 'inputs/'+year+'/GEN/'+proc+'_NanoGEN_2016_selection_particle.root' 
        filein2 = TFile('inputs/'+year+'/GEN/'+proc+'_NanoGEN_'+year+'_selection_particle.root')
        event2 = filein.Get('events').GetEntriesFast()
        print 'GenMatched '+proc+' -> '+str(event2)

	if event!=event2: print 'Number of events in GenMatched and Reco are different'
