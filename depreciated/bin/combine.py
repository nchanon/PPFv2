#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.directory_manager import *
from tools.sample_manager import *
from tools.style_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack 
from ROOT import TLegend, TApplication, TRatioPlot, TPad

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
observable = args.observable
year = args.year

bin_values = observable_values(observable)[0]

TH1.SetDefaultSumw2(1)

newfile = TFile(results_path(year,'combine','combine_'+observable+year+'.root'),'RECREATE')
newfile.Close()

################################################################################
## function
################################################################################

def buildObservable(name, min,max):
    canvas = TCanvas(name, name)
    hist = TH1F("","", bin_values[0], bin_values[1], bin_values[2])
    rootfile = []
    for i in range(min,max):
        namefile = sample_list['MC'][year][i]+'_'+observable+'.root'
        rootfile.append(TFile(results_path(year,'TH1/MC',namefile)))
        foo = rootfile[i-min].Get(observable)
        hist.Add(foo)
    del rootfile
    rootfile = TFile(results_path(year,'combine','combine_'+observable+year+'.root'), 'UPDATE')
    hist.Draw()
    hist.SetName(name)
    hist.Write()
    rootfile.Close()

def buildSystematics(name, min,max, syst, syst_index):
    canvas = TCanvas(name+systematic_list[syst_index]+syst, 
                     name+systematic_list[syst_index]+syst)
    hist = TH1F("","", bin_values[0], bin_values[1], bin_values[2])
    rootfile = []
    for i in range(min,max):
        namefile = sample_list['MC'][year][i]+'_'+observable+'_'+systematic_list[syst_index]+syst+'.root'
        rootfile.append(TFile(results_path(year,'TH1/SYST',namefile)))
        foo = rootfile[i-min].Get(observable)
        hist.Add(foo)
    del rootfile
    rootfile =  TFile(results_path(year,'combine','combine_'+observable+year+'.root'), 'UPDATE')
    hist.Draw()
    hist.SetName(name+systematic_list[syst_index]+syst)
    hist.Write()
    rootfile.Close()

###################
#signal
###################

canvas_data = TCanvas('data_obs', 'data_obs', 1000, 800)
hist_data = TH1F('data_obs','data_obs', bin_values[0], bin_values[1], bin_values[2])
rootfile = []
for i in range(len(sample_list['DATA'][year])):
    namefile = sample_list['DATA'][year][i]+'_'+observable+'.root'
    rootfile.append(TFile(results_path(year,'TH1/DATA',namefile)))
    foo = rootfile[i].Get(observable)
    hist_data.Add(foo)
del rootfile
rootfile =  TFile(results_path(year,'combine','combine_'+observable+year+'.root'),  'UPDATE')
hist_data.Draw()
hist_data.Write()
rootfile.Close()




buildObservable('sig',      0,3)
buildObservable('TTX',      3,11)
buildObservable('ST',       11,16)
buildObservable('dibosons', 16,18)
buildObservable('Wjets',    18,20)
buildObservable('Zjets',    20,23)

for i in range(len(systematic_list)):
    buildSystematics('sig',      0,3,   'Up',i)
    buildSystematics('sig',      0,3,   'Down',i)
    buildSystematics('TTX',      3,11,   'Up',i)
    buildSystematics('TTX',      3,11,   'Down',i)
    buildSystematics('ST',       11,16,   'Up',i) 
    buildSystematics('ST',       11,16,   'Down',i) 
    buildSystematics('dibosons', 16,18,  'Up',i)
    buildSystematics('dibosons', 16,18,  'Down',i)
    buildSystematics('Wjets',    18,20, 'Up',i)
    buildSystematics('Wjets',    18,20, 'Down',i)
    buildSystematics('Zjets',    20,23, 'Up',i)
    buildSystematics('Zjets',    20,23, 'Down',i)


###################
###################