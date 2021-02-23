import os, sys
import argparse 
import numpy as np
sys.path.append('../../')

from ROOT import TH1F, TFile, TCanvas, TLegend

from tools.style_manager import *

###################
## Initialisation
###################

parser = argparse.ArgumentParser()
parser.add_argument('inputname', help='display your input file')
parser.add_argument('observable', help='display your observable')

args = parser.parse_args()
inputname = args.inputname
observable = args.observable

inputfile = TFile(inputname)

###################
## Inputs
###################

#observable = raw_input('Choose your observable : ')
#print observable


###################
## Body
###################

histogram = inputfile.Get(observable)

if observable.find('m_dilep') != -1:
    title = 'dilepton mass' 
elif  observable.find('n_bjets') != -1:
    title = 'b-jets multiplicity'
if inputname.find('_2017') != -1:
    header = 'L = 41.53 fb^{-1} (13TeV)'
elif inputname.find('_2016') != -1:
    header = '35.9 fb^{-1} (13TeV)'
if observable.find('_uncorrelated') != -1:
    leg = 'sig uncorrelated / bkgd correlated'
elif observable.find('_correlated') != -1:
    leg = 'sig / bkgd correlated'



# cosmetic
histogram.SetName('')
histogram.SetTitle(title)

histogram.SetStats(0)
style_labels_counting(histogram, 'r (signal strengh)', 'sidereal time (h)')


legend = TLegend()
legend.SetHeader(header, 'C')
legend.AddEntry(histogram, leg, 'lep')

canvas = TCanvas('observable', observable, 1000, 600)
histogram.Draw()
legend.Draw('SAME')

raw_input('exit')

canvas.SaveAs('./out/'+observable+'.png')