import os, sys
import argparse

sys.path.append('./')

from ROOT import TFile, TH1, TH2, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame
from ROOT import TLine,TGraphAsymmErrors
from ROOT import TStyle, gStyle, TColor, TLatex

#from collections import OrderedDict

import math

from tools.style_manager import *
from tools.sample_manager import *
import array

import tools.tdrstyle as tdr
tdr.setTDRStyle()

###################
## Initialisation
###################

parser = argparse.ArgumentParser()
#parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
#parser.add_argument('asimov',nargs='?', help='set if asimov test', default='')
#parser.add_argument('title', help='display your observable title')

args = parser.parse_args()
#observable = args.observable
year = args.year

list_lumi = []

fileDATA = TFile('inputs/'+year+'/DATA/NtupleProducer_allData'+year+'.root')
tree = fileDATA.Get('events')
nEv=tree.GetEntriesFast()
#nEv=100
for ientry in range(nEv):
    if ientry %  10000 == 0:
	print ientry
    tree.GetEvent(ientry)
    #Apply selection
    if tree.n_jets<2 and tree.n_bjets<1:
	continue
    #Trigger needs to be applied
    keep=False
    for trig_path in triggers[year]:
	if tree.GetLeaf(trig_path).GetValue()==1:
	    keep=True
    #print keep
    if keep==False:
	continue
    #print('RunNumber='+str(tree.run)+' lumisection='+str(tree.lumi))
    list_lumi.append([tree.run, tree.lumi])

#print('list_lumi '+str(len(list_lumi))+' elements: '+str(list_lumi))
#print list_lumi[0]
#print list_lumi[5][1]
#list_lumi_sorted = sorted(list_lumi,key=lambda l:l[0])
#print list_lumi_sorted

set_lumi = set(tuple(element) for element in list_lumi)
list_lumi_duplicatesrm = [list(t) for t in set_lumi]
#list_lumi_duplicatesrm = list(set(list_lumi))
#list_lumi_duplicatesrm = list(OrderedDict.fromkeys(list_lumi))
#print list_lumi_duplicatesrm

list_lumi_sorted = sorted(list_lumi_duplicatesrm,key=lambda l:(l[0],l[1]))
#print ('list_lumi_sorted '+str(len(list_lumi_sorted) )+' elements: '+str(list_lumi_sorted))

print('list_lumi '+str(len(list_lumi))+' elements, list_lumi_sorted '+str(len(list_lumi_sorted) )+' elements')

