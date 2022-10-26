#!/usr/bin/env python

import os
import sys
sys.path.append('./')

import argparse

from tools.sample_manager import *
from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('analysis', help='inclusive, differential or sme')
parser.add_argument('timebin', help='display the time bin')

args = parser.parse_args()
observable = args.observable
year = args.year
analysis = args.analysis
timebin = int(args.timebin)

TH1.SetDefaultSumw2(1)

stimebin="";
if (timebin==-1):
     stimebin = "_puold";
if (timebin==-2):
     stimebin = "_punew";
if (timebin==-3):
     stimebin = "_puinc";
if (timebin>=0):
     stimebin = "_put"+str(timebin);

nbin = 24

################################################################################
## Check combine inputs
################################################################################

fInput = []
ttbar_list_ext = ttbar_list
wilsonList = ["cLXX","cLXY","cLXZ","cLYZ","cRXX","cRXY","cRXZ","cRYZ","cXX","cXY","cXZ","cYZ","dXX","dXY","dXZ","dYZ"]

if analysis=='inclusive':
    fInput.append(TFile('./combine/'+year+'/inclusive/inputs/'+observable+'_inclusive'+stimebin+'.root'))
    itimebin = 0
if analysis=='differential':
    for n in range(nbin):
	fInput.append(TFile('./combine/'+year+'/one_bin/inputs/'+observable+'_'+str(nbin)+'_'+str(n)+'.root'))
    itimebin = timebin
if analysis=='sme':
    fInput.append(TFile('./combine/'+year+'/sme/inputs/'+observable+'_sme_all.root'))
    itimebin = 0
    for wilson in wilsonList:
        ttbar_list_ext.append(wilson)

#Get nominal histos
h_nom = []
for l in fInput[itimebin].GetListOfKeys():
    for proc in ttbar_list_ext:
	if proc==l.GetName():
	    h_nom.append(fInput[itimebin].Get(l.GetName()))

#Get list of nuisances
list_nuis = []
for l in fInput[itimebin].GetListOfKeys():

    isNuis=True
    if l.GetName()=='data_obs': 
	isNuis=False
    for proc in ttbar_list_ext:
	if l.GetName()==proc:
	    isNuis=False

    if isNuis==True:
	#print l.GetName()
	alreadyCaught=False
	for nuis in list_nuis:
	    if l.GetName().find(nuis)!=-1:
		alreadyCaught=True
	if alreadyCaught==False:
	    found = l.GetName().rfind('Up')
	    if (found==-1):
		found = l.GetName().find('Down')
	    #print l.GetName()[:found]
	    list_nuis.append(l.GetName()[:found])

print list_nuis


#Get central, up, down histos
h_central = []
h_up = []
h_down = []
for nuis in list_nuis:
    #print nuis
    for proc in ttbar_list_ext:
	if nuis.find(proc)!=-1:
	    h_central.append(fInput[itimebin].Get(proc))
    h_up.append(fInput[itimebin].Get(nuis+"Up"))
    h_down.append(fInput[itimebin].Get(nuis+"Down"))
    #print h_central[-1].GetName()+' '+h_up[-1].GetName()+' '+h_down[-1].GetName()
    if (abs(h_up[-1].Integral()/h_central[-1].Integral()-1)>0.0001 or abs(h_down[-1].Integral()/h_central[-1].Integral()-1)>0.0001):
        print 'Ratio '+ nuis+' up='+str(h_up[-1].Integral()/h_central[-1].Integral())+' down='+str(h_down[-1].Integral()/h_central[-1].Integral())



