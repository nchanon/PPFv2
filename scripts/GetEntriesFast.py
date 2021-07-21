#!/usr/bin/env python

import os
import sys
sys.path.append('./')

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad

nbin = 24

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('directory', help='display your directory')

args = parser.parse_args()
directory = args.directory

TH1.SetDefaultSumw2(1)

listdir = os.listdir(directory)
listdir.sort()
for l in listdir:
    try :
        filein = TFile(directory+'/'+l+'/tree.root')
        event = filein.Get('events').GetEntriesFast()
        print l+' -> '+str(event)
    except:
        pass