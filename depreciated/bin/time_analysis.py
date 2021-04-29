#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.directory_manager import *
from tools.sample_manager import *
from tools.style_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TMath
from ROOT import TLegend, TApplication, TRatioPlot, TPad

################################################################################
## Initialisation stuff
################################################################################

nbin_time = 12

ntemoin = 20

temoin = TH1F("essaie", "essaie", ntemoin, 0, 10)
for i in range(ntemoin):
    temoin.SetBinContent(i+1,TMath.Exp(i/3.14))


time_histo = TH1F("time", "time", nbin_time*ntemoin, 0, nbin_time)
for i in range(nbin_time):
    for j in range(ntemoin):
        time_histo.SetBinContent(ntemoin*i + j+1, temoin.GetBinContent(j+1))
time_histo.Draw()

exit = raw_input("Press key to quit : ") 
