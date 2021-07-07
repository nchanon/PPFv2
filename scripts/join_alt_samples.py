import os 

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
observable = args.observable
year = args.year


namefile = './results/'+year+'/flattree/'+observable+'_alt.root'


infile = TFile(name, 'UPDATE')
for l in rootfile_input.GetListOfKeys():
