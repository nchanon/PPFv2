#!/usr/bin/env python

import os
import sys
sys.path.append('./')

from tools.directory_manager import *
from tools.sample_manager import *
from tools.style_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack 
from ROOT import TLegend, TApplication, TRatioPlot, TPad

nbin = int(raw_input("number of bin : "))


################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
observable = args.observable
year = args.year

binning = observable_values(observable)[0]


TH1.SetDefaultSumw2(1)

directory = 'unrolled'

################################################################################
## function
################################################################################

def card_generator(filename, observed):
    foo =  'imax 1 number of bins\n'
    foo += 'jmax 5 number of background processes\n'
    foo += 'kmax * number of nuisance parameters\n'
    foo += '-------------------------------------------------------------\n'
    foo += 'shapes * * ./inputs/'+filename+' $PROCESS $PROCESS_$SYSTEMATIC\n'
    foo += 'shapes sig * ./inputs/'+filename+' $PROCESS $PROCESS_$SYSTEMATIC\n'
    foo += '--------------------------------------------------------------\n'
    foo += 'bin '+str(observable)+'\n'
    foo += 'observation '+str(observed)+'\n'
    foo += '--------------------------------------------------------------------------------------------- \n'
    foo += 'bin                          '+str(observable)+'   '+str(observable)+'   '+str(observable)+'     '+str(observable)+'     '+str(observable)+'     '+str(observable)+'  \n'
    foo += 'process                      signal    ttx       singletop   dibosons    wjets       zjets    \n'
    foo += 'process                      0         1         2           3           4           5        \n'
    foo += 'rate                        -1        -1        -1          -1          -1          -1        \n'
    foo += '--------------------------------------------------------------------------------------------- \n'
    foo += 'syst_muon_iso          shape 1         1         1           1           1           1        \n'
    foo += 'syst_muon_id           shape 1         1         1           1           1           1        \n'
    foo += 'syst_em_trig           shape 1         1         1           1           1           1        \n'
    foo += 'syst_elec_id           shape 1         1         1           1           1           1        \n'
    foo += 'syst_elec_reco         shape 1         1         1           1           1           1        \n'
    foo += 'lumi                   lnN   1.023     1.023     1.023       1.023       1.023       1.023    \n'
    return foo

################################################################################
## loop
################################################################################

output_name = results_path(year,'combine/'+directory,'combine_'+observable+'_'+str(nbin)+'_unrolled')

newfile = TFile(output_name+'.root','RECREATE')
newfile.Close()

canvas_data = TCanvas('data_obs', 'data_obs', 1000, 800)
hist_data = TH1F('data_obs','data_obs', binning[0]*nbin,  0, nbin)


rootfile = []
for ibin in range(nbin):
    namefile = 'combine_'+observable+'_'+str(nbin)+'_'+str(ibin)+'.root'
    rootfile.append(TFile(results_path(year,'combine/one_bin',namefile)))
    foo = rootfile[ibin].Get('data_obs')
    for j in range(binning[0]):
        hist_data.SetBinContent(j + ibin*binning[0] + 1 , foo.GetBinContent(j))
del rootfile
n_observed = hist_data.Integral()
print 'data = '+str(hist_data.Integral())
rootfile =  TFile(output_name+'.root',  'UPDATE')
hist_data.Draw()
hist_data.Write()
rootfile.Close()


file = open(output_name+'_datacard.txt','w') 
file.write(card_generator('combine_'+observable+'_'+str(nbin)+'_unrolled.root',n_observed)) 
file.close() 

# MC 
for l in ttbar_list:
    canvas = TCanvas(l, l)
    hist = TH1F(l,l, binning[0]*nbin,  0, nbin)
    rootfile = []
    for ibin in range(nbin):
        rootfile.append(TFile(results_path(year,'combine/one_bin',namefile)))
        foo = rootfile[ibin].Get(l)
        for j in range(binning[0]):
            hist.SetBinContent(j + ibin*binning[0] + 1 , foo.GetBinContent(j))
    del rootfile
    rootfile = TFile(output_name+'.root', 'UPDATE')
    print l+' '+str(hist.Integral())
    hist.SetName(l)
    hist.Draw()
    hist.Write()
    rootfile.Close()


up_down = ['Up', 'Down']

#SYST 
for i, s in enumerate(up_down):
    for syst in systematic_list:
        for l in ttbar_list:
            rootfile = []
            canvas = TCanvas(l+'_'+syst+s, l+'_'+syst+s)
            hist = TH1F(l+'_'+syst+s, l+'_'+syst+s, binning[0]*nbin,  0, nbin)
            for ibin in range(nbin):
                rootfile.append(TFile(results_path(year,'combine/one_bin',namefile)))
                foo = rootfile[ibin].Get(l+'_'+syst+s)
                for j in range(binning[0]):
                    hist.SetBinContent(j + ibin*binning[0] + 1 , foo.GetBinContent(j))
            del rootfile
            rootfile =  TFile(output_name+'.root', 'UPDATE')
            hist.SetName(l+'_'+syst+s)
            hist.Draw()
            hist.Write()
            rootfile.Close()
