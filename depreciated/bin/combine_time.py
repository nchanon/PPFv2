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

directory = 'one_bin'

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
    foo += 'syst_pu                shape 1         1         1           1           1           1        \n'
    foo += 'lumi                   lnN   1.023     1.023     1.023       1.023       1.023       1.023    \n'
    return foo

################################################################################
## loop
################################################################################

for ibin in range(nbin):

    n_observed = 0

    output_name = results_path(year,'combine/'+directory,'combine_'+observable+'_'+str(nbin)+'_'+str(ibin))

    newfile = TFile(output_name+'.root','RECREATE')
    newfile.Close()

    # DATA
    canvas_data = TCanvas('data_obs', 'data_obs', 1000, 800)
    hist_data = TH1F('data_obs','data_obs', binning[0], binning[1], binning[2])
    rootfile = []
    for i in range(len(sample_list['DATA'][year])):
        namefile = sample_list['DATA'][year][i]+'_'+observable+'_time'+str(nbin)+'_'+str(ibin)+'.root'
        rootfile.append(TFile(results_path(year,'TH1/DATA/time',namefile)))
        foo = rootfile[i].Get(observable+'_time'+str(nbin)+'_'+str(ibin))
        hist_data.Add(foo)
    del rootfile
    n_observed = hist_data.Integral()
    print 'data = '+str(hist_data.Integral())
    rootfile =  TFile(output_name+'.root',  'UPDATE')
    hist_data.Draw()
    hist_data.Write()
    rootfile.Close()

    file = open(output_name+'_datacard.txt','w') 
    file.write(card_generator('combine_'+observable+'_'+str(nbin)+'_'+str(ibin)+'.root',n_observed)) 
    file.close() 

    # MC 
    for l in ttbar_list:
        canvas = TCanvas(l, l)
        namefile = l+'_'+observable
        rfile = TFile(results_path(year,'groups/MC',namefile+'.root'))
        hist = rfile.Get(namefile)
        hist.Scale(1./nbin)
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
                namefile = l+'_'+observable+'_'+syst+s
                canvas = TCanvas(namefile, namefile)
                rfile = TFile(results_path(year,'groups/SYST',namefile+'.root'))
                hist = rfile.Get(namefile)
                hist.Scale(1./nbin)
                rootfile =  TFile(output_name+'.root', 'UPDATE')
                hist.SetName(l+'_'+syst+s)
                hist.Draw()
                hist.Write()
                rootfile.Close()
