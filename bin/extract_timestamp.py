#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.sample_manager import *
import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack 
from ROOT import TLegend, TApplication, TRatioPlot, TPad


################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
year = args.year

timestamp = []
samples   = []
n_events  = []

################################################################################
## Functions stuff
################################################################################

def trigger_passed(tree, trigger_list):
    trig = 0
    for l in trigger_list:
        if(tree.GetLeaf(l).GetValue() == 1):
            trig += 1
    if(trig>1):
        return True
    return False

def read_DATA_variable(year, tree, variable, file_name):
    number = 0
    if(trigger_passed(tree, triggers[year])):
        if (file_name.find('MuonEG')!=-1):
            number = tree.GetLeaf(variable).GetValue()
            return number
            if (file_name.find('SingleMuon')!=-1):
                number = tree.GetLeaf(variable).GetValue()
                return number
                if (file_name.find('SingleElectron')!=-1):
                    number = tree.GetLeaf(variable).GetValue()
                    return number
    return number

def complete(word, size):
    foo = word
    remains = size - len(word)
    for i in range(remains):
        foo += ' '
    return foo


################################################################################
## Code body
################################################################################

for sample in sample_DATA[year]:
    rfile = TFile('./inputs/'+year+'/DATA/'+sample+'/tree.root')
    tree  = rfile.Get('events')
    for i in range(tree.GetEntriesFast()):
        tree.GetEntry(i)
        timestamp.append(read_DATA_variable(year, tree, 'unix_time', sample))
        samples.append(sample)
        n_events.append(i)
        if(i%100000 == 0):
            print '100 000 events passed in '+sample

with open('./results/'+year+'/timestamp.txt', 'w') as f:
    for i in range(len(timestamp)):
        print >> f, complete(str(timestamp[i]), 15), complete(str(samples[i]), 30), n_events[i] 