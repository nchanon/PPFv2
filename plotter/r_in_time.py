import os, sys
import argparse 
import numpy as np
sys.path.append('./')

from ROOT import TH1F, TFile, TCanvas, TLegend, THStack

from tools.style_manager import *
import tools.tdrstyle as tdr
tdr.setTDRStyle()

###################
## Initialisation
###################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('asimov',nargs='?', help='set if asimov test', default='')

args = parser.parse_args()
observable = args.observable
year = args.year
asimov = args.asimov

if asimov == 'asimov':
    asimov = '_asimov'
    head = 'Asimov '
else:
    asimov = ''
    head = ''

corr = '_uncorrelated'
file_name = './combine/'+year+'/one_bin/results/histograms_24bins_time_'+year+'.root'
get_name  = 'diff_'+observable+asimov+corr


###################
## Inputs
###################

#observable = raw_input('Choose your observable : ')
#print observable


###################
## Body
###################

canvas = TCanvas(observable, observable, 600, 600)
inputfile = TFile(file_name)
histogram = inputfile.Get(get_name)


if observable.find('m_dilep') != -1:
    title = 'dilepton mass' 
elif  observable.find('n_bjets') != -1:
    title = 'b-jets multiplicity'

leg = 't#bar{t} signal strength'
header = head+year

# cosmetic
histogram.SetName('')
histogram.SetTitle(title)

legend = TLegend(0.5,0.9,0.9,0.75)
legend.SetHeader(header, 'C')
legend.AddEntry(histogram, leg, 'lep')

stack = THStack()
stack.Add(histogram)
stack.Draw('E HIST')
legend.Draw('SAME')

stack.SetMinimum(0.95)
stack.SetMaximum(1.05)

style_histo(histogram, 1, 1, 0, 3001, 1, 20)
if(year == '2017'):
    tdr.cmsPrel(41530., 13.)
elif(year == '2016'):
    tdr.cmsPrel(35900., 13.)
style_labels_counting(stack, 'r (signal strengh)', 'sidereal time (h)')


raw_input('exit')

canvas.SaveAs('./results/'+year+'/stats/'+observable+corr+'_in_time_'+year+'.png')