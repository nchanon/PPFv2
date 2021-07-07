#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import argparse
import numpy as np


from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame
from ROOT import TGraphAsymmErrors

import tools.tdrstyle as tdr
tdr.setTDRStyle()

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
#parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
#observable = args.observable
year = args.year

def integral_complete(histo, max_bin):
    return histo.Integral()+histo.GetBinContent(int(max_bin+1))+histo.GetBinContent(0)


nbin = 0
min_bin = 0
max_bin = 0


legend_coordinates = [0.65, 0.75, 0.87, 0.87] 
TH1.SetDefaultSumw2(1)
signal_integral = 0
background_integral_i = []
background_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('differential measurment','differential measurment', 800, 700)
canvas.UseCurrentStyle()


r_2016 = [
    [1.000,0.270,0.367],
    [1.000,0.203,0.292],
    [1.000,0.207,0.302],
    [1.000,0.197,0.277],
    [1.000,0.199,0.282],
    [1.000,0.210,0.312],
    [1.000,0.218,0.334],
    [1.000,0.187,0.252],
    [1.000,0.193,0.268],
    [1.000,0.203,0.292],
    [1.000,0.198,0.279],
    [1.000,0.217,0.332],
    [1.000,0.212,0.316],
    [1.000,0.214,0.322],
    [1.000,0.208,0.305],
    [1.000,0.216,0.329],
    [1.000,0.207,0.303],
    [1.000,0.205,0.298],
    [1.000,0.390,1.010],
    [1.000,0.210,0.311],
    [1.000,0.206,0.300], 
    [1.000,0.210,0.312],
    [1.000,0.222,0.345],
    [1.000,0.208,0.305]
]

r_2017 = [
    [1.000,0.207,0.303],
    [1.000,0.205,0.298],
    [1.000,0.390,1.010],
    [1.000,0.210,0.311],
    [1.000,0.206,0.300],
    [1.000,0.210,0.312],
    [1.000,0.222,0.345],
    [1.000,0.208,0.305],
    [1.000,0.208,0.305],
    [1.000,0.216,0.329],
    [1.000,0.199,0.282],
    [1.000,0.210,0.312],
    [1.000,0.207,0.302],
    [1.000,0.197,0.277],
    [1.000,0.193,0.268],
    [1.000,0.203,0.292],
    [1.000,0.218,0.334],
    [1.000,0.187,0.252],
    [1.000,0.270,0.367],
    [1.000,0.203,0.292],
    [1.000,0.214,0.322],
    [1.000,0.212,0.316],
    [1.000,0.217,0.332],
    [1.000,0.198,0.279]
]

rate = []
if year=='2016':
    rate = r_2016
elif year=='2017':
    rate = r_2017
else:
    print 'ERROR year'
    exit()


def linearized(input, option):
    opt = 0 
    if option=='up' : 
        opt = 2
    elif option=='down' : 
        opt = 1
    else:
        print 'error option'
    foo = []
    for l in input:
        foo.append(l[opt])
    return foo

def squared(list_value):
    somme = 0
    for l in list_value:
        somme += l*l
    return np.sqrt(somme/float(len(list_value)))

################################################################################
## Create Histo 
################################################################################

x = np.array([1 for i in range(24)], dtype='double')
y = np.array([i+1 for i in range(24)], dtype='double')

error_left = np.array([0.01 for i in range(24)], dtype='double')
error_right = np.array([0.01 for i in range(24)], dtype='double')

error_up = np.array(linearized(rate, 'up'), dtype='double')
error_down = np.array(linearized(rate, 'down'), dtype='double')

print 'error + : '+str(squared(error_up))
print 'error - : '+str(squared(error_down))


s = ''
for i in range(1,13):
    s += '&'+str(i)
print s
s = ''
for i in range(13,25):
    s += '&'+str(i)
print s
s = ''
count =0
for l in error_up:
    s += '&'+str(l)
    count += 1
    if(count==12):
        s += '\n' 

print s
print '========'
s = ''
count =0
for l in error_down:
    s += '&'+str(l)
    count += 1
    if(count==12):
        s += '\n' 
print s

hist  = TGraphAsymmErrors(23, y, x , 
                          error_left, error_right, 
                          error_down, error_up) 


################################################################################
## Legend stuff
################################################################################



legend = TLegend(0.5,0.9,0.9,0.75)
legend.SetHeader('Asimov test '+year, 'C')
legend.AddEntry(hist, 't#bar{t} signal strength', 'lep')


#legend_args = (0.645, 0.79, 0.985, 0.85, '', 'NDC')
#legend = TLegend(*legend_args)
#legend.AddEntry(hist, "Data", "l")
#legend.AddEntry(hist_mc, "MC", "l")
#legend_box(legend, legend_coordinates)

################################################################################
## Draw stuff
################################################################################

hist.Draw("ap")
legend.Draw("SAME")

################################################################################
## Set Style
################################################################################

is_center=True

hist.GetYaxis().SetTitle('signal strength #it{#mu}')
hist.GetYaxis().SetRange(0,2)
hist.GetYaxis().SetTitleSize(0.04)
hist.GetYaxis().SetLabelSize(0.04)

hist.GetXaxis().SetTitle('sidereal time ')
hist.GetXaxis().SetRangeUser(0,24)
hist.GetXaxis().SetTitleSize(0.04)
hist.GetXaxis().SetLabelSize(0.04)

if(is_center):
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()

# line_color, line_width, fill_style, marker_style
#style_histo(hist_data, 2, 5, 1, 3001, 1, 1)
#style_histo(hist_mc,   4, 5, 1, 3001, 1, 1)
style_histo(hist,   4, 1, 4, 3005, 1,20)
hist.SetMarkerColor(223)

if(year=='2016'):
    tdr.cmsPrel(35900., 13.)
elif(year=='2017'):
    tdr.cmsPrel(41530., 13.)

################################################################################
## Save
################################################################################

resultname = './results/'+year+'/other/differential_time_'+year

rootfile_output = TFile(resultname+'.root', "RECREATE")
canvas.Write()
canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
rootfile_output.Close()

raw_input('exit')