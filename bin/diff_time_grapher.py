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
    [1.000,0.311,0.165],
    [1.000,0.172,0.200],
    [1.000,0.311,0.194],
    [1.000,0.195,0.230],
    [1.000,0.311,0.194],
    [1.000,0.192,0.227],
    [1.000,0.311,0.188],
    [1.000,0.203,0.242],
    [1.000,0.191,0.226],
    [1.000,0.314,0.174],
    [1.000,0.313,0.218],
    [1.000,0.311,0.193],
    [1.000,0.312,0.170],
    [1.000,0.205,0.244],
    [1.000,0.194,0.229],
    [1.000,0.312,0.176],
    [1.000,0.312,0.199],
    [1.000,0.311,0.194],
    [1.000,0.212,0.253],
    [1.000,0.312,0.213],
    [1.000,0.312,0.210],
    [1.000,0.312,0.210],
    [1.000,0.311,0.190],
    [1.000,0.313,0.219]
]

r_2017 = [
    [1.000,0.263,0.317],
    [1.000,0.302,0.254],
    [1.000,0.294,0.263],
    [1.000,0.301,0.241],
    [1.000,0.301,0.245],
    [1.000,0.294,0.271],
    [1.000,0.295,0.290],
    [1.000,0.299,0.220],
    [1.000,0.292,0.234],
    [1.000,0.302,0.254],
    [1.000,0.301,0.243],
    [1.000,0.304,0.288],
    [1.000,0.295,0.274],
    [1.000,0.295,0.280],
    [1.000,0.302,0.265],
    [1.000,0.303,0.285],
    [1.000,0.303,0.263],
    [1.000,0.302,0.259],
    [1.000,0.405,0.478],
    [1.000,0.294,0.270],
    [1.000,0.302,0.260],
    [1.000,0.294,0.271],
    [1.000,0.304,0.298],
    [1.000,0.294,0.265]

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