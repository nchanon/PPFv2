#!/usr/bin/env python

import sys
sys.path.append('./')
import math

from tools.style_manager import *
from tools.sample_manager import *
import array

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString, TColor, TGraphAsymmErrors
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame

import tools.tdrstyle as tdr
tdr.setTDRStyle()

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('title', help='display your observable title') # number of b jets
parser.add_argument('timebin', help='display the time bin')

args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title
timebin = int(args.timebin)

def integral_complete(histo, max_bin): 
     return histo.Integral()
#    return histo.Integral()+histo.GetBinContent(int(max_bin+1))+histo.GetBinContent(0)


nbin = 0
min_bin = 0
max_bin = 0
#legend_coordinates = observable_values(observable)[1]
TH1.SetDefaultSumw2(1)
signal_integral = 0
background_integral_i = []
background_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('stack_'+observable,'stack_'+observable, 800, 800)
canvas.UseCurrentStyle()

datafile_input = TFile('./results/'+year+'/flattree/'+observable+'_data_forComp.root')
stimebin="";
if (timebin==-1):
     stimebin = "_puold";
if (timebin==-2):
     stimebin = "_punew";
if (timebin==-3):
     stimebin = "_puinc";
if (timebin>=0):
     stimebin = "_put"+str(timebin);
rootfile_input = TFile('./results/'+year+'/flattree/'+observable+'_forComp'+stimebin+'.root')

################################################################################
## Create Histo 
################################################################################

doLog = False

r = 0.3
epsilon = 0.1

pad1 = TPad("pad1", "pad1", 0, r-epsilon, 1, 1)
pad1.SetBottomMargin(epsilon)
canvas.cd()
if (doLog): pad1.SetLogy()
pad1.Draw()
pad1.cd()


###########
# data part
###########
hist_data = datafile_input.Get('data_obs')
data_integral = integral_complete(hist_data, max_bin)


# convenient variables
nbin    = hist_data.GetNbinsX()
min_bin = hist_data.GetXaxis().GetXmin()
max_bin = hist_data.GetXaxis().GetXmax()
width_bin = (max_bin-min_bin)/nbin

###########
# mc part
###########

# signal
hist_signal = rootfile_input.Get('signal')
signal_integral = integral_complete(hist_signal, max_bin)

# backgrounds
hist_singletop = rootfile_input.Get('singletop')
hist_ttx = rootfile_input.Get('ttx')
hist_dibosons = rootfile_input.Get('dibosons')
hist_vjets = rootfile_input.Get('vjets')

hist_nom = []
hist_nom.append(hist_signal)
hist_nom.append(hist_singletop)
hist_nom.append(hist_ttx)
hist_nom.append(hist_dibosons)
hist_nom.append(hist_vjets)

background_integral = 0
hist_background  = TH1F("background", "background", nbin, min_bin, max_bin)
for l in rootfile_input.GetListOfKeys():
    for s in ttbar_list:
        if(l.GetName() == s and not l.GetName() == 'signal'):
            #print l.GetName()
            hist_background.Add(rootfile_input.Get(l.GetName()))
            background_integral_i.append([integral_complete(rootfile_input.Get(l.GetName()), max_bin), l.GetName()])
            background_integral += integral_complete(rootfile_input.Get(l.GetName()), max_bin)

##################
# Uncertainty band
##################

# The uncertainty band is wrong: please use pre-fit plot in combine-ttbar instead
'''
rootfile_forsyst = TFile('./combine/'+year+'/one_bin/inputs/'+observable+'_24_0.root')

hist_ref  = []
hist_systUp = []
hist_systDown = []

for l in rootfile_forsyst.GetListOfKeys():
    for proc in ttbar_list:
        if(l.GetName() == proc): #main histo
	     print('Ref '+l.GetName())
	     #h = TH1F(l.GetName(),l.GetName(),nbin, min_bin, max_bin)
	     #h.Add(rootfile_forsyst.Get(l.GetName()))
	     hist_ref.append(rootfile_forsyst.Get(l.GetName()))
	elif (l.GetName() != proc and l.GetName().find(proc)!= -1): #syst_histo
	     print('Syst '+l.GetName())
	     if (l.GetName().find("Up")!= -1): hist_systUp.append(rootfile_forsyst.Get(l.GetName()))
	     if (l.GetName().find("Down")!= -1): hist_systDown.append(rootfile_forsyst.Get(l.GetName()))

#print hist_systUp
#print hist_systDown

hist_totUp = []
hist_totDown = []

for hist_proc in hist_ref:
    print(hist_proc.GetName()+" : " +str(hist_proc.GetBinContent(1))+" "+str(hist_proc.GetBinContent(2))+" "+str(hist_proc.GetBinContent(3))+" "+str(hist_proc.GetBinContent(4))+" "+str(hist_proc.GetBinContent(5))+" "+str(hist_proc.GetBinContent(6))+" "+str(hist_proc.GetBinContent(7)))
    name_totUp = hist_proc.GetName() + '_TotUp'
    hTotUp = TH1F(name_totUp,name_totUp,nbin, min_bin, max_bin)
    name_totDown = hist_proc.GetName() + '_TotDown'
    hTotDown = TH1F(name_totDown,name_totDown,nbin, min_bin, max_bin)
    for hist_up in hist_systUp:
        if (hist_up.GetName().find(hist_proc.GetName())!=-1):
	    #print hist_up.GetName()
	    for i in range(nbin):
		diff = math.sqrt((hist_up.GetBinContent(1+i)-hist_proc.GetBinContent(1+i))*(hist_up.GetBinContent(1+i)-hist_proc.GetBinContent(1+i)))
		val = math.sqrt(hTotUp.GetBinContent(1+i)*hTotUp.GetBinContent(1+i)+diff*diff)
		#print("diff="+str(diff)+" val="+str(val))
		hTotUp.SetBinContent(1+i, val)
    for i in range(nbin): 
    	hTotUp.SetBinContent(i+1, hTotUp.GetBinContent(i+1)/hist_proc.GetBinContent(i+1))
    hist_totUp.append(hTotUp)
    print(hist_proc.GetName()+"_totUp : " +str(hist_totUp[-1].GetBinContent(1))+" "+str(hist_totUp[-1].GetBinContent(2))+" "+str(hist_totUp[-1].GetBinContent(3))+" "+str(hist_totUp[-1].GetBinContent(4))+" "+str(hist_totUp[-1].GetBinContent(5))+" "+str(hist_totUp[-1].GetBinContent(6))+" "+str(hist_totUp[-1].GetBinContent(7)))
    for hist_down in hist_systDown:
        if (hist_down.GetName().find(hist_proc.GetName())!=-1):
            #print hist_down.GetName()
            for i in range(nbin):
                diff = math.sqrt((hist_down.GetBinContent(1+i)-hist_proc.GetBinContent(1+i))*(hist_down.GetBinContent(1+i)-hist_proc.GetBinContent(1+i)))
                val = math.sqrt(hTotDown.GetBinContent(1+i)*hTotDown.GetBinContent(1+i)+diff*diff)
                #print("diff="+str(diff)+" val="+str(val))
                hTotDown.SetBinContent(1+i, val)
    for i in range(nbin): 
        hTotDown.SetBinContent(i+1, hTotDown.GetBinContent(i+1)/hist_proc.GetBinContent(i+1))
    hist_totDown.append(hTotDown)
    print(hist_proc.GetName()+"_totDown : " +str(hist_totDown[-1].GetBinContent(1))+" "+str(hist_totDown[-1].GetBinContent(2))+" "+str(hist_totDown[-1].GetBinContent(3))+" "+str(hist_totDown[-1].GetBinContent(4))+" "+str(hist_totDown[-1].GetBinContent(5))+" "+str(hist_totDown[-1].GetBinContent(6))+" "+str(hist_totDown[-1].GetBinContent(7)))


hist_signal_up = hist_signal.Clone()
hist_signal_up.SetName("hist_signal_up")
hist_singletop_up = hist_singletop.Clone()
hist_singletop_up.SetName("hist_singletop_up")
hist_ttx_up = hist_ttx.Clone()
hist_ttx_up.SetName("hist_ttx_up")
hist_dibosons_up = hist_dibosons.Clone()
hist_dibosons_up.SetName("hist_dibosons_up")
hist_vjets_up = hist_vjets.Clone()
hist_vjets_up.SetName("hist_vjets_up")

hist_signal_down = hist_signal.Clone()
hist_signal_down.SetName("hist_signal_down")
hist_singletop_down = hist_singletop.Clone()
hist_singletop_down.SetName("hist_singletop_down")
hist_ttx_down = hist_ttx.Clone()
hist_ttx_down.SetName("hist_ttx_down")
hist_dibosons_down = hist_dibosons.Clone()
hist_dibosons_down.SetName("hist_dibosons_down")
hist_vjets_down = hist_vjets.Clone()
hist_vjets_down.SetName("hist_vjets_down")


for hist_proc in hist_nom:
    print(hist_proc.GetName()+" : " +str(hist_proc.GetBinContent(1))+" "+str(hist_proc.GetBinContent(2))+" "+str(hist_proc.GetBinContent(3))+" "+str(hist_proc.GetBinContent(4))+" "+str(hist_proc.GetBinContent(5))+" "+str(hist_proc.GetBinContent(6))+" "+str(hist_proc.GetBinContent(7)))

sum_Down = []
sum_Up = []
sum_MCstat = []
h_tmp = TH1F("tmp","tmp",nbin,min_bin,max_bin)
for i in range(nbin):
    sum_i_Up  = 0
    sum_i_Down = 0
    val_i_Nom = 0
    val_i_MCstat = 0
    sum_i_MCstat = 0
    for hist_Down in hist_totDown:
	if (hist_Down.GetName().find("signal")!=-1): h_tmp = hist_signal
	if (hist_Down.GetName().find("hist_singletop")!=-1): h_tmp = hist_singletop
	if (hist_Down.GetName().find("ttx")!=-1): h_tmp = hist_ttx
	if (hist_Down.GetName().find("dibosons")!=-1): h_tmp = hist_dibosons
	if (hist_Down.GetName().find("vjets")!=-1): h_tmp = hist_vjets
	val_i_Nom = h_tmp.GetBinContent(i+1)
	val_i_MCstat = h_tmp.GetBinError(i+1)
	sum_i_MCstat = math.sqrt(sum_i_MCstat*sum_i_MCstat+val_i_MCstat*val_i_MCstat)
        sum_i_Down = math.sqrt(sum_i_Down*sum_i_Down+val_i_Nom*val_i_Nom*hist_Down.GetBinContent(i+1)*hist_Down.GetBinContent(i+1)+val_i_MCstat*val_i_MCstat)
    sum_Down.append(sum_i_Down)
for i in range(nbin):
    sum_i_Up  = 0
    sum_i_Down = 0
    val_i_Nom = 0
    val_i_MCstat = 0
    sum_i_MCstat = 0
    for hist_Up in hist_totUp:
        if (hist_Up.GetName().find("signal")!=-1): h_tmp = hist_signal
        if (hist_Up.GetName().find("hist_singletop")!=-1): h_tmp = hist_singletop
        if (hist_Up.GetName().find("ttx")!=-1): h_tmp = hist_ttx
        if (hist_Up.GetName().find("dibosons")!=-1): h_tmp = hist_dibosons
        if (hist_Up.GetName().find("vjets")!=-1): h_tmp = hist_vjets
        val_i_Nom = h_tmp.GetBinContent(i+1)
        val_i_MCstat = h_tmp.GetBinError(i+1)
        sum_i_MCstat = math.sqrt(sum_i_MCstat*sum_i_MCstat+val_i_MCstat*val_i_MCstat)
        sum_i_Up = math.sqrt(sum_i_Up*sum_i_Up+val_i_Nom*val_i_Nom*hist_Up.GetBinContent(i+1)*hist_Up.GetBinContent(i+1)+val_i_MCstat*val_i_MCstat)
    sum_Up.append(sum_i_Up)
    sum_MCstat.append(sum_i_MCstat) #Need to do it only for up, not for down (since up and down MCstat uncertainty are the same)

x = []
ex_left = []
ex_right =  []
y = []
ey_up = []
ey_down = []
yratio = []
ey_ratio_down = []
ey_ratio_up = []

for i in range(nbin):
    x.append(min_bin+width_bin/2.+width_bin*i)
    ex_left.append(width_bin/2.)
    ex_right.append(width_bin/2.)
    y.append(hist_background.GetBinContent(i+1)+hist_signal.GetBinContent(i+1))
    print("Bin"+str(i+1)+" : Up="+ str(sum_Up[i])+" Down="+str(sum_Down[i])+" (with MCStat="+str(sum_MCstat[i]/y[i])+" of bin content)")
    ey_up.append(sum_Up[i])
    ey_down.append(sum_Down[i])
    yratio.append(1)
    ey_ratio_up.append(sum_Up[i]/y[i])
    ey_ratio_down.append(sum_Down[i]/y[i])

UncertaintyBand = TGraphAsymmErrors(len(x),array.array('d', x),array.array('d', y),array.array('d', ex_left),array.array('d', ex_right),array.array('d', ey_down),array.array('d', ey_up))

UncertaintyBandRatio = TGraphAsymmErrors(len(x),array.array('d', x),array.array('d', yratio),array.array('d', ex_left),array.array('d', ex_right),array.array('d', ey_ratio_down),array.array('d', ey_ratio_up))
'''
################################################################################
## Legend stuff
################################################################################

legend_args = (0.73, 0.65, 0.93, 0.92, '', 'NDC')
legend = TLegend(*legend_args)
legend.AddEntry(hist_signal, "t#bar{t} SM", "f")
#legend.AddEntry(hist_background, "non-t#bar{t}", "f")
legend.AddEntry(hist_singletop, "single top", "f")
legend.AddEntry(hist_vjets, "W/Z+jets", "f")
legend.AddEntry(hist_dibosons, "Dibosons", "f")
legend.AddEntry(hist_ttx, "t#bar{t}+X", "f")
legend.AddEntry(hist_data, "data")
#legend_box(legend, legend_coordinates)

################################################################################
## Draw stuff
################################################################################

stack = THStack()
#stack.Add(hist_background)
#stack.Add(hist_signal)
#stack.Add(hist_singletop)
#stack.Add(hist_ttx)
#stack.Add(hist_dibosons)
stack.Add(hist_ttx)
stack.Add(hist_dibosons)
stack.Add(hist_vjets)
#stack.Add(hist_dibosons)
#stack.Add(hist_ttx)
stack.Add(hist_singletop)
stack.Add(hist_signal)
if (doLog): stack.SetMinimum(10)

#stack.GetHistogram().GetXaxis().SetRange(min_bin,max_bin)
#UncertaintyBand.GetXaxis().SetRangeUser(min_bin,max_bin)
#UncertaintyBand.SetMinimum(0)
#if (doLog): UncertaintyBand.SetMinimum(10)

stack.Draw()
#UncertaintyBand.Draw("2AP SAME")
stack.Draw("HIST SAME")
hist_data.Draw("E SAME")
legend.Draw("SAME")

################################################################################
## Set Style
################################################################################

# line_color, line_width, fill_color, fill_style, marker_size, marker_style=1
style_histo(hist_signal, 2, 1, 2, 3004, 0)
style_histo(hist_singletop, 4, 1, 4, 3005, 0)
style_histo(hist_ttx, 8, 1, 8, 3005, 0)
style_histo(hist_dibosons, 42, 1, 42, 3005, 0)
style_histo(hist_vjets, 619, 1, 619, 3005, 0)
#style_histo(hist_background, 4, 1, 4, 3005, 0)
style_histo(hist_data, 1, 1, 0, 3001, 1, 20)

#style_histo(UncertaintyBand, 1, 1, 1, 3005, 0)
#style_labels_counting(UncertaintyBand, 'Events', title)
#UncertaintyBand.GetXaxis().SetLabelSize(0)
#UncertaintyBand.GetXaxis().SetTitleSize(0)

style_labels_counting(stack, 'Events', title)
stack.GetXaxis().SetLabelSize(0)
stack.GetXaxis().SetTitleSize(0)

if(year=='2016'):
    tdr.cmsPrel(35900., 13.,simOnly=False,thisIsPrelim=True)
elif(year=='2017'):
   tdr.cmsPrel(41500., 13.,simOnly=False,thisIsPrelim=True)

################################################################################
## Ratio
################################################################################

pad2 = TPad("pad2", "pad2", 0, 0, 1, r*(1-epsilon))
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0.4)
pad2.SetFillStyle(0)
canvas.cd()
pad2.Draw()
pad2.cd()

ratio_coef = 0.3

h_one = TH1F("one", "one", 1, min_bin, max_bin)
h_one.SetBinContent(1, 1)
h_one.SetLineWidth(1)
h_one.SetLineColor(15)
h_num = hist_data.Clone()
h_denom = hist_signal+hist_background
h_num.Divide(h_denom)
h_num.GetXaxis().SetTitle("aksjd")
ratio = THStack()
ratio.Add(h_num)

ratio.SetMaximum(1+ratio_coef)
ratio.SetMinimum(1-ratio_coef)
ratio.Draw()
h_one.Draw("SAME")

#UncertaintyBandRatio.Draw("2A SAME")
#h_num.Draw("SAME")
#h_one.Draw("SAME")

style_labels_counting(ratio, 'Data/mc', title)
ratio.GetYaxis().SetLabelSize(0.1)
ratio.GetYaxis().SetTitleSize(0.1)
ratio.GetYaxis().SetTitleOffset(0.5)
ratio.GetXaxis().SetLabelSize(0.15)
ratio.GetXaxis().SetTitleSize(0.17)
ratio.GetXaxis().SetLabelOffset(0.01)

#style_histo(UncertaintyBandRatio, 1, 1, 1, 3005, 0)
#UncertaintyBandRatio.GetXaxis().SetRangeUser(min_bin,max_bin)
#UncertaintyBandRatio.SetMinimum(1-ratio_coef)
#UncertaintyBandRatio.SetMaximum(1+ratio_coef)

#style_labels_counting(UncertaintyBandRatio, 'Ratio data/mc', title)
#UncertaintyBandRatio.GetYaxis().SetLabelSize(0.1)
#UncertaintyBandRatio.GetYaxis().SetTitleSize(0.1)
#UncertaintyBandRatio.GetYaxis().SetTitleOffset(0.5)
#UncertaintyBandRatio.GetXaxis().SetLabelSize(0.15)
#UncertaintyBandRatio.GetXaxis().SetTitleSize(0.17)
#UncertaintyBandRatio.GetXaxis().SetLabelOffset(0.01)


################################################################################
## Save
################################################################################

resultname = './results/'+year+'/comparaison/'+observable+'_'+year+'_splitted'

rootfile_output = TFile(resultname+stimebin+'.root', "RECREATE")
canvas.Write()
canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
rootfile_output.Close()


print 'Signal integral     : ','%.2f'%signal_integral
for i in background_integral_i:
    print i[1], '      : ', '%.2f'%i[0]
print 'Total Background integral : ', '%.2f'%background_integral
print 'Total MC integral : ', '%.2f'%(signal_integral+background_integral)
print 'Data integral       : ', data_integral
print 'Data/MC agreement  : ', '%.1f'%(100*(signal_integral+background_integral-data_integral)/(signal_integral+background_integral)), '%'

#raw_input('exit')
