#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import os
import argparse
import numpy as np
import array

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame
from ROOT import TGraphAsymmErrors
from ROOT import TStyle, gStyle, TColor, TLatex, TLine

import tools.tdrstyleNew as tdr
tdr.setTDRStyle()

#python bin/sme_time_reco.py n_bjets cLXX cLXY cLXZ cLYZ signal

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
#parser.add_argument('year', help='year of samples')
parser.add_argument('wilson1', help='wilson coefficent (cL, cR, ...)', default='cL')
parser.add_argument('wilson2', help='wilson coefficent (cL, cR, ...)', default='cL')
parser.add_argument('wilson3', help='wilson coefficent (cL, cR, ...)', default='cL')
parser.add_argument('wilson4', help='wilson coefficent (cL, cR, ...)', default='cL')
parser.add_argument('singletop', help='year of samples', default='')

args = parser.parse_args()
observable = args.observable
#year = args.year
wilson1 = args.wilson1
wilson2 = args.wilson2
wilson3 = args.wilson3
wilson4 = args.wilson4
singletop = args.singletop

wilsonList = []
wilsonList.append(wilson1)
wilsonList.append(wilson2)
wilsonList.append(wilson3)
wilsonList.append(wilson4)

cmunu = 0.01

def textWilson(wilson):
    if (wilson=="cLXX"): modwilson = "c_{L,XX}=#minusc_{L,YY}=" + str(cmunu)
    if (wilson=="cLXY"): modwilson = "c_{L,XY}=c_{L,YX}=" + str(cmunu)
    if (wilson=="cLXZ"): modwilson = "c_{L,XZ}=c_{L,ZX}=" + str(cmunu)
    if (wilson=="cLYZ"): modwilson = "c_{L,YZ}=c_{L,ZY}=" + str(cmunu)

    if (wilson=="cRXX"): modwilson = "c_{R,XX}=#minusc_{R,YY}=" + str(cmunu)
    if (wilson=="cRXY"): modwilson = "c_{R,XY}=c_{R,YX}=" + str(cmunu)
    if (wilson=="cRXZ"): modwilson = "c_{R,XZ}=c_{R,ZX}=" + str(cmunu)
    if (wilson=="cRYZ"): modwilson = "c_{R,YZ}=c_{R,ZY}=" + str(cmunu)

    if (wilson=="cXX"): modwilson = "c_{XX}=#minusc_{YY}=" + str(cmunu)
    if (wilson=="cXY"): modwilson = "c_{XY}=c_{YX}=" + str(cmunu)
    if (wilson=="cXZ"): modwilson = "c_{XZ}=c_{ZX}=" + str(cmunu)
    if (wilson=="cYZ"): modwilson = "c_{YZ}=c_{ZY}=" + str(cmunu)

    if (wilson=="dXX"): modwilson = "d_{XX}=#minusd_{YY}=" + str(cmunu)
    if (wilson=="dXY"): modwilson = "d_{XY}=d_{YX}=" + str(cmunu)
    if (wilson=="dXZ"): modwilson = "d_{XZ}=d_{ZX}=" + str(cmunu)
    if (wilson=="dYZ"): modwilson = "d_{YZ}=d_{ZY}=" + str(cmunu)

    return modwilson

sme_file_2016 = TFile('combine/2016/sme/inputs/n_bjets_sme_all.root')
sme_file_2017 = TFile('combine/2017/sme/inputs/n_bjets_sme_all.root')

def mergeHisto(hist_list):
    nbin = hist_list[0].GetNbinsX()
    print hist_list[0].GetName()
    hist_allbins = TH1F(hist_list[0].GetName()+"_allbins", hist_list[0].GetName()+"_allbins", nbin*2, 0, nbin*2)
    for i in range(len(hist_list)):
        for j in range(nbin):
            hist_allbins.SetBinContent(j+1+i*nbin, hist_list[i].GetBinContent(j+1))
            hist_allbins.SetBinError(j+1+i*nbin, hist_list[i].GetBinError(j+1))
    return hist_allbins



histSM_2016 = sme_file_2016.Get('signal')
histSM_2017 = sme_file_2017.Get('signal')
histSM = mergeHisto([histSM_2016, histSM_2017])

histSME_2016 = []
histSME_2017 = []
histSME = []
histSMEratio = []
for i in range(4):
    histSME_2016.append(sme_file_2016.Get(wilsonList[i]))
    histSME_2017.append(sme_file_2017.Get(wilsonList[i]))
    histSME.append(mergeHisto([histSME_2016[-1], histSME_2017[-1]]))
    histSMEratio.append(histSME[-1].Clone())
    histSMEratio[-1].Divide(histSM)

#histSMEratio[0].Draw()
#raw_input('exit')

histSME_2016_up = []
histSME_2017_up = []
histSME_2016_down = []
histSME_2017_down = []
histSME_up = []
histSME_down = []
histSMEratio_up = []
histSMEratio_down = []


def mergeHisto(hist_list):
    nbin = 24 * 4
    print hist_list[0].GetName()
    hist_allbins = TH1F(hist_list[0].GetName()+"_allbins", hist_list[0].GetName()+"_allbins", nbin*len(hist_list), 0, nbin*len(hist_list))
    for i in range(len(hist_list)):
        for j in range(nbin):
            hist_allbins.SetBinContent(j+1+i*nbin, hist_list[i].GetBinContent(j+1))
            hist_allbins.SetBinError(j+1+i*nbin, hist_list[i].GetBinError(j+1))
    return hist_allbins

for ic in range(4):
    histSME_2016_up.append(TH1F('wilson'+str(ic)+'_2016_up','wilson'+str(ic)+'_2016_up',4*24,0,4*24))
    histSME_2016_down.append(TH1F('wilson'+str(ic)+'_2016_down','wilson'+str(ic)+'_2016_down',4*24,0,4*24))
    histSME_2017_up.append(TH1F('wilson'+str(ic)+'_2017_up','wilson'+str(ic)+'_2017_up',4*24,0,4*24))
    histSME_2017_down.append(TH1F('wilson'+str(ic)+'_2017_down','wilson'+str(ic)+'_2017_down',4*24,0,4*24))

    for ib in range(4):
        histSME_2016_up_ib = sme_file_2016.Get(wilsonList[ic]+'_MCstat_binobs'+str(ib)+'_sme_2016Up')
        histSME_2016_down_ib = sme_file_2016.Get(wilsonList[ic]+'_MCstat_binobs'+str(ib)+'_sme_2016Down')
        histSME_2017_up_ib = sme_file_2017.Get(wilsonList[ic]+'_MCstat_binobs'+str(ib)+'_sme_2017Up')
        histSME_2017_down_ib = sme_file_2017.Get(wilsonList[ic]+'_MCstat_binobs'+str(ib)+'_sme_2017Down')

	for it in range(24):
	    histSME_2016_up[-1].SetBinContent(1+it*4+ib, histSME_2016_up_ib.GetBinContent(1+it*4+ib))
            histSME_2016_down[-1].SetBinContent(1+it*4+ib, histSME_2016_down_ib.GetBinContent(1+it*4+ib))
            histSME_2017_up[-1].SetBinContent(1+it*4+ib, histSME_2017_up_ib.GetBinContent(1+it*4+ib))
            histSME_2017_down[-1].SetBinContent(1+it*4+ib, histSME_2017_down_ib.GetBinContent(1+it*4+ib))

    histSME_up.append(mergeHisto([histSME_2016_up[-1],histSME_2017_up[-1]])) 
    histSME_down.append(mergeHisto([histSME_2016_down[-1],histSME_2017_down[-1]]))
    histSMEratio_up.append(histSME_up[-1].Clone())
    histSMEratio_up[-1].Divide(histSM)
    histSMEratio_down.append(histSME_down[-1].Clone())
    histSMEratio_down[-1].Divide(histSM)

#histSMEratio_up[0].Draw()
#histSMEratio_down[0].Draw()
#histSME_2016_up[0].Divide(histSME_2016[0])
#histSME_2016_up[0].Draw()
#raw_input('exit')

if (observable=="n_bjets"):
    nbin = 4
    min_bin = 0 #1
    max_bin = 4 #5

width_bin = (max_bin-min_bin)/nbin

def getUncertaintyBandGraph(hist,histUp,histDown):
    x = []
    ex_left = []
    ex_right =  []
    y =  []
    ey_up  = []
    ey_down = []
    for i in range(hist.GetNbinsX()):
    #for i in range(nbin):
        x.append(min_bin+width_bin/2.+width_bin*i)
        ex_left.append(width_bin/2.)
        ex_right.append(width_bin/2.)
        y.append(hist.GetBinContent(i+1))
        ey_up.append(abs(histUp.GetBinContent(i+1)-hist.GetBinContent(i+1)))
        ey_down.append(abs(histDown.GetBinContent(i+1)-hist.GetBinContent(i+1)))
    graph_new = TGraphAsymmErrors(len(x),array.array('d', x),array.array('d', y),array.array('d', ex_left),array.array('d', ex_right),array.array('d', ey_down),array.array('d', ey_up))
    graph_new.SetName("uncertainty_band_"+hist.GetName())
    graph_new.SetTitle("uncertainty_band_"+hist.GetName())
    return graph_new

color = []
color.append(2)
color.append(4)
color.append(8)
color.append(797)

UncertaintyBand = []
for ic in range(4):
    print 'UncertaintyBand '+wilsonList[ic]
    UncertaintyBand.append(getUncertaintyBandGraph(histSMEratio[ic], histSMEratio_up[ic], histSMEratio_down[ic]))
    #UncertaintyBand[-1].SetFillColor(color(ic))
    #UncertaintyBand[-1].GetXaxis().SetMinimum(0)
    #UncertaintyBand[-1].GetXaxis().SetMaximum(2*24*4)
    UncertaintyBand[-1].GetXaxis().SetLimits(0,2*24*4)

canvas = TCanvas('SME_comparison','SME_comparison',2400,800)
canvas.UseCurrentStyle()

pad1 = TPad("pad1", "pad1", 0, 0, 1, 1)
#pad1.SetBottomMargin(epsilon)
pad1.SetLeftMargin(0.2)
pad1.SetRightMargin(0.01)

tm = gStyle.GetPadTopMargin()
print 'TopMargin: '+str(tm)+' -> '+str(1.5*tm)
gStyle.SetPadTopMargin(1.5*tm)
pad1.SetTopMargin(1.5*tm)

canvas.cd()
pad1.Draw()
pad1.cd()

for ic in range(4):
    style_histo(UncertaintyBand[ic], color[ic], 1, color[ic], 3001, 0)
    #UncertaintyBand[ic].SetFillColor(color[ic])
    histSMEratio[ic].SetLineColor(color[ic])


#style_histo(UncertaintyBand[0], 2, 1, 2, 3001, 0)

#style_histo(UncertaintyBand[0], 1, 1, 1, 3002, 0)
#histSMEratio[0].Draw("HIST")
#UncertaintyBand[0].Draw("2AP same")
#histSMEratio[0].Draw("HIST same")

Ymin = histSMEratio_up[0].GetMinimum()*(1-0.002)
Ymax = histSMEratio_up[0].GetMaximum()*1.002

print 'Ymin='+str(Ymin)+' Ymax='+str(Ymax)

#UncertaintyBand[0].GetYaxis().SetLimits(Ymin, Ymax)
UncertaintyBand[0].SetMinimum(Ymin)
UncertaintyBand[0].SetMaximum(Ymax)
UncertaintyBand[0].GetYaxis().SetLabelSize(0.03)
UncertaintyBand[0].GetYaxis().SetTitle('N_{t#bar{t},SME}/N_{t#bar{t},SM}')
UncertaintyBand[0].GetYaxis().CenterTitle()
UncertaintyBand[0].GetYaxis().SetTitleOffset(0.8)
UncertaintyBand[0].GetYaxis().SetTitleSize(0.05)

UncertaintyBand[0].GetXaxis().SetTitle("Sidereal hour + 0.25 #times (number of b jets -1)")
#UncertaintyBand[0].GetXaxis().SetTitle("Number of b jets #times sidereal time bin (h)")
UncertaintyBand[0].GetXaxis().SetMaxDigits(0)
UncertaintyBand[0].GetXaxis().CenterTitle()
UncertaintyBand[0].GetXaxis().SetLabelSize(0)#0.15)
UncertaintyBand[0].GetXaxis().SetTitleSize(0.05)
UncertaintyBand[0].GetXaxis().SetLabelOffset(0.01)

UncertaintyBand[0].GetXaxis().SetNdivisions(48,False)
UncertaintyBand[0].GetXaxis().SetTickLength(0)#0.05)

legend_args = (0.005, 0.74, 0.079, 0.92, 'SME model', 'NDC')
#legend_args = (0.005, 0.77, 0.079, 0.92, 'SME model', 'NDC')
legend = TLegend(*legend_args)
#legend.SetTextSize(0.025)
legend.SetTextSize(0.03)
legend.SetBorderSize(0)
legend.AddEntry(UncertaintyBand[0], textWilson(wilsonList[0]), "fl")
legend.AddEntry(UncertaintyBand[1], textWilson(wilsonList[1]), "fl")
legend.AddEntry(UncertaintyBand[2], textWilson(wilsonList[2]), "fl")
legend.AddEntry(UncertaintyBand[3], textWilson(wilsonList[3]), "fl")

UncertaintyBand[0].Draw("2AP")
histSMEratio[0].Draw("HIST same")
UncertaintyBand[1].Draw("2 same")
histSMEratio[1].Draw("HIST same")
UncertaintyBand[2].Draw("2 same")
histSMEratio[2].Draw("HIST same")
UncertaintyBand[3].Draw("2 same")
histSMEratio[3].Draw("HIST same")

legend.Draw("SAME")

line_year = TLine(24*nbin,histSMEratio_up[0].GetMinimum()*(1-0.0015),24*nbin,histSMEratio_up[0].GetMaximum()*1.0015)
line_year.SetLineStyle(9)
line_year.SetLineWidth(2)
line_year.SetLineColor(15)
line_year.Draw("SAME")
latex = TLatex()
latex.SetTextSize(1.2*gStyle.GetPadTopMargin())
latex.SetTextSize(0.03)
latex.DrawLatex(11*nbin,histSMEratio_up[0].GetMaximum()*1.0008,"2016")
latex.DrawLatex(35*nbin,histSMEratio_up[0].GetMaximum()*1.0008,"2017")


line_axis = []
text_axis = []
for i in range(0,48):
    inew = i
    if i>=24:
	inew=i-24
    #else:
    #    inew=i
    print str(i)
    line_axis.append(TLine(i*nbin,Ymin-(Ymax-Ymin)/50.,i*nbin,Ymin+(Ymax-Ymin)/50.))
    #line_axis[-1].SetNDC()
    line_axis[-1].Draw("SAME")
    text_axis.append(TLatex())
    text_axis[-1].SetTextSize(0.025)
    if inew<10:
        text_axis[-1].DrawLatex(i*nbin-0.7,Ymin-(Ymax-Ymin)/50.*2.5,str(inew))
    if inew>=10:
        text_axis[-1].DrawLatex(i*nbin-1.5,Ymin-(Ymax-Ymin)/50.*2.5,str(inew))

#tdr.cmsPrel(77400., 13.,simOnly=True,thisIsPrelim=True)
tdr.cmsPrel(0, 13.,simOnly=True,thisIsPrelim=False)

canvas.Update()

resultname = './results/Comb/other/'+observable+'_sme_time_'+wilsonList[0]+'_'+wilsonList[1]+'_'+wilsonList[2]+'_'+wilsonList[3]+'_fullreco_Comb_'+singletop
print resultname
canvas.SaveAs(resultname+'.pdf')

raw_input('exit')

#def displaySMEcomparison():

