#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import os
import argparse
import numpy as np


from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame
from ROOT import TGraphAsymmErrors
from ROOT import TStyle, gStyle, TColor, TLatex

import tools.tdrstyle as tdr
tdr.setTDRStyle()

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
#parser.add_argument('wilson', help='wilson coefficent (cLXX, cRXX, ...)', default='cLXX')
parser.add_argument('process', help='signal or singletop or signal_dilep_LO', default='signal')
parser.add_argument('timebin', help='display the time bin')


args = parser.parse_args()
observable = args.observable
year = args.year
#wilson = args.wilson
process = args.process
timebin = int(args.timebin)

stimebin="";
if (timebin==-1):
     stimebin = "_puold";
if (timebin==-2):
     stimebin = "_punew";
if (timebin==-3):
     stimebin = "_puinc";
if (timebin>=0):
     stimebin = "_put"+str(timebin);

if process=='signal' or process=='singletop':
    f = TFile("results/"+year+"/flattree/"+observable+stimebin+".root")
elif process=='signal_dilep_LO':
    f = TFile("results/"+year+"/flattree/"+observable+"_alt"+stimebin+".root")

def plotResponseMatrix(SMEparam):

    if SMEparam=="":
        h = f.Get(process+"_responseMatrix")
    else:
	h = f.Get(process+"_responseMatrix_"+SMEparam)

    nbinX = h.GetNbinsX()
    nbinY = h.GetNbinsY()

    doNormalizeByColumn = True

    if doNormalizeByColumn:
	for ix in range(nbinX):
	    area = h.Integral(1+ix, 1+ix, 1, nbinY)
	    for iy in range(nbinY):
		bincontent = h.GetBinContent(ix+1,iy+1)
		h.SetBinContent(ix+1, iy+1, bincontent/area)

    canvas = TCanvas("responseMatrix","responseMatrix",600,800)

    pad = TPad("pad","pad",0,0,1,1)
    pad.SetLeftMargin(0.15)
    pad.SetBottomMargin(0.1)
    pad.SetRightMargin(0.14)
    pad.Draw()
    pad.cd()

    h.GetXaxis().SetBinLabel(1,"1")
    h.GetXaxis().SetBinLabel(2,"2")
    h.GetXaxis().SetBinLabel(3,"3")
    h.GetXaxis().SetBinLabel(4,"#geq4")


    h.GetYaxis().SetBinLabel(1,"OOA")
    h.GetYaxis().SetBinLabel(2,"0")
    h.GetYaxis().SetBinLabel(3,"1")
    h.GetYaxis().SetBinLabel(4,"2")
    h.GetYaxis().SetBinLabel(5,"3")
    h.GetYaxis().SetBinLabel(6,"#geq4")

    gStyle.SetPaintTextFormat("4.2f");

    h.GetXaxis().SetLabelSize(0.06)
    h.GetYaxis().SetLabelSize(0.06)
    #h.GetZaxis().SetLabelSize(0.025)
    h.SetMarkerSize(1.5)

    h.GetYaxis().SetTitleOffset(1.5)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    #h.SetTitle("(normalized by column)")
    h.SetYTitle("Number of generator-level b jets")
    h.SetXTitle("Number of reconstructed b jets")
    h.Draw("COLZTEXT")

    latex = TLatex()
    latex.SetTextSize(0.55*gStyle.GetPadTopMargin());

    latex.SetNDC();
    latex.DrawLatex(0.75,0.98,"Fraction of events");
    latex.DrawLatex(0.71,0.96,"(normalized by column)");


    latex.SetTextSize(0.85*gStyle.GetPadTopMargin());
    latex.DrawLatex(0.5,0.96,year);

    if SMEparam=="":
        canvas.SaveAs("results/"+year+"/other/"+observable+"_responseMatrix_"+process+"_"+year+".pdf")
    else:
	canvas.SaveAs("results/"+year+"/other/"+observable+"_responseMatrix_"+process+"_"+year+"_"+SMEparam+".pdf")

plotResponseMatrix("") #SM

if process=='signal_dilep_LO':
    scoeff_list = ["cL","cR","c","d"]
    sdir_list = ["XX","XY","XZ","YZ"]
    for scoeff in scoeff_list:
	for sdir in sdir_list:
	    for it in range(24):
		print scoeff+"_"+sdir+"_"+str(it)    
	        plotResponseMatrix(scoeff+"_"+sdir+"_"+str(it))

raw_input()
exit()
