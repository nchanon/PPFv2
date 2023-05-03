#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
#from tools.sample_manager import *

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
parser.add_argument('wilson', help='wilson coefficent (cLXX, cRXX, ...)', default='cLXX')
parser.add_argument('singletop', help='signal or singletop', default='signal')
parser.add_argument('timebin', help='display the time bin')

args = parser.parse_args()
observable = args.observable
year = args.year
wilson = args.wilson
singletop = args.singletop
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

print observable+' '+year+' '+wilson+' '+singletop+' '+str(timebin)

legend_coordinates = [0.65, 0.75, 0.87, 0.87] 
TH1.SetDefaultSumw2(1)
signal_integral = 0
background_integral_i = []
background_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('sme modulations','sme modulations', 700, 700)
canvas.UseCurrentStyle()
gStyle.SetPalette(55);

#doResponseMatrix = False
doResponseMatrix = True

doResponseMatrixLO = True

doLOparticle = True #from dedicated MiniAOD
#doLOparticle = False

#doLOreco =  True #from dedicated MiniAOD
doLOreco = False


################################################################################
## Create Histo 
################################################################################

prefix=''
suffix=''

if singletop=='singletop':
    prefix='singletop_'
    suffix='_singletop'

if doLOreco:
    prefix = 'reco_'+year+'_'

if doLOparticle:
    prefix = 'signal_miniAOD_LO_'+year+'_'

cmunu = 0.01

sme_file = TFile('./results/'+year+'/flattree/'+observable+'_sme.root')
#hist = sme_file.Get(prefix+wilson+'_details')
#hist.Scale(cmunu)

smassbin = [20,60,100,140,180,220,260]
sptemubin = [0,35,70,105,140,175,210,245]
snbjetsbin = [0, 1, 2, 3, 4, 5]

if observable=='m_dilep':
    svar = 'm_{ll}'
    skinbin = smassbin
elif observable=='pt_emu':
    svar = 'p_{T,ll}'
    skinbin = sptemubin
elif observable=='n_bjets':
    svar = 'n_{bjets}'
    skinbin = snbjetsbin

hist_kinbin = []
for i in range(8):
    if (observable=='m_dilep' and i>0):
        hist_kinbin.append(sme_file.Get(prefix+wilson+'_details'+str(i)))
	hist_kinbin[-1].Scale(cmunu)
    if (observable=='pt_emu' and i<7):
        hist_kinbin.append(sme_file.Get(prefix+wilson+'_details'+str(i)))
        hist_kinbin[-1].Scale(cmunu)
    if (observable=='n_bjets' and i<6 and doLOreco==False):
        hist_kinbin.append(sme_file.Get(prefix+wilson+'_details'+str(i)))
        hist_kinbin[-1].Scale(cmunu)
    if (observable=='n_bjets' and i>=1 and i<=4 and doLOreco==True):
	hist_kinbin.append(sme_file.Get(prefix+wilson+'_details'+str(i)))
	hist_kinbin[-1].Scale(cmunu)


nhist_max = 7
if (observable=='n_bjets'):
    nhist_max = 6

################################################################################
## Apply response matrix 
################################################################################

doNormalizeByColumn = True

hist_kinbin_rec = []

if doResponseMatrix:
    if doResponseMatrixLO==False:
        fResponseMatrix = TFile("results/"+year+"/flattree/"+observable+stimebin+".root")
        hResponseMatrix = fResponseMatrix.Get(singletop+"_responseMatrix")
    if doResponseMatrixLO==True:
	fResponseMatrix = TFile("results/"+year+"/flattree/"+observable+"_alt"+stimebin+".root")
        hResponseMatrix = fResponseMatrix.Get("signal_dilep_LO_responseMatrix")

    nbinX = hResponseMatrix.GetNbinsX()
    nbinY = hResponseMatrix.GetNbinsY()
    #print('ResponseMatrix nbinX='+str(nbinX)+' nbinY='+str(nbinY))

    if doNormalizeByColumn:
        for ix in range(nbinX):
            area = hResponseMatrix.Integral(1+ix, 1+ix, 1, nbinY)
            for iy in range(nbinY):
                bincontent = hResponseMatrix.GetBinContent(ix+1,iy+1)
                hResponseMatrix.SetBinContent(ix+1, iy+1, bincontent/area)

    for ix in range(nbinX):
	hist_kinbin_rec.append(hist_kinbin[ix].Clone())
	for it in range(hist_kinbin[ix].GetNbinsX()):
	    recbincontent = 0
	    for iy in range(nbinY):
		if iy==0:
		    genbin = nbinY-1
		else:
		    genbin = iy-1
		#print('ix='+str(ix)+' iy='+str(iy)+' genbin='+str(genbin))
	        recbincontent += hResponseMatrix.GetBinContent(1+ix,1+iy)*hist_kinbin[genbin].GetBinContent(1+it)	    
	    hist_kinbin_rec[ix].SetBinContent(1+it, recbincontent)
	hist_kinbin[ix] = hist_kinbin_rec[ix]

if (observable=='n_bjets' and (doResponseMatrix or doLOreco)):
    nhist_max = 4

################################################################################
## Legend stuff
################################################################################

slegendtitle = ''
if doResponseMatrix==False and doLOreco==False:
    slegendtitle = 'Generator level'
if doResponseMatrix==True or doLOreco==True:
    slegendtitle = 'Reconstructed level'


legend = TLegend(0.68,0.64,0.9,0.94, slegendtitle)
legend.SetBorderSize(0)
if doResponseMatrix==False:
    legend.SetTextSize(0.023)
else:
    legend.SetTextSize(0.027)


if (observable!='n_bjets'):
    #legend.AddEntry(hist, 'inclusive ('+svar+'>'+str(skinbin[0])+' GeV)', 'l')
    for i in range(nhist_max):
        legend.AddEntry(hist_kinbin[i], str(skinbin[i])+'<'+svar+'<'+str(skinbin[i+1])+' GeV', 'l')
    legend.AddEntry(hist_kinbin[6], svar+'>'+str(skinbin[6])+' GeV', 'l')

if (observable=='n_bjets'):
    if doResponseMatrix==False and doLOreco==False:
	#legend.AddEntry(hist, 'inclusive ('+svar+'>='+str(skinbin[0])+')', 'l')
        for i in range(nhist_max-2):
            legend.AddEntry(hist_kinbin[i], svar+'='+str(skinbin[i]), 'l')
        legend.AddEntry(hist_kinbin[nhist_max-2], svar+'#geq'+str(skinbin[nhist_max-2]), 'l')
	legend.AddEntry(hist_kinbin[nhist_max-1], 'out of acceptance', 'l')

    if doResponseMatrix==True or doLOreco:
        for i in range(nhist_max-1):
            legend.AddEntry(hist_kinbin[i], svar+'='+str(skinbin[i+1]), 'l')
        legend.AddEntry(hist_kinbin[nhist_max-1], svar+'#geq'+str(skinbin[nhist_max]), 'l')



################################################################################
## Draw stuff
################################################################################

min = 99.
max = -99.

#if doResponseMatrix==False:
#    min = hist.GetMinimum()
#    max = hist.GetMaximum()
#    hist.SetLineColor(1)
#    hist.SetLineWidth(2)
#    hist.Draw("hist")

for i in range(nhist_max):
    if (hist_kinbin[i].GetMinimum() < min): min = hist_kinbin[i].GetMinimum()
    if (hist_kinbin[i].GetMaximum() > max): max = hist_kinbin[i].GetMaximum()
    if (observable!='n_bjets'):
	hist_kinbin[i].SetLineColor(gStyle.GetColorPalette((i+1)*255/5))
    else:
	if doResponseMatrix==False:
	    hist_kinbin[i].SetLineColor(gStyle.GetColorPalette(255-(i+1)*255/7))
	if doResponseMatrix==True or doLOreco:
	    hist_kinbin[i].SetLineColor(gStyle.GetColorPalette((i+1)*255/5))
    hist_kinbin[i].SetLineWidth(2)
    #if doResponseMatrix==True:
    if i==0:
	hist_kinbin[i].Draw("hist")
    else:
        hist_kinbin[i].Draw("histsame")

#if doResponseMatrix==False:
#    min = hist.GetMinimum()
#    max = hist.GetMaximum()
#    hist.SetLineColor(1)
#    hist.SetLineWidth(2)
#    hist.Draw("histsame")

legend.Draw("SAME")

################################################################################
## Set Style
################################################################################

is_center=True

hist_kinbin[0].GetYaxis().SetRangeUser(min*1.2,max*1.2)
hist_kinbin[0].GetYaxis().SetTitle("f(t) = #sigma_{SME} / #sigma_{SM} - 1")
hist_kinbin[0].GetYaxis().SetTitleSize(0.04)
hist_kinbin[0].GetYaxis().SetLabelSize(0.04)

hist_kinbin[0].GetXaxis().SetRangeUser(0,24)
hist_kinbin[0].GetXaxis().SetTitle('Sidereal time (h)')
hist_kinbin[0].GetXaxis().SetTitleSize(0.04)
hist_kinbin[0].GetXaxis().SetLabelSize(0.04)

if(is_center):
    hist_kinbin[0].GetXaxis().CenterTitle()
    hist_kinbin[0].GetYaxis().CenterTitle()

if doResponseMatrix==False and doLOreco==False:
    tdr.cmsPrel(-1,13.)
if doResponseMatrix==True or doLOreco==True:
    if(year=='2016'):
        tdr.cmsPrel(35900., 13.,simOnly=True,thisIsPrelim=True)
    elif(year=='2017'):
        tdr.cmsPrel(41500., 13.,simOnly=True,thisIsPrelim=True)


latex = TLatex()
latex.SetTextSize(0.65*gStyle.GetPadTopMargin())
latex.SetNDC()

if (wilson=="cLXX"): modwilson = "c_{L,XX}=-c_{L,YY}=" + str(cmunu)
if (wilson=="cLXY"): modwilson = "c_{L,XY}=c_{L,YX}=" + str(cmunu)
if (wilson=="cLXZ"): modwilson = "c_{L,XZ}=c_{L,ZX}=" + str(cmunu)
if (wilson=="cLYZ"): modwilson = "c_{L,YZ}=c_{L,ZY}=" + str(cmunu)

if (wilson=="cRXX"): modwilson = "c_{R,XX}=-c_{R,YY}=" + str(cmunu)
if (wilson=="cRXY"): modwilson = "c_{R,XY}=c_{R,YX}=" + str(cmunu)
if (wilson=="cRXZ"): modwilson = "c_{R,XZ}=c_{R,ZX}=" + str(cmunu)
if (wilson=="cRYZ"): modwilson = "c_{R,YZ}=c_{R,ZY}=" + str(cmunu)

if (wilson=="cXX"): modwilson = "c_{XX}=-c_{YY}=" + str(cmunu)
if (wilson=="cXY"): modwilson = "c_{XY}=c_{YX}=" + str(cmunu)
if (wilson=="cXZ"): modwilson = "c_{XZ}=c_{ZX}=" + str(cmunu)
if (wilson=="cYZ"): modwilson = "c_{YZ}=c_{ZY}=" + str(cmunu)

if (wilson=="dXX"): modwilson = "d_{XX}=-d_{YY}=" + str(cmunu)
if (wilson=="dXY"): modwilson = "d_{XY}=d_{YX}=" + str(cmunu)
if (wilson=="dXZ"): modwilson = "d_{XZ}=d_{ZX}=" + str(cmunu)
if (wilson=="dYZ"): modwilson = "d_{YZ}=d_{ZY}=" + str(cmunu)

latex.DrawLatex(0.25,0.9,modwilson)

if doLOreco:
    latex.DrawLatex(0.25,0.85,"LO reco sample "+year)
if doLOparticle and doResponseMatrix==False:
    latex.DrawLatex(0.25,0.85,"LO particle sample "+year)
if doLOparticle and doResponseMatrix==True:
    latex.DrawLatex(0.25,0.85,"LO particle sample "+year+" + response matrix")


#if(year=='2016'):
#    tdr.cmsPrel(35900.,13.)
#elif(year=='2017'):
#    tdr.cmsPrel(41530.,13.)

################################################################################
## Save
################################################################################

if doResponseMatrix==False:
    resultname = './results/'+year+'/other/'+observable+'_sme_time_'+wilson+suffix
if doResponseMatrix and doLOreco==False and doLOparticle==False:
    resultname = './results/'+year+'/other/'+observable+'_sme_time_'+wilson+suffix+'_reco_'+year
if doLOreco:
    resultname = './results/'+year+'/other/'+observable+'_sme_time_'+wilson+suffix+'_directreco_'+year 
if doLOparticle:
    resultname = './results/'+year+'/other/'+observable+'_sme_time_'+wilson+suffix+'_miniAODparticle_'+year
if doResponseMatrix and doLOreco==False and doLOparticle==True:
    resultname = './results/'+year+'/other/'+observable+'_sme_time_'+wilson+suffix+'_miniAODparticleResponseMatrix_'+year
if doResponseMatrix and doResponseMatrixLO and doLOreco==False and doLOparticle==True:
    resultname = './results/'+year+'/other/'+observable+'_sme_time_'+wilson+suffix+'_miniAODparticleResponseMatrixLO_'+year


#rootfile_output = TFile(resultname+'.root', "RECREATE")
#canvas.Write()
#canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
#rootfile_output.Close()


#raw_input('exit')
