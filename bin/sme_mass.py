#!/usr/bin/env python

import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import argparse

from ROOT import TFile, TH1, TCanvas, TH1F, THStack, TString
from ROOT import TLegend, TApplication, TRatioPlot, TPad, TFrame
from ROOT import TStyle, gStyle, TColor, TLatex

import tools.tdrstyle as tdr
tdr.setTDRStyle()

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')
parser.add_argument('title', help='display your observable title')
parser.add_argument('wilson',nargs='?', help='wilson coefficent (cLXX, cRXX, ...)', default='cLXX')
parser.add_argument('timebin', help='display the time bin')


args = parser.parse_args()
observable = args.observable
year = args.year
title = args.title
wilson = args.wilson
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


doResponseMatrix = True
cmunu = 0.01


def integral_complete(histo, max_bin): 
     return histo.Integral()
#    return histo.Integral()+histo.GetBinContent(int(max_bin+1))+histo.GetBinContent(0)

################################################################################
## Create Canvas 
################################################################################

nbin = 0
min_bin = 0
max_bin = 0
legend_coordinates = observable_values(observable)[1]
TH1.SetDefaultSumw2(1)
background_integral_i = []
background_integral = 0
data_integral = 0
syst_up_integral = 0
syst_down_integral = 0
canvas = TCanvas('stack_'+observable,'stack_'+observable, 800, 800)
canvas.UseCurrentStyle()

r = 0.3
epsilon = 0.1

pad1 = TPad("pad1", "pad1", 0, r-epsilon, 1, 1)
pad1.SetBottomMargin(epsilon)
canvas.cd()
#pad1.SetLogy()
pad1.Draw()
pad1.cd()


################################################################################
## Get particle-level time modulations 
################################################################################

rootfile_input = TFile('./results/'+year+'/flattree/'+observable+stimebin+'.root')
sme_file = TFile('./results/'+year+'/flattree/'+observable+'_sme.root')


# signal observable
hist_signal = rootfile_input.Get('signal')
hist_signal.Scale(1./24)

nbin    = hist_signal.GetNbinsX()
min_bin = hist_signal.GetXaxis().GetXmin()
max_bin = hist_signal.GetXaxis().GetXmax()

# sme time
sme_inclusive = sme_file.Get(wilson)
sme_inclusive.Scale(cmunu)
#binmax_sme_inclusive = sme_inclusive.GetMaximumBin()

sme_time_kinbin = []
if observable=="m_dilep":
    for i in range(8):
        if (i>0):
            sme_time_kinbin.append(sme_file.Get(wilson+"_"+str(i)))
            sme_time_kinbin[-1].Scale(cmunu)
if observable=="pt_emu":
    for i in range(7):
        sme_time_kinbin.append(sme_file.Get(wilson+"_"+str(i)))
        sme_time_kinbin[-1].Scale(cmunu)
if observable=="n_bjets":
    for i in range(6):
        sme_time_kinbin.append(sme_file.Get(wilson+"_"+str(i)))
        sme_time_kinbin[-1].Scale(cmunu)


################################################################################
## Apply response matrix 
################################################################################

doNormalizeByColumn = True

sme_time_kinbin_rec = []

if doResponseMatrix:
    fResponseMatrix = TFile("results/"+year+"/flattree/"+observable+stimebin+".root")
    hResponseMatrix = fResponseMatrix.Get("signal_responseMatrix") #singletop+"_responseMatrix")

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
        sme_time_kinbin_rec.append(sme_time_kinbin[ix].Clone())
        for it in range(sme_time_kinbin[ix].GetNbinsX()):
            recbincontent = 0
            for iy in range(nbinY):
                if iy==0:
                    genbin = nbinY-1
                else:
                    genbin = iy-1
                #print('ix='+str(ix)+' iy='+str(iy)+' genbin='+str(genbin))
                recbincontent += hResponseMatrix.GetBinContent(1+ix,1+iy)*sme_time_kinbin[genbin].GetBinContent(1+it)
            sme_time_kinbin_rec[ix].SetBinContent(1+it, recbincontent)
        sme_time_kinbin[ix] = sme_time_kinbin_rec[ix]

if (observable=='n_bjets' and doResponseMatrix):
    nhist_max = 4

#######################
# Observable histograms
#######################


hist_sme_inclusive = hist_signal.Clone()
hist_sme_inclusive.SetName("sme_inclusive")
hist_sme_inclusive.SetTitle("sme_inclusive")
hist_sme_inclusive.Scale(1+sme_inclusive.GetMaximum()) #show SM kin scaled by maximum SME variation

min = hist_sme_inclusive.GetMinimum()
max = hist_sme_inclusive.GetMaximum()

hist_sme_kin = []
for i in range(24):
    hist_sme_kin.append(hist_signal.Clone())
    hist_sme_kin[-1].SetName("sme_"+str(i))
    hist_sme_kin[-1].SetTitle("sme_"+str(i))
    for m in range(nbinX):
        val_t = hist_signal.GetBinContent(m+1)
	ft_m = sme_time_kinbin[m].GetBinContent(i+1)
	hist_sme_kin[-1].SetBinContent(m+1, val_t*(1+ft_m))
	if (max < hist_sme_kin[-1].GetMaximum()): max = hist_sme_kin[-1].GetMaximum()

usedbins = [2, 5, 9, 12]

# background
#background_integral = 0
#hist_background  = TH1F("background", "background", nbin, min_bin, max_bin)
#for l in rootfile_input.GetListOfKeys():
#    for s in ttbar_list:
#        if(l.GetName() == s and not l.GetName() == 'signal'):
#            print l.GetName()
#            hist_background.Add(rootfile_input.Get(l.GetName()))
#            background_integral_i.append([integral_complete(rootfile_input.Get(l.GetName()), max_bin), l.GetName()])
#            background_integral += integral_complete(rootfile_input.Get(l.GetName()), max_bin)



################################################################################
## Legend stuff
################################################################################

skin=''
if observable=='m_dilep':
    skin='m_{ll}'
if observable=='pt_emu':
    skin='p_{T,ll}'
if observable=='n_bjets':
    skin='n_{bjets}'

legend_args = (0.57, 0.7, 0.95, 0.91, '', 'NDC')
legend = TLegend(*legend_args)
legend.AddEntry(hist_signal, "t#bar{t} SM any time bin", "l")
legend.SetTextSize(0.03)
#legend.AddEntry(hist_sme_inclusive, "SME "+skin+"-inclusive, max time bin", "l")
legend.AddEntry(hist_sme_kin[usedbins[0]], "SME "+skin+"-dependent, "+str(usedbins[0])+"-"+str(usedbins[0]+1)+"h", "l")
legend.AddEntry(hist_sme_kin[usedbins[1]], "SME "+skin+"-dependent, "+str(usedbins[1])+"-"+str(usedbins[1]+1)+"h", "l")
legend.AddEntry(hist_sme_kin[usedbins[2]], "SME "+skin+"-dependent, "+str(usedbins[2])+"-"+str(usedbins[2]+1)+"h", "l")
legend.AddEntry(hist_sme_kin[usedbins[3]], "SME "+skin+"-dependent, "+str(usedbins[3])+"-"+str(usedbins[3]+1)+"h", "l")


#legend.AddEntry(hist_background, "non-t#bar{t}", "f")
#legend.AddEntry(hist_data, "data")
#legend_box(legend, legend_coordinates)

################################################################################
## Draw stuff
################################################################################

#stack = THStack()
#stack.Add(hist_background)
#stack.Add(hist_signal)
#stack.Draw("E HIST")
#hist_data.Draw("E SAME")

gStyle.SetPalette(55)

hist_signal.GetXaxis().SetLabelSize(0);
hist_signal.GetYaxis().SetRangeUser(0, max*1.2)
hist_signal.SetLineWidth(2)
hist_signal.SetLineColor(1)
hist_signal.Draw("E HIST")

#hist_sme_inclusive.SetLineWidth(2)
#hist_sme_inclusive.SetLineColor(2)
#hist_sme_inclusive.Draw("histsame")

hist_sme_kin[usedbins[0]].SetLineWidth(2)
hist_sme_kin[usedbins[0]].SetLineColor(2)
hist_sme_kin[usedbins[0]].Draw("histsame")

hist_sme_kin[usedbins[1]].SetLineWidth(2)
hist_sme_kin[usedbins[1]].SetLineColor(9)
hist_sme_kin[usedbins[1]].Draw("histsame")

hist_sme_kin[usedbins[2]].SetLineWidth(2)
hist_sme_kin[usedbins[2]].SetLineColor(800)
hist_sme_kin[usedbins[2]].Draw("histsame")

hist_sme_kin[usedbins[3]].SetLineWidth(2)
hist_sme_kin[usedbins[3]].SetLineColor(823)
hist_sme_kin[usedbins[3]].Draw("histsame")

legend.Draw("SAME")


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


################################################################################
## Set Style
################################################################################

# line_color, line_width, fill_style, marker_style
#style_histo(hist_signal, 2, 1, 2, 3004, 0)
#style_histo(hist_background, 4, 1, 4, 3005, 0)
#style_histo(hist_data, 1, 1, 0, 3001, 1, 20)

#style_labels_counting(stack, 'Events', title)
#stack.GetXaxis().SetLabelSize(0)
#stack.GetXaxis().SetTitleSize(0)

if(year=='2016'):
    tdr.cmsPrel(35900., 13.,simOnly=True)
elif(year=='2017'):
   tdr.cmsPrel(41500., 13.,simOnly=True)

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
h_one.SetLineWidth(2)
h_one.SetLineColor(1)
#h_one.Draw("SAME")

h_denom = hist_signal
ratio = THStack()

#h_inc = hist_sme_inclusive.Clone()
#h_inc.Divide(h_denom)
#ratio.Add(h_inc)
h_0 = hist_sme_kin[usedbins[0]].Clone()
h_0.SetName("h_0")
h_0.Divide(h_denom)
ratio.Add(h_0)
h_1 = hist_sme_kin[usedbins[1]].Clone()
h_1.SetName("h_1")
h_1.Divide(h_denom)
ratio.Add(h_1)
h_2 = hist_sme_kin[usedbins[2]].Clone()
h_2.SetName("h_2")
h_2.Divide(h_denom)
ratio.Add(h_2)
h_3 = hist_sme_kin[usedbins[3]].Clone()
h_3.SetName("h_3")
h_3.Divide(h_denom)
ratio.Add(h_3)

ratio.SetMaximum(1+ratio_coef)
ratio.SetMinimum(1-ratio_coef)
ratio.Draw()
#h_one.Draw("SAME")

style_labels_counting(ratio, 'Ratio SME/SM', title)
ratio.GetYaxis().SetLabelSize(0.1)
ratio.GetYaxis().SetTitleSize(0.1)
ratio.GetYaxis().SetTitleOffset(0.5)

ratio.GetXaxis().SetLabelSize(0.15)
ratio.GetXaxis().SetTitleSize(0.17)
ratio.GetXaxis().SetLabelOffset(0.01)

ratio.Draw("HIST")
h_one.Draw("histSAME")
#h_inc.Draw("histSAME")
h_0.Draw("histSAME")
h_1.Draw("histSAME")
h_2.Draw("histSAME")
h_3.Draw("histSAME")

################################################################################
## Save
################################################################################

resultname = './results/'+year+'/other/sme_'+observable+'_'+wilson
if doResponseMatrix==True:
    resultname += '_reco_'+year

#rootfile_output = TFile(resultname+'.root', "RECREATE")
#canvas.Write()
#hist_signal.Write()
#hist_sme_inclusive.Write()
#hist_sme_kin[usedbins[0]].Write()
#hist_sme_kin[usedbins[1]].Write()
#hist_sme_kin[usedbins[2]].Write()
#hist_sme_kin[usedbins[3]].Write()
#canvas.SaveAs(resultname+'.png')
canvas.SaveAs(resultname+'.pdf')
#rootfile_output.Close()

#raw_input('exit')

#print 'Signal integral     : ','%.2f'%signal_integral
#for i in background_integral_i:
#    print i[1], '      : ', '%.2f'%i[0]
#print 'Total Background integral : ', '%.2f'%background_integral
#print 'Total MC integral : ', '%.2f'%(signal_integral+background_integral)
#print 'Data integral       : ', data_integral
#print 'Data/MC agreement  : ', '%.1f'%(100*(signal_integral+background_integral-data_integral)/(signal_integral+background_integral)), '%'

#raw_input('exit')
