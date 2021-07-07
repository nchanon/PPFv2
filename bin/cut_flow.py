import sys
sys.path.append('./')

from tools.sample_manager import *

import numpy as np
import argparse

from ROOT import TFile

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
year = args.year



################################################################################

def results_path(year, directory, file):
    return './results/'+year+'/'+directory+'/'+file

def efficiency_input(year, nature, sample, filter):
    return './inputs/'+year+'/'+nature+'/'+sample+'/'+filter+'/efficiency.txt'


def percent(x, total):
    return round(100.* x/total, 2)


def total_list(list):
    foo = 0
    for l in list:
        foo += l
    return foo

def event_filtered(file, is_after = True):
    event = []
    f = open(file, 'r')
    for i in f.read().split():
        try:
            event.append(int(i.strip()))
        except:
            pass
    if(is_after):
        return event[1]
    else:
        return event[0]

def cut_flow_human_readable(filter, year):
    filter_h = []
    if(year == '2016'):
        filter_h.append(filter[0]
                       +filter[1]
                       +filter[2]) #signal
        filter_h.append(filter[3]
                       +filter[4]
                       +filter[5]
                       +filter[6]
                       +filter[7]
                       +filter[8]
                       +filter[9]) #TTX
        filter_h.append(#filter[11]
                       filter[10]
                       +filter[11]
                       +filter[12]
                       #+filter[15]
                       +filter[13]
                       #+filter[15]
                       +filter[14]) #ST
        filter_h.append(#filter[19]
                       filter[15]
                       +filter[16]
                       +filter[17]
                       +filter[18]
                       +filter[19]
                       +filter[20]) #diboson
        filter_h.append(filter[21]
                       +filter[22]) #Wjets

        #filter_h.append(filter[13]+filter[14]) #Zjets
    elif(year == '2017'):
        filter_h.append(filter[0]
                       +filter[1]
                       +filter[2]) #signal
        filter_h.append(filter[3]
                       +filter[4]
                       +filter[5]
                       +filter[6]
                       +filter[7]
                       +filter[8]
                       +filter[9]
                       +filter[10]) #TTX
        filter_h.append(#filter[11]
                       filter[11]
                       +filter[12]
                       +filter[13]
                       #+filter[15]
                       +filter[14]
                       #+filter[15]
                       +filter[15]) #ST
        filter_h.append(#filter[19]
                       filter[16]
                       +filter[17]) #diboson
        filter_h.append(filter[18]
                       +filter[19]) #Wjets
        filter_h.append(filter[20]
                       +filter[21]
                       +filter[22]
                       ) #Zjets
    return filter_h

def latex_table(name, filter1, filter2, filter3, filter4):
    foo = r'''
    \begin{tabular}{|c|c|c|c|c|}\hline
     channel & no cuts & 2 leptons & 2 jets & 1 b-jets \\ \hline \hline
    '''
    for i in range(len(name)):
        foo += '    '+name[i].replace('_',' ')+' & '\
            +str(round(filter1[i]))+'$\pm$'+str(round(np.sqrt(filter1[i])))+\
            ' ['+str(percent(filter1[i],total_list(filter1)))+'\%] & '\
            +str(round(filter2[i]))+'$\pm$'+str(round(np.sqrt(filter2[i])))+\
            ' ['+str(percent(filter2[i],total_list(filter2)))+'\%] & '\
            +str(round(filter3[i]))+'$\pm$'+str(round(np.sqrt(filter3[i])))+\
            ' ['+str(percent(filter3[i],total_list(filter3)))+'\%] &'\
            +str(round(filter4[i]))+'$\pm$'+str(round(np.sqrt(filter4[i])))+\
            ' ['+str(percent(filter4[i],total_list(filter4)))+'\%]  '\
            +'  \\\\ \n    '
        if(i==0):
            foo += '\hline \n    '
    foo += '\hline \n    '
    foo += 'total MC & '+str(round(total_list(filter1)))\
            +'$\pm$'+str(round(np.sqrt(total_list(filter1))))\
            +' & '+str(round(total_list(filter2)))\
            +'$\pm$'+str(round(np.sqrt(total_list(filter2))))\
            +' & '+str(round(total_list(filter3))) \
            +'$\pm$'+str(round(np.sqrt(total_list(filter3))))\
            +' & '+str(round(total_list(filter4))) \
            +'$\pm$'+str(round(np.sqrt(total_list(filter4))))\
            +'\\\\ \n '
    foo += 'total data & '\
            +' & '\
            +' & '\
            +' & XXXX \\\\ \n'\
            +'\hline'\
            +' \n '
    foo += r'''\end{tabular}'''
    return foo

################################################################################

onebjet = 'OneBJets'
onedilep = 'OneDilepton'
twojets = 'TwoJets'


full         = []
one_bjet     = []
one_dilepton = []
two_jets     = []


for i in sample_MC[year]:
    full.append(event_filtered(efficiency_input(year,'MC', i, onedilep), False))
    one_dilepton.append(event_filtered(efficiency_input(year,'MC', i, onedilep), True))
    two_jets.append(event_filtered(efficiency_input(year,'MC', i, twojets), True))
    one_bjet.append(event_filtered(efficiency_input(year,'MC', i, onebjet), True))


integrals = []

rootfile = []
for i in range(len(sample_MC[year])):
    namefile = 'm_dilep'+'.root'
    print results_path(year,'flattree',namefile)
    rootfile.append(TFile(results_path(year,'flattree',namefile)))
    foo = rootfile[i].Get('m_dilep')
    integrals.append(foo.Integral())
del rootfile


name_channel_h = [
    '$t\\bar{t}$ signal',
    'TTX',
    'single top',
    'dibosons',
    'W$+$Jets',
    'Z$+$Jets'
]

integrals_h    = cut_flow_human_readable(integrals, year)
full_h         = cut_flow_human_readable(full, year)
one_bjet_h     = cut_flow_human_readable(one_bjet, year)
one_dilepton_h = cut_flow_human_readable(one_dilepton, year)
two_jets_h     = cut_flow_human_readable(two_jets, year)

to_one_dilep = []
to_two_jets = []
to_one_bjets = []

for i in range(len(full_h)):
    print i,one_bjet_h[i]
    to_one_bjets.append(two_jets_h[i]/float(one_bjet_h[i]))
    to_two_jets.append(one_dilepton_h[i]/float(one_bjet_h[i]))
    to_one_dilep.append(full_h[i]/float(one_bjet_h[i]))

for i in range(len(integrals_h)):
    print name_channel_h[i]+' full > '+str(to_one_dilep[i]*integrals_h[i])+" +- "+str(np.sqrt(to_one_dilep[i]*integrals_h[i]))
    print name_channel_h[i]+' one dilep > '+str(to_two_jets[i]*integrals_h[i])+" +- "+str(np.sqrt(to_two_jets[i]*integrals_h[i]))
    print name_channel_h[i]+' two jets > '+str(to_one_bjets[i]*integrals_h[i])+" +- "+str(np.sqrt(to_one_bjets[i]*integrals_h[i]))
    print name_channel_h[i]+' one bjets > '+str(integrals_h[i])+" +- "+str(np.sqrt(integrals_h[i]))
    print '---'

latex_full = []
latex_one_dilep = []
latex_two_jets = []
latex_one_bjet = []
for i in range(len(integrals_h)):
    latex_full.append(to_one_dilep[i]*integrals_h[i])
    latex_one_dilep.append(to_two_jets[i]*integrals_h[i])
    latex_two_jets.append(to_one_bjets[i]*integrals_h[i])
    latex_one_bjet.append(integrals_h[i])

## Latex part

header = r'''\documentclass{standalone}
\begin{document}
'''
main =  latex_table(name_channel_h, latex_full, latex_one_dilep, latex_two_jets, latex_one_bjet)
footer = r'''
\end{document}
''' 

content = header+main+footer

with open('./results/cutflow/cutflow'+year+'.tex','w') as f:
    f.write(content)

