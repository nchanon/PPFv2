import sys
sys.path.append('./')

from tools.style_manager import *
from tools.sample_manager import *

import numpy as np
import argparse

from ROOT import TFile, TH1F, TH1

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
year = args.year

sample = 'MC_signal_dilep'
observable = 'm_dilep'

################################################################################

def form(key, dictionnary):
    return ('%12s %10s %s5' % (key, dictionnary[key], int(np.sqrt(dictionnary[key]))))

def table(integrals):
    for key in integrals:
        if key != 'total' and key != 'data_obs':
            print  form(key, integrals)
    print ''
    print form('total', integrals)
    print form('data_obs', integrals)

def path(year, sample, isInput=True):
    if isInput:
        return './inputs/'+year+'/MC/'+sample+'/'
    else:
        return './inputs/'+year+'/'+sample+'/'

def event_filtered(file, is_final = True):
    event = []
    f = open(file, 'r')
    for i in f.read().split():
        try:
            event.append(int(i.strip()))
        except:
            pass
    if is_final:
        return event[1]
    else:
        return event[0]

def get_events(year, sample, filter, is_final=True):
    p = path(year, sample)+filter+'/efficiency.txt'
    return event_filtered(p,is_final)

def cut_flow(year, sample, weight=1.):
    #print 'Get cut-flow ',year,' ',sample
    return {
        'full'     : weight*get_events(year, sample, 'OneDilepton', False),
        'onedilep' : weight*get_events(year, sample, 'OneDilepton'),
        'twojets'  : weight*get_events(year, sample, 'TwoJets'),
        'onebjet'  : weight*get_events(year, sample, 'OneBJets')
    }

def filtre(samples):
    foo = {}
#    for l in ttbar_list:
#        foo.update({l : })
    print samples
    for l in samples.keys():
	print l
        for g in ttbar_list:
            if l.find(g) != -1:
                foo[g] += int(sample[l])
            if l.find('jets') != -1:
                foo[g] += int(sample[l])
    return foo

'''
ratio = {
    'full'     : 1.,
    'onedilep' : get_ratio(year, sample, 'OneDilepton'),
    'twojets'  : get_ratio(year, sample, 'TwoJets'),
    'onebjet'  : get_ratio(year, sample, 'OneBJets')
}
'''

effective_N0 = generate_eventN0(year, sample_MC[year])
mc_rescale   = rescaling(year, effective_N0)
print 'N0=',effective_N0
print 'mc_rescale=',mc_rescale

#values = {}
#for l in range(len(sample_MC[year])):
    #values.append(
#    values.update(
#        { sample_MC[year][l] : cut_flow(year, sample_MC[year][l], mc_rescale[l])}
#        { sample_MC[year][l] : cut_flow(year, sample_MC[year][l], 77*mc_rescale[l])}
#    )

print ttbar_list

cutflowtable = []

for i in range(len(ttbar_list)): 
  proc = [ttbar_list[i], 0, 0, 0, 0]
  cutflowtable.append(proc)
  #cutflowtable[i].append(0)
  #cutflowtable[i].append(0)
  #cutflowtable[i].append(0)
  #cutflowtable[i].append(0) 
print cutflowtable

 
for l in range(len(sample_MC[year])):
	print sample_MC[year][l]
	print 'nocut'     , mc_rescale[l]*get_events(year, sample_MC[year][l], 'OneDilepton', False)
	print 'onedilep' , mc_rescale[l]*get_events(year, sample_MC[year][l], 'OneDilepton')
	print 'twojets'  , mc_rescale[l]*get_events(year, sample_MC[year][l], 'TwoJets')
	print 'onebjet'  , mc_rescale[l]*get_events(year, sample_MC[year][l], 'OneBJets')
	n_nocut = mc_rescale[l]*get_events(year, sample_MC[year][l], 'OneDilepton', False)
	n_onedilep = mc_rescale[l]*get_events(year, sample_MC[year][l], 'OneDilepton')
	n_twojets = mc_rescale[l]*get_events(year, sample_MC[year][l], 'TwoJets')
	n_onebjet = mc_rescale[l]*get_events(year, sample_MC[year][l], 'OneBJets')
	for i in range(len(ttbar_list)):
	  if (sample_MC[year][l].find(ttbar_list[i])!=-1): 
		cutflowtable[i][1] += int(n_nocut)
                cutflowtable[i][2] += int(n_onedilep)
                cutflowtable[i][3] += int(n_twojets)
                cutflowtable[i][4] += int(n_onebjet)
	  if (sample_MC[year][l].find("wjets")!=-1 or sample_MC[year][l].find("zjets")!=-1):
		cutflowtable[-1][1] += int(n_nocut)
                cutflowtable[-1][2] += int(n_onedilep)
                cutflowtable[-1][3] += int(n_twojets)
                cutflowtable[-1][4] += int(n_onebjet)
		

print 'Process		nocut		onedilep		twojets		onebjet'
for i in [0, 1, 2, 3, 4]:
#print ttbar_list[0], '             ',ttbar_list[1], '             ',ttbar_list[2], '             ',ttbar_list[3], '             ',ttbar_list[4]
#print cutflowtable[0][1], '             ',cutflowtable[1][1],'             ',cutflowtable[2][1],'             ',cutflowtable[3][1],'             ',cutflowtable[4][1]
  print ttbar_list[i], '		', cutflowtable[i][1], '	',cutflowtable[i][2],'		',cutflowtable[i][3],'		',cutflowtable[i][4]
#integrals_cut = filtre(values)

#print integrals_cut

integrals = {}

nbin = 0
min_bin = 0
max_bin = 0
TH1.SetDefaultSumw2(1)

datafile_input = TFile('./results/'+year+'/flattree/'+observable+'_data.root')
rootfile_input = TFile('./results/'+year+'/flattree/'+observable+'_forComp.root')

integrals.update({'data_obs' : int(datafile_input.Get('data_obs').Integral())})
integrals.update({'signal' : int(rootfile_input.Get('signal').Integral())})
total_mc = integrals['signal']

background_integral = 0
for l in rootfile_input.GetListOfKeys():
    for s in ttbar_list:
        if(l.GetName() == s and not l.GetName() == 'signal'):
            integrals.update({l.GetName() : int(rootfile_input.Get(l.GetName()).Integral()) })
            total_mc += integrals[l.GetName()]
            #ratio.update({l.GetName() : cut_flow})

integrals.update({'total' : total_mc})

table(integrals) 

'''
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


with open('./results/cutflow/cutflow'+year+'.tex','w') as f:
    f.write(content)

'''
