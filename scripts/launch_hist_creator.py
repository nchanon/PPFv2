import os
import argparse

###################
## Initialisation
###################

observables = [
    'm_dilep',
    'pt_elec',
    'eta_elec',
    'pt_muon',
    'eta_muon',
    'n_jets',
    'n_bjets',
    'pt_lead',
    'pt_sublead',
    'b1_pt',
    'j1_pt',
    'j1_eta',
#    'b1_eta'
##    'b2_pt',
    'j2_pt',
    'j2_eta',
##    'b2_eta',
    'pt_emu'
##    'dijet_m'
]

years = [
    '2016',
    '2017'
]

for y in years:
    for o in observables:
        os.system('./bin/histograms_creator '+o+' '+y+' forComp'+ ' -2' )
