import os
import argparse

###################
## Initialisation
###################

observables = [
    'pt_muon',
    'eta_muon',

    'pt_elec',
    'eta_elec',

    'm_dilep',
    'pt_emu'

    'j1_pt',
    'j1_eta',

    'j2_pt',
    'j2_eta',

    'b1_pt',
    'b1_eta',

    'n_jets',
    'n_bjets'

#    'm_dilep',
#    'pt_elec',
#    'eta_elec',
#    'pt_muon',
#    'eta_muon',
#    'n_jets',
#    'n_bjets',
#    'pt_lead',
#    'pt_sublead',
#    'b1_pt',
#    'j1_pt',
#    'j1_eta',
#    'b1_eta',
#    'j2_pt',
#    'j2_eta',
#    'pt_emu'

##    'b2_pt',
##    'b2_eta',
##    'dijet_m'
]

years = [
    '2016',
    '2017'
]

for y in years:
    for o in observables:
        os.system('./bin/histograms_creator '+o+' '+y+' forComp'+ ' -2' )
