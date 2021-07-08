import os
import argparse

###################
## Initialisation
###################

observables = [
    'm_dilep',
    'n_bjets',
    'pt_lead',
    'pt_sublead',
    'pt_elec',
    'pt_muon',
    'b1_pt',
    'j1_pt',
    'eta_elec',
    'eta_muon',
    'j1_eta',
    'b1_eta'
]

years = [
    '2016',
    '2017'
]

for y in years:
    for o in observables:
        os.system('./bin/histograms_creator '+o+' '+y)
