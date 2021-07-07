import os
import argparse

###################
## Initialisation
###################

observables = [
    'm_dilep',
    #'n_bjets',
#    'pt_lead'
]

years = [
    '2016',
#    '2017'
]

for o in observables:
    for y in years:
        os.system('bash scripts/launch_inclusive_stuff.sh '+o+' '+y)
        os.system('bash scripts/launch_one_bin_stuff.sh '+o+' '+y)
        os.system('bash scripts/launch_unrolled_stuff.sh '+o+' '+y)

        os.system('python scripts/launch_sme_all.py '+o+' '+y)
