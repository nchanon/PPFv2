import os
import argparse

###################
## Initialisation
###################

parser = argparse.ArgumentParser()
parser.add_argument('observable', help='display your observable')
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
observable = args.observable
year = args.year


wilson_list = [
    'fXX_L',
    'fXY_L',
    'fXZ_L',
    'fYZ_L',
    'fXX_R',
    'fXY_R',
    'fXZ_R',
    'fYZ_R',
    'fXX_C',
    'fXY_C',
    'fXZ_C',
    'fYZ_C',
    'fXX_D',
    'fXY_D',
    'fXZ_D',
    'fYZ_D'
]


################################################################################
## Code body
################################################################################

for l in wilson_list:
    print  "  1) rootfile creation for combine"
    os.system('python ./bin/combine_unrolled.py '+observable+' '+year+' '+l)
    print  "  2) datacard creation for combine"
    os.system('./bin/card_creator '+observable+' '+year+' SME '+l)

print "  3) export to lyoserv"
os.system('python ./scripts/export_combine.py sme '+year)