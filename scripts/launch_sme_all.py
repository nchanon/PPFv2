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
    'cLXX',
    'cLXY',
    'cLXZ',
    'cLYZ',
    'cRXX',
    'cRXY',
    'cRXZ',
    'cRYZ',
    'cXX',
    'cXY',
    'cXZ',
    'cYZ',
    'dXX',
    'dXY',
    'dXZ',
    'dYZ',
    'sme_all'
]


################################################################################
## Code body
################################################################################

for l in wilson_list:
    print  "  1) rootfile creation for combine"
    os.system('python ./bin/combine_unrolled.py '+observable+' '+year+' '+l)
    #print  "  2) datacard creation for combine"
    #os.system('./bin/card_creator '+observable+' '+year+' SME '+l)

print "  3) export to lyoserv"
os.system('python ./scripts/export_combine.py sme '+year)
