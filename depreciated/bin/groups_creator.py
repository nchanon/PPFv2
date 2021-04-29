import time
t_start = time.clock()

import sys
sys.path.append('./')
from tools.generator_manager import *

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')
parser.add_argument('nature', help='MC, SYST')

args = parser.parse_args()
year = args.year
nature = args.nature

################################################################################
# Variables
################################################################################

################################################################################
# Generate pure TH1 without fancy style 
################################################################################

if nature == 'MC':
    generate_MC_groups(year, variables, ttbar_list)

elif nature == 'SYST':
    for syst in systematic_list :
        generate_SYST_groups(year, variables, syst, ttbar_list)


t_end = time.clock()

print "elapsed time : ", (t_end - t_start)/60., "min"