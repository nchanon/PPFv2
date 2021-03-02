#!/usr/bin/env python

import os
import sys
sys.path.append('./')


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('directory', help='export directory')
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
directory = args.directory
year = args.year

path_out = './combine/'+year+'/'+directory+'/results/'
path_in  = 'acarle@lyoserv.in2p3.fr:/gridgroup/cms/acarle/CMSSW_8_1_0/src/Analyse/'+directory+'/results/'

os.system('scp '+path_in+'*'+year+'* '+path_out)