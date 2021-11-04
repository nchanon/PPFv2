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

path_in = './combine/'+year+'/'+directory+'/inputs/'
path_out = '/gridgroup/cms/nchanon/CMSSW_10_2_13/src/combine-ttbar/'+directory+'/inputs/'+year+'/'

os.system('scp '+path_in+'* '+path_out)
