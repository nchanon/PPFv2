import os 

observable = [
    'm_dilep',
    'n_bjets',
    'met',
    'pt_lead',
    'pt_sublead',
    'n_jets',
]

year = [
    '2016',
    '2017'
]

for y in year:
    for o in observable:
        cmd = './bin/histograms_creator '+o+' '+y
        os.system(cmd)