import os 

observable = [
    ['m_dilep', '\"dilepton mass\"'],
    ['n_jets', '\"jets multiplicity\"'],
    ['met', '\"met\"'],
    ['pt_lead', '\"pt_leading lepton\"'],
    ['pt_sublead', '\"pt_subleading lepton\"'],
    ['n_bjets', '\"b-jets multiplicity\"']

]

year = [
    '2016',
    '2017'
]

for y in year:
    for o in observable:
        cmd = 'python ./bin/data_mc_comparaison.py '+o[0]+' '+y+' '+o[1]
        os.system(cmd)