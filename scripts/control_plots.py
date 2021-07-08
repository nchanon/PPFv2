import os 


#    ['n_jets', '\"jets multiplicity\"'],
#['met', '\"met\"'],
observable = [
    ['m_dilep', '\"dilepton mass\"'],
    ['pt_lead', '\"pt leading lepton\"'],
    ['pt_sublead', '\"pt subleading lepton\"'],
    ['n_bjets', '\"b-jets multiplicity\"'],
    ['pt_elec', '\"pt electron\"'],
    ['pt_muon', '\"pt muon\"'],
    ['j1_pt', '\"pt leading jet\"'],
    ['b1_pt', '\"pt leading b-jet\"'],
    ['eta_elec', '\"eta electron\"'],
    ['eta_muon', '\"eta muon\"'],
    ['j1_eta', '\"eta leading jet\"'],
    ['b1_eta', '\"eta leading b-jet\"']
]

year = [
    '2016',
    '2017'
]

for y in year:
    for o in observable:
        cmd = 'python ./bin/data_mc_comparaison.py '+o[0]+' '+y+' '+o[1]
        os.system(cmd)