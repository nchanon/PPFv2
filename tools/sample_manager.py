
##########

def succeeds_job(percent):
    return 100./(percent)

################################################################################
## 2016 
################################################################################

cross_sec_2016 = [ # alphabetic order
    118.7,                # dibosons WW          
    118.7,                # dibosons WW2         
    47.13,                # dibosons WZ          
    47.13,                # dibosons WZ2         
    16.523,               # dibosons ZZ          
    16.523,               # dibosons ZZ2         

    89.048226,            # signal dilep       
    377.96006,            # signal hadronic    
    366.91429,            # signal semilep     

    10.32,                # singletop STs
    80.95,                # singletop STt antitop
    136.02,               # singletop STt top
    35.85,                # singletop tW  antitop
    35.85,                # singletop tW  top

    0.2043,               # TTX TTW         
    0.2043,               # TTX TTW2        
    0.4062,               # TTX TTW3        
    0.2529,               # TTX TTZ         
    0.2529,               # TTX TTZ2        
    0.2529,               # TTX TTZ3        
    0.5297,               # TTX TTZ4        
    #3.697,                # TTX  TTG        
    #3.697,                # TTX  TTG2        

    61526.7,              # wjets WJets       
    61526.7,              # wjets WJets2      

    6225.4,               # zjets DY        
    #22635.14              # zjets DY 10-50   
    #22635.14              # zjets DY 10-50   
    #22635.14              # zjets DY 10-50   
]

effective_data_event_2016 = [
    succeeds_job(100.),    #MuonEG_B
    succeeds_job(100.),    #MuonEG_C
    succeeds_job(100.),    #MuonEG_D
    succeeds_job(100.),    #MuonEG_E
    succeeds_job(100.),    #MuonEG_F
    succeeds_job(100.),    #MuonEG_G
    succeeds_job(100.),    #MuonEG_H
    succeeds_job(100.),    #SingleElectron_B
    succeeds_job(100.),    #SingleElectron_C
    succeeds_job(100.),    #SingleElectron_D
    succeeds_job(100.),    #SingleElectron_E
    succeeds_job(100.),    #SingleElectron_F
    succeeds_job(100.),    #SingleElectron_G
    succeeds_job(100.),    #SingleElectron_H
    succeeds_job(100.),    #SingleMuon_B
    succeeds_job(100.),    #SingleMuon_C
    succeeds_job(100.),    #SingleMuon_D
    succeeds_job(100.),    #SingleMuon_E
    succeeds_job(100.),    #SingleMuon_F
    succeeds_job(100.),    #SingleMuon_G
    succeeds_job(100.)     #SingleMuon_H
]

systematic_rate_2016 = [
    1.3,
    1.2,
    1.5,
    1.5,
    1.5
]

################################################################################
## 2017
################################################################################

cross_sec_2017 = [
    118.7,    # dibosons WW
    47.13,    # dibosons WZ
    16.523,   # dibosons ZZ 

    89.05,    # signal dilep
    380.11,   # signal hadronic
    364.31,   # signal semilep

    10.32,    # singletop STs
    10.32,    # singletop STs2
    35.5,     # singletop tW antitop
    35.5,     # singletop tw2 antitop
    35.5,     # singletop tW  top
    35.5,     # singletop tW2 top
    80.95,    # singletop STt antitop
    136.02,   # singletop STt top

    0.2043,   # TTX TTW         
    0.2043,   # TTX TTW2        
    0.4062,   # TTX TTW3        
    0.2529,   # TTX TTZ         
    0.2529,   # TTX TTZ2        
    0.5297,   # TTX TTZ3        

    0.4062,   # wjets WJets       
    0.4062,   # wjets Wjets2      

    22635.1,  # zjets DY 10-50     
    6225.4,   # zjets DY        
    6225.4    # zjets DY2         
]

effective_data_event_2017 = [
    succeeds_job(100.),    #MuonEG_B
    succeeds_job(100.),    #MuonEG_C
    succeeds_job(100.),    #MuonEG_D
    succeeds_job(95.8),    #MuonEG_E
    succeeds_job(100.),    #MuonEG_F
    succeeds_job(87.6),    #SingleElectron_B
    succeeds_job(100.),    #SingleElectron_C
    succeeds_job(100.),    #SingleElectron_D
    succeeds_job(83.3),    #SingleElectron_E
    succeeds_job(76.5),    #SingleElectron_F
    succeeds_job(100.),    #SingleMuon_B
    succeeds_job(100.),    #SingleMuon_C
    succeeds_job(100.),    #SingleMuon_D
    succeeds_job(100.),    #SingleMuon_E
    succeeds_job(100.)     #SingleMuon_F
]

systematic_rate_2017 = [
    1.3,
    1.2,
    1.5,
    1.5,
    1.5
]

################################################################################
## 2018 
################################################################################

#
#

################################################################################
################################################################################

effective_data_event = {
    '2016' : effective_data_event_2016,
    '2017' : effective_data_event_2017,
    '2018' : 35.9,
}

data_name = {
    '2016' : ['Run2016'],
    '2017' : ['Run2017'],
    '2018' : 0
}

luminosity = {
    '2016' : 35.9,
    '2017' : 41.53,
    '2018' : 0
}

cross_sec = {
    '2016' : cross_sec_2016,
    '2017' : cross_sec_2017,
    '2018' : 0
}


ttbar_list = [ 
    'signal',
    'ttx',
    'singletop',
    'dibosons',
    'wjets',
    'zjets'
]

systematic_list = [
    'syst_elec_reco',
    'syst_elec_id',
    'syst_muon_id',
    'syst_muon_iso',
    'syst_em_trig',
    'syst_pu'
]

systematic_time_list = [
    'lumi_inclusive',
    'lumi_stability',
    'lumi_linearity'
]

trig_2016 = [
    'trg_muon_electron_mu23ele12_fired',
    'trg_muon_electron_mu23ele12DZ_fired',
    'trg_muon_electron_mu8ele23_fired',
    'trg_muon_electron_mu8ele23DZ_fired',
    'trg_muon_mu24_fired',
    'trg_muon_mutk24_fired',
    'trg_electron_ele27_fired'
]

trig_2017 = [
    'trg_muon_electron_mu8ele23DZ_fired',
    'trg_muon_electron_mu23ele12_fired',
    'trg_muon_mu27_fired',
    'trg_electron_ele35_fired'
]

triggers = {
    '2016' : trig_2016,
    '2017' : trig_2017,
#    '2018' : trig_2018
}

systematic_rate = {
    '2016' : systematic_rate_2016,
    '2017' : systematic_rate_2017,
#    '2018' : trig_2018
}