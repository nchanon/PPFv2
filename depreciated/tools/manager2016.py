
def succeeds_job(percent):
    return 100./(percent)

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


################################################################################
# Trigger and systematics name list 
################################################################################


elecmu_trig_2016 = [
    'trg_muon_electron_mu23ele12_fired',
    'trg_muon_electron_mu23ele12DZ_fired',
    #'trg_muon_electron_mu12ele23_fired',
    #'trg_muon_electron_mu12ele23DZ_fired',
    'trg_muon_electron_mu8ele23_fired',
    'trg_muon_electron_mu8ele23DZ_fired'
]

mu_trig_2016 = [
    #'trg_muon_mu22eta21_fired',
    #'trg_muon_mutk22eta21_fired',
    'trg_muon_mu24_fired',
    'trg_muon_mutk24_fired'
]

elec_trig_2016 = [
    'trg_electron_ele27_fired',
    #'trg_electron_ele25eta21_fired',
    #'trg_electron_ele32eta21_fired'
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

