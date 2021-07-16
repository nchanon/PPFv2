#!/usr/bin/env python

import sys, os
##########

def succeeds_job(percent):
    return 100./(percent)

def inputs_name(path, year, nature):
    foo = os.listdir('./'+path+'/'+year+'/'+nature+'/')
    foo.sort()
    if(not foo):
        print 'Error : dataset list is empty'
    return foo


################################################################################
## 2016 
################################################################################

cross_section_2016 = {
    'MC_dibosons_WW'             :118.7,                # dibosons WW          
    'MC_dibosons_WW2'            :118.7,                # dibosons WW2         
    'MC_dibosons_WZ'             :47.13,                # dibosons WZ          
    'MC_dibosons_WZ2'            :47.13,                # dibosons WZ2         
    'MC_dibosons_ZZ'             :16.523,               # dibosons ZZ          
    'MC_dibosons_ZZ2'            :16.523,               # dibosons ZZ2         
    'MC_signal_dilep'            :89.048226,            # signal dilep       
    'MC_signal_hadronic'         :377.96006,            # signal hadronic    
    'MC_signal_semilep'          :366.91429,            # signal semilep     
    'MC_singletop_ST_s'          :10.32,                # singletop STs
    'MC_singletop_ST_tW_antitop' :35.85,                # singletop STt antitop
    'MC_singletop_ST_tW_top'     :35.85,                # singletop STt top
    'MC_singletop_ST_t_antitop'  :80.95,                # singletop tW  antitop
    'MC_singletop_ST_t_top'      :36.02,                # singletop tW  top
    'MC_ttx_TTW'                 :0.2043,               # TTX TTW         
    'MC_ttx_TTW2'                :0.2043,               # TTX TTW2        
    'MC_ttx_TTW3'                :0.4062,               # TTX TTW3        
    'MC_ttx_TTZ'                 :0.2529,               # TTX TTZ         
    'MC_ttx_TTZ2'                :0.2529,               # TTX TTZ2        
    'MC_ttx_TTZ3'                :0.2529,               # TTX TTZ3        
    'MC_ttx_TTZ4'                :0.5297,               # TTX TTZ4        
     #3.697,                # TTX  TTG        
     #3.697,                # TTX  TTG2        
    'MC_wjets_WJets'             :61526.7,              # wjets WJets       
    'MC_wjets_WJets2'            :61526.7,              # wjets WJets2      
    'MC_zjets_DY_50'             :6225.4,               # zjets DY        
    'MC_zjets_DY_1050'           :22635.14,              # zjets DY 10-50   
    'MC_zjets_DY_10502'          :22635.14,              # zjets DY 10-50   
    'MC_zjets_DY_10503'          :22635.14              # zjets DY 10-50       
}


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
    35.85,                # singletop tW  antitop
    35.85,                # singletop tW  top
    80.95,                # singletop STt antitop
    136.02,               # singletop STt top

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
    22635.14,              # zjets DY 10-50   
    22635.14,              # zjets DY 10-50   
    22635.14              # zjets DY 10-50   
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

systematic_rate_2016 = [ # AN2019-228 v5, TOP-20-006
    1.05,
    1.2, #ttx
    1.3, #singletop
    1.3, #dibosons
    1.3, #wjets
    1.2  #zjets
]

alt_list_2016 = [
'alt_MC_CP5Down',              
'alt_MC_CP5Down_pmx',          
'alt_MC_CP5Up',          
'alt_MC_CP5Up_pmx',      
'alt_MC_erdOn',
'alt_MC_erdOn_ext',
'alt_MC_GluonMove',
'alt_MC_hdampDown',
'alt_MC_hdampDown_pmx',
'alt_MC_hdampUp',
'alt_MC_hdampUp_pmx',
'alt_MC_mtop169p5',
'alt_MC_mtop175p5',
'alt_MC_QCDbased'
]

################################################################################
## 2017
################################################################################

cross_section_2017 = {
    'MC_dibosons_WW'              :118.7,    
    'MC_dibosons_WZ'              :47.13,    
    'MC_dibosons_ZZ'              :16.523,   
    'MC_signal_dilep'             :89.05,    
    'MC_signal_hadronic'          :377.96006,   
    'MC_signal_semilep'           :366.91429,   
    'MC_singletop_ST_s'           :10.32,    
    'MC_singletop_ST_s2'          :10.32,    
    'MC_singletop_ST_tW_antitop'  :35.5,     
    'MC_singletop_ST_tW_antitop2' :35.5,     
    'MC_singletop_ST_tW_top'      :35.5,     
    'MC_singletop_ST_tW_top2'     :35.5,     
    'MC_singletop_ST_t_antitop'   :80.95,    
    'MC_singletop_ST_t_top'       :136.02,   
    'MC_ttx_TTW'                  :0.2043,   
    'MC_ttx_TTW2'                 :0.2043,   
    'MC_ttx_TTW3'                 :0.4062,   
    'MC_ttx_TTZ'                  :0.2529,   
    'MC_ttx_TTZ2'                 :0.2529,   
    'MC_ttx_TTZ3'                 :0.5297,   
    'MC_wjets_WJets'              :61526.7,   
    'MC_wjets_WJets2'             :61526.7,   
    'MC_zjets_DY_1050'            :22635.1,  
    'MC_zjets_DY_50'              :6225.4,   
    #'MC_zjets_DY_502'             :6225.4    

}

cross_sec_2017 = [
    118.7,    # dibosons WW
    47.13,    # dibosons WZ
    16.523,   # dibosons ZZ 

    89.05,    # signal dilep
    377.96006,   # signal hadronic
    366.91429,   # signal semilep

    10.32,    # singletop STs
    10.32,    # singletop STs2
    35.85,     # singletop tW antitop
    35.85,     # singletop tw2 antitop
    35.85,     # singletop tW  top
    35.85,     # singletop tW2 top
    80.95,    # singletop STt antitop
    136.02,   # singletop STt top

    0.2043,   # TTX TTW         
    0.2043,   # TTX TTW2        
    0.4062,   # TTX TTW3        
    0.2529,   # TTX TTZ         
    0.2529,   # TTX TTZ2        
    0.5297,   # TTX TTZ3        

    61526.7,   # wjets WJets       
    61526.7,   # wjets Wjets2      

    22635.1,  # zjets DY 10-50     
    6225.4,   # zjets DY        
    #6225.4    # zjets DY2         
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

systematic_rate_2017 = [ # AN2019-228 v5, TOP-20-006
    1.05,
    1.2, #ttx
    1.3, #singletop
    1.3, #dibosons
    1.3, #wjets
    1.2  #zjets
]

alt_list_2017 = [
    'alt_MC_CP5Down',        
    'alt_MC_CP5Down_pmx',    
    'alt_MC_CP5Up',          
    'alt_MC_CP5Up_pmx',      
    'alt_MC_erdOn',        
    'alt_MC_GluonMove',    
    'alt_MC_hdampDown',
    'alt_MC_hdampDown_pmx',
    'alt_MC_hdampUp',
    'alt_MC_hdampUp_pmx',
    'alt_MC_mtop169p5',
    'alt_MC_mtop169p5_pmx',
    'alt_MC_mtop175p5',
    'alt_MC_mtop175p5_pmx',
    'alt_MC_QCDbased',
    'alt_MC_QCDbased_ext'
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


alt_list = {
    '2016' : alt_list_2016,
    '2017' : alt_list_2017,
    '2018' : 35.9
}

data_name = {
    '2016' : ['Run2016'],
    '2017' : ['Run2017'],
    '2018' : 0
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


jec_list = [ 
    'TotalUp',
    'TotalDown',
]


systematic_list = [
    'syst_elec_reco',
    'syst_elec_id',
    'syst_muon_id',
    'syst_muon_iso',
    'syst_pu',
    'syst_b',
    #
    # 'syst_c',
    #'syst_l',
    'syst_pt_top'
]

alt_syst_list = [
    'CP5Up',
    'CP5Down',
    'hdampUp',
    'hdampDown',
    'mtop169',
    'mtop175',
    'erdOn',
    'QCD',
    'GluonMove'
]

systematic_time_list = [
    'lumi_flat',
    'lumi_stability',
    'lumi_linearity',
    'emu_trig'
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

sample_MC = {
    '2016' : inputs_name('inputs', '2016', 'MC'),
    '2017' : inputs_name('inputs', '2017', 'MC'),
#    '2018' : trig_2018    
}

sample_ALT = {
    '2016' : inputs_name('inputs', '2016', 'ALT'),
    '2017' : inputs_name('inputs', '2017', 'ALT'),
#    '2018' : trig_2018    
}

sample_DATA = {
    '2016' : inputs_name('inputs', '2016', 'DATA'),
    '2017' : inputs_name('inputs', '2017', 'DATA'),
#    '2018' : trig_2018    
}